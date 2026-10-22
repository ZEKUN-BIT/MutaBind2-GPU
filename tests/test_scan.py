import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
from unittest.mock import patch

import pytest
from mutabind2_tool.cli import main
from mutabind2_tool.core.runtime import TOOL_ROOT, RESOURCES, load_config
from mutabind2_tool.core.structure import inspect_job, prepare_structure
from mutabind2_tool.scan import pdb, plan, provean, scheduler, task


def options(tmp_path, **changes):
    values=dict(output=tmp_path/'scan',pdb=TOOL_ROOT/'examples/1cse.pdb',partner1='E',partner2='I',
                mutations=TOOL_ROOT/'examples/mutations.tsv',sites=None,mutated_partner=None,
                fragment_mode='auto',fragments=None,gap_threshold=2.0,cap_mode='charged',
                surface_backend='local',provean_mode='fragment',provean_cache=None,max_per_batch=19,
                workers=2,foldx_workers=2,pb_executable=None,vdw_executable=None,backend='gpu',
                config=TOOL_ROOT/'config/tools.example.json',no_annotations=False,resume=False)
    values.update(changes)
    return argparse.Namespace(**values)


def make_cache(path, sequence='ACDE', chain='X', mutation='C2A', score=-1.25):
    path.mkdir(parents=True)
    (path/'query.fasta').write_text('>full\n'+sequence+'\n')
    (path/'scores.tsv').write_text(f'Chain\tMutation\tPROVEAN\n{chain}\t{mutation}\t{score}\n')
    state=dict(stage='complete',chain=chain,score_count=1,
               query_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
               query_file_sha256=provean.sha256(path/'query.fasta'),
               scores_sha256=provean.sha256(path/'scores.tsv'))
    (path/'status.json').write_text(json.dumps(state))
    return path


def test_general_plan_preserves_original_identity_and_detects_gaps(tmp_path):
    args=options(tmp_path)
    m,_=plan.build(args)
    assert m['mutation_identity']==[['1','I','L45P']]
    assert m['mutated_partner']==2 and m['backend']=='gpu'
    assert m['plan']['fragments']==2
    rows,_=inspect_job(args.output/'batches/batch_0001/batch.input','single')
    assert rows[0]['MutChain']=='B_1'
    assert not (args.output/'batches/batch_0001/fragments.pdb').is_symlink()
    out=tmp_path/'prepared'
    prepare_structure(args.output/'plan/mutations.input',out,'single')
    cleaned=scheduler.read_table(args.output/'plan/mutations.input.cleaned')
    assert cleaned[0]['Mutation_cleaned']=='LB38P'


def test_range_expansion_and_partner_one(tmp_path):
    args=options(tmp_path,mutations=None,sites='I:45-46,I:48',partner1='I',partner2='E')
    m,_=plan.build(args)
    assert m['plan']['mutations']==57 and m['plan']['sites']==3
    assert m['mutated_partner']==1
    assert {r[2][1:-1] for r in m['mutation_identity']}=={'45','46','48'}


@pytest.mark.parametrize('value',['A:9-2',';touch /tmp/unsafe','A:1,$(id)','A:',''])
def test_unsafe_site_selection_rejected(value):
    with pytest.raises(ValueError):pdb.parse_sites(value)


@pytest.mark.parametrize('changed',['plan/mutations.input','plan/mutations.tsv','plan/fragments.pdb','batches.json'])
def test_changed_plan_rejected_without_rewriting(tmp_path,changed):
    a=options(tmp_path);plan.build(a)
    p=a.output/changed;p.write_text(p.read_text()+'\n')
    previous=p.read_bytes()
    a.resume=True
    with pytest.raises(ValueError):plan.build(a)
    assert p.read_bytes()==previous


def test_configuration_relative_paths_resolve_against_config(tmp_path):
    p=tmp_path/'config.json';p.write_text(json.dumps({'foldx':'./missing','blast_db':'./db/nr'}))
    c=load_config(p)
    assert c['foldx']==str(tmp_path/'missing') and c['blast_db']==str(tmp_path/'db/nr')


def test_duplicate_mutations_rejected():
    s=pdb.parse_pdb(TOOL_ROOT/'examples/1cse.pdb')
    row={'chain':'I','mutation':'L45P'}
    with pytest.raises(ValueError,match='Duplicate'):pdb.validate_mutations(s,[row,row])


def test_no_fragment_mode_rejects_gapped_chains():
    s=pdb.parse_pdb(TOOL_ROOT/'examples/1cse.pdb')
    assert len(pdb.build_fragments(s,['E'],['I'],mode='none')[0]) == 2
    # Move the peptide N far away to make a genuine coordinate break.
    residue=s.residue('E','57');line=residue.atoms['N']
    residue.atoms['N']=line[:30]+f'{999.0:8.3f}'+line[38:]
    with pytest.raises(ValueError,match='breaks'):pdb.build_fragments(s,['E'],['I'],mode='none')


def test_terminal_flags_do_not_mark_every_natural_end(tmp_path):
    s=pdb.parse_pdb(TOOL_ROOT/'examples/1cse.pdb')
    f,_,_,_=pdb.build_fragments(s,['E'],['I'])
    a=pdb.annotate_sites(s,f,2)
    assert a[('I','8')]['is_truncated_terminus']==1
    assert a[('I','70')]['is_truncated_terminus']==0


def test_cached_conservation_uses_original_identity_but_engine_local_number(tmp_path):
    cache=make_cache(tmp_path/'cache')
    rows=[dict(Result_Id='17',Mutation_PDB='C99A',Mutation_cleaned='CB1A',MutChain='Z_1',PDBid='fragment')]
    out=tmp_path/'out';out.mkdir()
    provean.write_engine_outputs(cache,rows,out,{'17':['X','C2A']})
    text=(out/'provean_FRAGMENT_Z_1.out').read_text()
    assert 'C1A\t-1.25' in text and 'C99A' not in text


@pytest.mark.parametrize('file',['query.fasta','scores.tsv'])
def test_modified_provean_cache_rejected(tmp_path,file):
    cache=make_cache(tmp_path/'cache')
    (cache/file).write_text((cache/file).read_text()+'\n')
    with pytest.raises(ValueError):provean.validated_scores(cache)


def test_same_mutations_different_full_sequence_not_reused(tmp_path):
    cache=make_cache(tmp_path/'cache')
    sequence=tmp_path/'other.fasta';sequence.write_text('>other\nACDF\n')
    variants=tmp_path/'mut.tsv';variants.write_text('chain\tmutation\nX\tC2A\n')
    with pytest.raises(ValueError,match='new or empty'):
        provean.build_cache(sequence,'X',variants,cache)


def test_stop_marker_prevents_any_native_start(tmp_path):
    a=options(tmp_path);plan.build(a);(a.output/'STOP_AFTER_BATCH').touch()
    with patch.object(scheduler,'execute_task',side_effect=AssertionError('native start')), \
         patch.object(scheduler,'prepare_scan_cache',side_effect=AssertionError('PROVEAN start')):
        assert main(['scan','run','--resume','--output',str(a.output)])==0
    assert json.loads((a.output/'batch_status.json').read_text())['stage']=='paused'


def test_plan_only_does_not_compute(tmp_path):
    with patch.object(scheduler,'execute_task',side_effect=AssertionError('native start')), \
         patch.object(scheduler,'prepare_scan_cache',side_effect=AssertionError('PROVEAN start')):
        assert main(['scan','run','--pdb',str(TOOL_ROOT/'examples/1cse.pdb'),
                     '--partner1','E','--partner2','I','--sites','I:45',
                     '--output',str(tmp_path/'out'),'--plan-only'])==0


def test_output_lock_excludes_duplicate_scheduler(tmp_path):
    import fcntl
    out=tmp_path/'out';out.mkdir()
    with (out/'scan.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        assert main(['scan','run','--resume','--output',str(out)])==2
    assert not (out/'batch_status.json').exists()


def test_owned_cleanup_rejects_other_directory_or_symlink(tmp_path):
    root=tmp_path/'scratch';root.mkdir()
    outside=tmp_path/'outside';outside.mkdir()
    for p in (outside,root/'unsafe'):
        with pytest.raises(ValueError):scheduler.owned_workspace(p,root)
    link=root/'batch_0001_attempt_001';link.symlink_to(outside)
    with pytest.raises(ValueError):scheduler.owned_workspace(link,root)
    assert outside.is_dir()


def test_baseline_damaged_cache_rejected(tmp_path):
    a=options(tmp_path);m,_=plan.build(a)
    baseline=tmp_path/'baseline';(baseline/'wt_cache').mkdir(parents=True)
    (baseline/'wt_cache/file').write_text('changed')
    rows=scheduler.read_table(a.output/'plan/mutations.input')
    state=dict(stage='complete',wt_only=True,cap_mode='charged',protocol_sha256='same',
        pdb_sha256=plan.sha256(a.output/'plan/fragments.pdb'),
        partners=[rows[0]['Partner1'],rows[0]['Partner2']],wt_cache_sha256={'file':'wrong'})
    (baseline/'status.json').write_text(json.dumps(state))
    with pytest.raises(ValueError,match='checksum'):
        task.validate_baseline(baseline,a.output/'plan/mutations.input','charged','same')


def test_model_bytes_are_original():
    expected=json.loads((TOOL_ROOT/'docs/upstream_models.sha256.json').read_text())
    for name,checksum in expected.items():
        assert plan.sha256(RESOURCES/name)==checksum


def test_completed_orphan_is_recovered_without_native_launch(tmp_path):
    output=tmp_path/'scan'; output.mkdir()
    root=scheduler.workspace_root(output,tmp_path/'scratch')
    workspace=root/'batch_0001_attempt_001'; workspace.mkdir()
    destination=output/'batch'; destination.mkdir()
    saved={'stage':'complete','protocol_sha256':'protocol'}
    (workspace/'status.json').write_text(json.dumps(saved))
    (destination/'active_process.json').write_text(json.dumps({'workspace':str(workspace),'pid':None}))
    with patch.object(scheduler.subprocess,'Popen',side_effect=AssertionError('must adopt completed work')):
        actual,state=scheduler.execute_task(output,destination,'batch_0001',tmp_path/'input',root,
                                             {'protocol_sha256':'protocol'})
    assert actual==workspace and state==saved
    with pytest.raises(ValueError,match='protocol differs'):
        scheduler.execute_task(output,destination,'batch_0001',tmp_path/'input',root,
                               {'protocol_sha256':'different'})


def test_archive_before_cleanup_and_tampered_result_rejected(tmp_path):
    output=tmp_path/'scan'; output.mkdir()
    root=scheduler.workspace_root(output,tmp_path/'scratch')
    workspace=root/'batch_0001_attempt_001'; workspace.mkdir()
    destination=output/'batch'; destination.mkdir()
    row=dict(Result_Id='1',PDBfile='fragments.pdb',Partner1='A',Partner2='B',
             MutChain='B_1',Mutation_PDB='L45P',isPI='0',SunDDG='4.86',Confidence='High')
    scheduler.write_table(workspace/'results.tsv',list(row),[row])
    (workspace/'job').mkdir(); (workspace/'jobout').mkdir()
    (workspace/'job/job.input.cleaned.outdata').write_text('retained scientific features\n')
    (workspace/'status.json').write_text(json.dumps({'stage':'complete','protocol_sha256':'p','result_ids':['1']}))
    scheduler._archive_batch(workspace,destination,root,[row],'p')
    assert not workspace.exists()
    assert (destination/'scientific_intermediates.tar.gz').stat().st_size>0
    assert scheduler._validate_batch(destination,[row],'p')==[row]
    result=destination/'results.tsv'; result.write_text(result.read_text().replace('4.86','4.85'))
    with pytest.raises(ValueError,match='checksum'):
        scheduler._validate_batch(destination,[row],'p')


@pytest.mark.parametrize('command',['plan','aggregate'])
def test_plan_and_aggregate_cannot_change_active_scan(tmp_path,command):
    import fcntl
    args=options(tmp_path); plan.build(args)
    before=(args.output/'scan.json').read_bytes()
    with (args.output/'scan.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        assert main(['scan',command,'--resume','--output',str(args.output)] if command=='plan'
                    else ['scan',command,'--output',str(args.output)])==2
    assert (args.output/'scan.json').read_bytes()==before
