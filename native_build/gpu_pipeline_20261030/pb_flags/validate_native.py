"""Compare frozen structures under the production energy lock and memory guards."""
from pathlib import Path
import fcntl
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
INSTALL = ROOT.parents[2]
sys.path.insert(0,str(INSTALL))
from mutabind2_tool.scan.optimize import energy


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    validation = ROOT/'validation'
    validation.mkdir(exist_ok=True)
    baseline = INSTALL/'projects/7mkj_128_255_20261030_optimized/baseline'
    wt = validation/'wt'
    if not wt.exists():
        (wt/'jobout').mkdir(parents=True)
        (wt/'sinputfiles').symlink_to(baseline/'sinputfiles',target_is_directory=True)
        for name in ('fragments_charmm.psf','fragments_min.crd','energy_fragments_p2.inp'):
            shutil.copy2(baseline/'wt_cache'/name,wt/'jobout'/name)
    mutant = INSTALL/'native_build/pb_reen_warp_20261030/validation/scan_mutant'
    reports = {}
    for name,workspace,stem,label in [('wt',wt,'fragments','fragments_p2'),
            ('mutant',mutant,'fragments_ve41a','fragments_ve41a')]:
        paths = [workspace/'jobout'/f'{stem}_{suffix}' for suffix in ('charmm.psf','min.crd')]
        paths += [workspace/'jobout'/f'energy_{label}.inp']
        hashes = {str(p):digest(p) for p in paths}
        deck = energy._make_decks(workspace,label)['pb']
        reference = ROOT.parent/'pb/validation'/f'{name}_cells/result.json'
        reference_result=json.loads(reference.read_text())
        assert sorted(reference_result['source_sha256'].values())==sorted(hashes.values())
        case_reports = {'cells':reference_result}
        # Reverse the order on the second structure to reduce ordering bias.
        backends = [('mstep',ROOT/'bin/charmm')]
        if name=='mutant':
            backends.reverse()
        for backend,executable in backends:
            case = validation/f'{name}_{backend}'
            case.mkdir(exist_ok=False)
            def save_phi(match):
                component=match[1].lower()
                return (match[0]+f'\nopen write file unit 40 name {case/(component+".phi")}\n'
                        'write phi unit 40\nclose unit 40')
            saved,count = re.subn(r'^SET (Gab|Ga|Gb) = \?enpb[^\n]*$',save_phi,deck,flags=re.M)
            assert count==3
            (case/'pb.inp').write_text(saved)
            print(json.dumps(dict(waiting_for_energy_lock=str(case))),flush=True)
            with (INSTALL/'runtime/energy.lock').open('a') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX)
                gpu = subprocess.check_output(['nvidia-smi','--query-gpu=memory.free',
                    '--format=csv,noheader,nounits'],text=True)
                if min(map(int,gpu.split()))<2048:
                    raise RuntimeError('Need 2 GiB free GPU memory')
                print(json.dumps(dict(starting=str(case))),flush=True)
                result = energy._run_process(executable,workspace,stem,case,'pb',1)
            result['audit'] = energy._audit_log(case/'pb.log','pb',workspace,stem,require_gpu=True)
            text = (case/'pb.log').read_text()
            result['stage_seconds'] = {}
            for stage,seconds in re.findall(r'PBEQ_STAGE\s+(\w+)\s+seconds=\s*([\d.Ee+-]+)',text):
                result['stage_seconds'].setdefault(stage,[]).append(float(seconds))
            result['reen_gpu_solves'] = len(re.findall(r'^\s*PBEQ_REEN_CUDA_SUCCESS\b',text,re.M))
            result['boundary_gpu_solves'] = len(re.findall(r'^\s*PBEQ_BOUNDARY_CUDA_SUCCESS\b',text,re.M))
            result['cell_entries'] = list(map(int,re.findall(r'PBEQ_REEN_CELL_CANDIDATES entries=(\d+)',text)))
            assert result['reen_gpu_solves']==result['boundary_gpu_solves']==3
            if backend=='mstep':
                assert text.count('PBEQ_MSTEP_CUDA_SUCCESS')==3,'MSTEP GPU fallback'
                assert len(result['cell_entries'])==3,'Spatial index unexpectedly fell back'
            result['phi_sha256'] = {p.name:digest(p) for p in sorted(case.glob('*.phi'))}
            assert len(result['phi_sha256'])==3
            assert all(digest(Path(p))==v for p,v in hashes.items())
            result['source_sha256']=hashes
            result['reference_run']=str(reference)
            result['runtime']=energy._runtime_identity(executable)
            (case/'result.json').write_text(json.dumps(result,indent=2)+'\n')
            case_reports[backend]=result
            print(json.dumps(dict(completed=str(case),wall_seconds=result['wall_seconds'],
                  neighbor_seconds=sum(result['stage_seconds']['REEN_CUDA_NEIGHBOR_PREP']))),flush=True)
        assert case_reports['mstep']['audit']==case_reports['cells']['audit'],name
        assert case_reports['mstep']['phi_sha256']==case_reports['cells']['phi_sha256'],name
        reports[name]=dict(bitwise_phi_equal=True,logged_results_equal=True,backends=case_reports)
        (validation/'comparison.json').write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps(dict(passed=True,cases=list(reports))),flush=True)


if __name__=='__main__':
    main()
