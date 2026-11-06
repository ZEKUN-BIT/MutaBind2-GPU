from pathlib import Path
import argparse,fcntl,hashlib,json,re,subprocess,sys,time
R=Path(__file__).resolve().parent;P=R.parents[2];sys.path.insert(0,str(P))
from mutabind2_tool.scan.optimize import energy
GIB=1024**3

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def run(name,workspace,stem,deck,backend,variant):
    case=R/'validation'/name;case.mkdir(parents=True,exist_ok=False);(case/'vdw.inp').write_text(deck)
    executable=R/variant/'bin'/('charmm_'+backend) if backend in ('gpu','cpu') else P/'bin/charmm-vdw-compact'
    paths=[workspace/'jobout'/f'{stem}_{s}' for s in ('charmm.psf','min.crd')];hashes={str(p):sha(p) for p in paths}
    print(json.dumps(dict(starting=name)),flush=True)
    result=energy._run_process(executable,workspace,stem,case,'vdw',1)
    result['audit']=energy._audit_log(case/'vdw.log','vdw',workspace,stem)
    text=(case/'vdw.log').read_text()
    result['hashes']={key:dict(count=int(n),sha256=h) for key,n,h in re.findall(r'MBVDW_SHA256\s+(\w+)\s+(\d+)\s+([0-9A-F]{64})',text)}
    result['energy_bits']={i:bits for i,bits in re.findall(r'MBVDW_ENERGY_BITS\s+(\d+)\s+([0-9A-F]{16})',text)}
    result['printed']={key:re.findall(r'Parameter:\s+'+key+r'\s+<-\s+"([^"]+)"',text)[-1] for key in ('VDWENER','ELECENER')}
    result['dense_gpu']=list(map(int,re.findall(r'MBVDW_DENSE_CUDA_SUCCESS pairs=(\d+)',text)))
    result['selected_gpu']=list(map(int,re.findall(r'MBVDW_SELECT_CUDA_SUCCESS pairs=(\d+)',text)))
    result['negative_pairs']=list(map(int,re.findall(r'MBVDW_NEGATIVE_PAIRS\s+(\d+)',text)))
    result['runtime']=energy._runtime_identity(executable);result['source_sha256']=hashes
    assert all(sha(p)==h for p,h in hashes.items())
    if variant=='audit' and backend!='installed':assert len(result['hashes'])>=6 and len(result['energy_bits'])==192
    (case/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(completed=name,seconds=result['wall_seconds'],dense=result['dense_gpu'],selected=result['selected_gpu'])),flush=True)
    return result

def pair(label,w,stem,deck,variant,expect_dense=True,expect_select=True):
    old=run(label+'_cpu',w,stem,deck,'cpu' if variant=='audit' else 'installed',variant)
    new=run(label+'_gpu',w,stem,deck,'gpu',variant)
    assert old['printed']==new['printed'] and old['audit']==new['audit'],label
    if variant=='audit':
        assert old['hashes']==new['hashes'],(label,'list/force hashes')
        assert old['energy_bits']==new['energy_bits'],(label,'energy bits')
    assert bool(new['dense_gpu'])==expect_dense,(label,'dense')
    assert bool(new['selected_gpu'])==expect_select,(label,'selected')
    return dict(label=label,passed=True,bitwise_equal=variant=='audit',backends=dict(cpu=old,gpu=new))
def main():
    parser=argparse.ArgumentParser();parser.add_argument('phase',choices=['small','large_audit','large_performance']);args=parser.parse_args()
    small=Path('/home/zek/mutabind2_runs/7mkj_rpod/optimization_adapter_small')
    deck=(small/'jobout/optimized_energy/mutabinds_p2/vdw.inp').read_text()
    reports=[]
    if args.phase=='small':
        cases=[('wt','mutabinds',deck,True,True),('mutant','mutabinds_lb38p',deck,True,True),
            ('cutoff','mutabinds',deck.replace('update cutnb 999 switch vswitch eps @EpsG','update cutnb 14 ctonnb 10 ctofnb 12 switch vswitch eps @EpsG'),False,True),
            ('covalent','mutabinds',deck.replace('define a select segid ch2 end\ndefine b select segid ch1 end','define a select bynu 1:1000 end\ndefine b select .not. a end'),True,True),
            ('overlap','mutabinds',deck.replace('define b select segid ch1 end','define b select all end'),True,False),
            ('partial','mutabinds',deck.replace('define b select segid ch1 end','define b select segid ch1 .and. resid 1:50 end'),True,False),
            ('keep','mutabinds',deck.replace('interaction sele A end sele B end','interaction keep sele A end sele B end'),True,False),
            ('shift','mutabinds',deck.replace('update cutnb 999 switch vswitch','update cutnb 999 shift vshift'),True,False),
            ('bycubes','mutabinds',deck.replace('update cutnb 999 switch vswitch','update cutnb 999 switch vswitch bycu'),False,False),
            ('fixed','mutabinds',deck.replace('update cutnb 999','cons fix sele bynu 1:100 end\nupdate cutnb 999'),False,True)]
        for name,stem,inp,dense,selected in cases:reports.append(pair('small_'+name,small,stem,inp,'audit',dense,selected))
        assert reports[3]['backends']['gpu']['negative_pairs'][0]>0
    else:
        variant='audit' if args.phase=='large_audit' else 'production'
        wt=P/'native_build/gpu_pipeline_20261030/pb/validation/wt'
        mutant=P/'native_build/pb_reen_warp_20261030/validation/scan_mutant'
        for name,w,stem,label in [('wt',wt,'fragments','fragments_p2'),('mutant',mutant,'fragments_ve41a','fragments_ve41a')]:
            inp=energy._make_decks(w,label)['vdw'];reports.append(pair(args.phase+'_'+name,w,stem,inp,variant))
    report=dict(passed=True,phase=args.phase,cases=reports)
    (R/'validation'/(args.phase+'_report.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(passed=True,phase=args.phase)),flush=True)
if __name__=='__main__':
    with (P/'runtime/energy.lock').open('a') as lock:
        print('Waiting for production energy lock',flush=True);fcntl.flock(lock,fcntl.LOCK_EX)
        gpu=subprocess.check_output(['nvidia-smi','--query-gpu=memory.free','--format=csv,noheader,nounits'],text=True)
        assert min(map(int,gpu.split()))>=2048
        main()
