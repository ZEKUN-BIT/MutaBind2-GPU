#!/usr/bin/env python3
"""Compare CPU/GPU energies on identical archived WT and first mutant geometries."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import tarfile

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from mutabind2_tool.scan import task
from mutabind2_tool.scan.optimize import energy
from mutabind2_tool.core.runtime import load_config


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True,help='new lowercase ASCII directory')
    p.add_argument('--batch',type=int,default=1)
    p.add_argument('--rebuilt-pb',type=Path)
    p.add_argument('--rebuilt-vdw',type=Path)
    p.add_argument('--config',type=Path,help='tool configuration for the CPU and VDW executables')
    args=p.parse_args()
    config=load_config(args.config)
    if not config.get('charmm'):
        p.error('Configure the matching CPU CHARMM executable')
    source=args.scan.resolve(); output=energy._ascii_path(args.output)
    if output.exists(): raise ValueError('Use a new comparison output directory')
    output.mkdir(parents=True)
    baseline=source/'baseline'
    wt=next((baseline/'jobout/optimized_energy').glob('*/manifest.json'))
    cases=[('wt',json.loads(wt.read_text()),None)]
    archive=source/f'batches/batch_{args.batch:04d}/scientific_intermediates.tar.gz'
    with tarfile.open(archive) as tar:
        name=next(n for n in tar.getnames() if n.startswith('jobout/optimized_energy/')
                  and n.endswith('/manifest.json') and n.split('/')[-2]!=wt.parent.name)
        cases.append(('mutant',json.loads(tar.extractfile(name).read()),archive))
    results=[]
    for label,gpu,archive in cases:
        stem=gpu['context']['stem']; energy_label=gpu['context']['label']
        workspace=output/label; (workspace/'jobout').mkdir(parents=True)
        shutil.copytree(baseline/'sinputfiles',workspace/'sinputfiles')
        names=[f'jobout/{stem}_charmm.psf',f'jobout/{stem}_min.crd',f'jobout/energy_{energy_label}.inp']
        if archive:
            with tarfile.open(archive) as tar:
                for name in names: (workspace/name).write_bytes(tar.extractfile(name).read())
        else:
            for name in names: shutil.copyfile(baseline/name,workspace/name)
        for name in names:
            if energy._hash(workspace/name)!=gpu['context']['inputs'][name]:
                raise ValueError('Geometry/deck differs from saved GPU calculation')
        state={}; variants={'cpu':Path(config['charmm'])}
        if args.rebuilt_pb: variants['rebuilt_gpu']=args.rebuilt_pb.resolve()
        record={'case':label,'stem':stem,'reference_gpu':gpu['audits'],
                'reference_gpu_seconds':{k:v.get('wall_seconds') for k,v in gpu['processes'].items()},
                'reference_manifest':gpu}
        for backend,exe in variants.items():
            vdw=(args.rebuilt_vdw.resolve() if backend=='rebuilt_gpu' and args.rebuilt_vdw
                 else Path(config['vdw_executable']) if config.get('vdw_executable') else exe)
            calculated=task.run_pb(workspace,state,backend,energy.calculate_energy,
                workspace,stem,energy_label,exe,1,vdw)
            a,b=gpu['audits'],calculated['audits']
            delta={key: b['pb'][key]-a['pb'][key] for key in ('GAB','GA','GB','DELTAGELEC')}
            delta['VDWENER']=b['vdw']['VDWENER']-a['vdw']['VDWENER']
            delta['solve_energies']=[y-x for x,y in zip(a['pb']['energies'],b['pb']['energies'],strict=True)]
            passed=(a['pb']['grids']==b['pb']['grids'] and
                    max(abs(x) for x in [*delta['solve_energies'],delta['DELTAGELEC'],delta['VDWENER']])<=1e-4)
            record[backend]={'audits':b,'delta_vs_gpu':delta,'within_1e_4_kcal_mol':passed,
                             'timing':state['pb_runs'][backend],'manifest':calculated}
            if not passed: raise ValueError(f'{label}/{backend} differs by more than validation tolerance')
        results.append(record)
        (output/'comparison.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in r.items() if k in ('case','stem')}
                      for r in results],indent=2))


if __name__=='__main__': main()
