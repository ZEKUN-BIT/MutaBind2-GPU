"""Observe one production result and verify deployed backends and frozen input."""
from pathlib import Path
import json,re,shutil,sys,time
import continue_scan as c
R=Path(__file__).resolve().parent

def main():
    output=R/'validation';output.mkdir(exist_ok=True)
    deadline=time.monotonic()+1800
    observed_maps={}
    while time.monotonic()<deadline:
        transition=c.read(R/'transition.json')
        status=c.read(c.NEW/'pipeline_status.json')
        if status.get('stage')=='failed':raise RuntimeError(status)
        entry=status.get('energy') or {}
        if entry.get('workspace'):
            workspace=Path(entry['workspace'])
            # Record actual mapped native libraries in addition to loader identity.
            for proc in Path('/proc').iterdir():
                if not proc.name.isdecimal():continue
                try:
                    command=(proc/'cmdline').read_bytes()
                    if bytes(str(workspace/'jobout'),'utf8') not in command:continue
                    maps=(proc/'maps').read_text()
                    for mode,name in [('pb','charmm-pb-cuda-mstep'),('vdw','charmm-vdw-cuda-lists')]:
                        lines=[line for line in maps.splitlines() if name+'/lib/libchmm.so' in line]
                        if lines:observed_maps[mode]=dict(pid=int(proc.name),lines=lines)
                except (FileNotFoundError,PermissionError,ProcessLookupError):pass
            for path in sorted((workspace/'jobout/optimized_energy').glob('*/manifest.json')):
                manifest=c.read(path)
                if manifest.get('status')!='complete':continue
                context=manifest['context']
                deployed=c.read(R/'deployment.json')['runtimes']
                if context['runtime']!=deployed['pb'] or context['vdw_runtime']!=deployed['vdw']:
                    continue
                label=context['label']
                if label.endswith('_p2'):continue
                for relative,expected in context['inputs'].items():
                    assert c.digest(workspace/relative)==expected,(relative,'changed input')
                source=path.parent
                pblog=(source/'pb.log').read_text();vdwlog=(source/'vdw.log').read_text()
                assert pblog.count('PBEQ_MSTEP_CUDA_SUCCESS')==3
                assert 'MBVDW_DENSE_CUDA_SUCCESS' in vdwlog and 'MBVDW_SELECT_CUDA_SUCCESS' in vdwlog
                c.pipeline.energy._audit_log(source/'vdw.log','vdw',workspace,context['stem'])
                c.pipeline.energy._audit_log(source/'pb.log','pb',workspace,context['stem'],require_gpu=True)
                checkpoint=c.read(workspace/'gpu_migration.json')
                coordinates={p:h for p,h in checkpoint['scientific_files'].items() if p.endswith('_min.crd')}
                assert all(c.digest(workspace/p)==h for p,h in coordinates.items())
                protocol=c.read(c.NEW/'protocol.json')['sha256']
                assert protocol==transition['protocol_sha256']==status['protocol_sha256']
                case=output/'production_first';case.mkdir(exist_ok=True)
                for name in ('manifest.json','pb.log','vdw.log'):shutil.copy2(source/name,case/name)
                report=dict(observed_utc=c.core.now(),mutation=entry.get('mutation'),label=label,
                    energy_manifest=str(path),status='complete',processes=manifest['processes'],
                    protocol_sha256=protocol,loaded_libraries=observed_maps,
                    frozen_coordinate_files_verified=len(coordinates),
                    mstep_seconds=list(map(float,re.findall(r'PBEQ_STAGE MAYER_MSTEP seconds=\s*([0-9.]+)',pblog))),
                    production_combined_energy_seconds=sum(v['wall_seconds'] for v in manifest['processes'].values()),
                    preserved_completed_predictions=c.collect()[1]['completed'],
                    evidence_sha256={p.name:c.digest(p) for p in case.iterdir()})
                c.core.atomic_json(output/'production_activation.json',report)
                print(json.dumps(report,ensure_ascii=False,indent=2));return
        time.sleep(1)
    raise TimeoutError('Production activation not observed within 30 minutes')
if __name__=='__main__':main()
