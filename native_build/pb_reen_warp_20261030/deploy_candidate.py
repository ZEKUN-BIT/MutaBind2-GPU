"""Install a separate immutable backend after numerical and timing validation."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
INSTALL = Path('/home/zek/mutabind2-tool')
sys.path.insert(0, str(INSTALL))
from mutabind2_tool.scan.optimize.energy import _runtime_identity


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


for name in ['cuda', 'errors', 'bridge', 'no_device']:
    assert json.loads((HERE / f'validation/{name}.json').read_text())['all_passed']
report = json.loads((HERE / 'validation/scan_mutant_comparison.json').read_text())
assert report['three_phi_files_bitwise_equal'] and report['logged_results_identical']
assert report['backends']['candidate']['wall_seconds'] < .95 * report['backends']['baseline']['wall_seconds']
build = json.loads((HERE / 'build_manifest.json').read_text())
assert digest(HERE / 'lib/libchmm.so') == build['library_sha256']
assert digest(HERE / 'pb_reen_cuda.cu') == build['source_sha256']
old = INSTALL / 'toolchain/external/charmm-pb-cuda-reen/lib/libchmm.so'
assert digest(old) == '2d3ba37325a46e472caf1624bded5c100348edb16517154941217123f836c16e'
target = INSTALL / 'toolchain/external/charmm-pb-cuda-reen-warp'
target.mkdir(exist_ok=False)
(target / 'lib').mkdir()
(target / 'source').mkdir()
shutil.copy2(HERE / 'lib/libchmm.so', target / 'lib/libchmm.so')
for name in ['pb_reen_cuda.cu', 'pb_reen_cuda.h', 'source_original.cu', 'build_candidate.py']:
    shutil.copy2(HERE / name, target / 'source' / name)
launcher = INSTALL / 'bin/charmm-pb-gpu-warp'
assert not launcher.exists()
launcher.write_text((HERE / 'bin/charmm').read_text().replace(str(HERE / 'lib'), str(target / 'lib')))
launcher.chmod(0o755)
runtime = _runtime_identity(launcher)
assert runtime['pb_backend'] == 'cuda'
assert runtime['files'][str(target / 'lib/libchmm.so')] == build['library_sha256']
assert digest(old) == '2d3ba37325a46e472caf1624bded5c100348edb16517154941217123f836c16e'
record = dict(deployed_utc=datetime.now(timezone.utc).isoformat(), runtime=runtime,
              source_build=str(HERE), validation=str(HERE / 'validation'),
              old_backend_unchanged=True)
(target / 'DEPLOYMENT.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'launcher': str(launcher), 'library_sha256': build['library_sha256']}, indent=2))
