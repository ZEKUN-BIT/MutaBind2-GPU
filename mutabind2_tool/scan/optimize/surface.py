"""Exact CHARMM SASA for the target residue, using a conservative local context.

This adapter changes neither the shared templates nor CHARMM.  Its local sphere
set contains every possible occluder of every target atom.  Original radii,
coordinates, analytic SURFAC algorithm, atom order and DRASAP arithmetic remain.
Only mutant outputs are produced: the upstream WT surface output is unused.
"""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import resource
import shutil
import signal
import subprocess
import threading
import time

GENERATION = "charmm-analytic-local-sasa-v1"
STREAM_SHA256 = "24a112d68321f166f3ac40b29bb3c19900bf3c6f0af0afbe3128a0616c0e6468"
TEMPLATE_SHA256 = "c76b9c635d58d9611e0ba06948dd46769c2cc0da42c105f5b54bb504ccd5c189"
GIB = 1024 ** 3
TIMEOUT = 180
_RUNTIME_MODULE = None


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _ascii(path):
    value = os.path.abspath(os.fspath(path))
    if not re.fullmatch(r"[A-Za-z0-9_./-]+", value):
        raise ValueError(f"CHARMM requires an ASCII path without spaces: {value}")
    return Path(value)


def _name(value):
    if not re.fullmatch(r"[A-Za-z0-9_]+", str(value)):
        raise ValueError(f"Invalid CHARMM identifier: {value}")
    return str(value)


def _available():
    return int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                    if line.startswith('MemAvailable:'))) * 1024


def _limits():
    resource.setrlimit(resource.RLIMIT_AS, (20 * GIB, 20 * GIB))


def _runtime(executable):
    # Reuse the established read-only ELF/library fingerprint, without changing
    # optimized_energy.py or any existing energy manifest.
    global _RUNTIME_MODULE
    if _RUNTIME_MODULE is None:
        from . import energy
        _RUNTIME_MODULE = energy
    project_value = os.environ.get("MUTABIND2_TOOLCHAIN_ROOT")
    if not project_value:
        from mutabind2_tool.core.runtime import load_config
        config = load_config()
        charmm = config.get("charmm")
        project_value = config.get("toolchain_root")
    project = Path(project_value) if project_value else Path()
    native = Path(os.environ.get('MUTABIND2_CHARMM_NATIVE') or (project / 'tools/charmm'))
    capped = (Path(__file__).resolve().parents[1] / 'charmm_capped.py').resolve()
    text = Path(executable).read_text() if Path(executable).read_bytes()[:4] != b'\x7fELF' else ''
    passthrough = '#!/bin/sh\nexec ' + str(native) + ' "$@"\n'
    cap_wrapper = '#!/bin/sh\nexec ' + str(capped) + ' "$@"\n'
    extra = {}
    target = executable
    if text in {passthrough, cap_wrapper}:
        extra[str(executable)] = _hash(executable)
        target = native
        if text == cap_wrapper:
            extra[str(capped)] = _hash(capped)
    identity = dict(_RUNTIME_MODULE._runtime_identity(target))
    identity['surface_launch_wrappers'] = extra
    return identity


def make_deck(template, resources, local=True, diagnostics=False, original=False):
    resources = _ascii(resources)
    text = Path(template).read_text()
    marker = 'stream @pathpara/surface_residue.str'
    stream = (resources / 'surface_residue.str').read_text()
    canonical = text
    for name, placeholder in [('a', 'p1'), ('b', 'p2')]:
        canonical = re.sub(r'define ' + name + r' select segid ch\d+(?: -\n \.or\. segid ch\d+)* end',
                           f'define {name} select segid {placeholder} end', canonical)
    supported = (hashlib.sha256(stream.encode()).hexdigest() == STREAM_SHA256
                 and hashlib.sha256(canonical.encode()).hexdigest() == TEMPLATE_SHA256
                 and text.count(marker) == 1 and not text.split(marker)[-1].strip()
                 and 'read coor card name @path/@protname_charmm.crd' in text)
    if not supported:
        raise ValueError('Surface template changed; use the explicit original backend after reviewing the change')
    prefix = text.split(marker)[0]
    def dump(name):
        # A Z edit descriptor writes the exact IEEE-754 bits of each selected
        # WMAIN value.  Keep this instrumentation out of production decks.
        return (f'\nopen write card unit 90 name target_{name.lower()}.hex\n'
                'format (Z16.16)\nscalar wmain write 90 select r end\nformat\nclose unit 90\n')
    if original:
        body = stream
        if diagnostics:
            for name, command in [('MSR', 'set MSr ?area'), ('ASAC', 'set ASAC ?STOT'), ('ASAMP', 'set ASAMP ?STOT')]:
                body = body.replace(command, command + dump(name))
        return prefix + body + '\nstop\n', 'original'
    body = '''set rsol 1.4
coor surf rprobe @rsol select r end
set MSr ?area
'''
    if diagnostics:
        body += dump('MSR')
    if local:
        # SCALAR RADIUS reads the same VDWR(ITC(IAC)) as analytic SURFAC.
        # Any occluding sphere lies within 2*(max radius + probe) of a target
        # atom.  The extra Angstrom makes rounded scalar/selection thresholds
        # strictly conservative.  All spheres keep their original atom order.
        body += '''scalar radius statistics select c end
set srmax ?SMAX
set srmin ?SMIN
calc surcut = 2.0 * ( @srmax + @rsol ) + 1.0
if @srmin lt 0.0 stop
if @srmax gt 100.0 stop
define near select ( r .around. @surcut ) .or. r end
define lc select c .and. near end
define lp select mp .and. near end
set nlocal ?NSEL
'''
    for selection, name in [(('lc' if local else 'c'), 'ASAC'),
                            (('lp' if local else 'mp'), 'ASAMP')]:
        body += f'coor surf rprobe @rsol select {selection} end\n'
        body += f'scalar wmain statistics select r end\nset {name} ?STOT\n'
        if diagnostics:
            body += dump(name)
    # Keep the upstream separate divisions and subtraction (not (ASAMP-ASAC)/MSr).
    body += '''calc rASAc = @ASAC/@MSr
calc rASAmp = @ASAMP/@MSr
calc drASAp = @rASAmp - @rASAc
stop
'''
    return prefix + body, ('local' if local else 'full_context')


def audit(text, local=False):
    if ('ABNORMAL TERMINATION' in text or not re.search(r'^\s*NORMAL TERMINATION\b', text, re.M)
            or 'PBEQ>' in text or 'INTE EXTERN>' in text):
        raise ValueError('Surface process did not terminate normally or ran an unexpected energy stage')
    names = ['MSR', 'ASAC', 'ASAMP', 'RASAC', 'RASAMP', 'DRASAP']
    if local:
        names += ['SRMIN', 'SRMAX', 'SURCUT']
    values = {}
    for name in names:
        matches = re.findall(r'Parameter:\s+' + name + r'\s+<-\s+"([^"\n]+)"', text)
        if len(matches) != 1 or not math.isfinite(float(matches[0])):
            raise ValueError(f'Missing or nonfinite surface value: {name}')
        values[name] = matches[0]
    if float(values['MSR']) <= 0 or any(float(values[k]) < 0 for k in ['ASAC', 'ASAMP', 'DRASAP']):
        raise ValueError('Invalid surface area or interface difference')
    if local and not (0 <= float(values['SRMIN']) <= float(values['SRMAX']) <= 100):
        raise ValueError('Unsupported SASA radii')
    values['interface'] = 'yes' if float(values['DRASAP']) > 0 else 'no'
    return values


def run_deck(executable, workspace, stem, chain, residue, partner, case, deck):
    """Run only a surface deck; public for independent oracle comparisons."""
    executable, workspace, case = map(_ascii, (executable, workspace, case))
    stem, chain, residue, partner = map(_name, (stem, chain, residue, partner))
    if partner.lower() not in {'a', 'b'} or not residue.isdecimal():
        raise ValueError('Invalid partner or residue')
    if _available() < 30 * GIB:
        raise RuntimeError('Surface stage requires 30 GiB MemAvailable at admission')
    case.mkdir(parents=True, exist_ok=True)
    inp, log = case / 'surface.inp', case / 'surface.log'
    inp.write_text(deck)
    argv = ['/usr/bin/time', '-v', '-o', str(case / 'time.txt'), str(executable),
            f'pathinput={workspace / "sinputfiles"}', f'path={workspace / "jobout"}',
            f'protname={stem}', f'chainm={chain}', f'resnumber={residue}', f'partner={partner}']
    env = dict(os.environ, OMP_NUM_THREADS='1', MUTABIND2_CHARMM_THREADS='1')
    previous = None
    proc = None
    started = time.monotonic()
    def interrupted(signum, frame):
        raise InterruptedError(f'Surface interrupted by signal {signum}')
    try:
        if threading.current_thread() is threading.main_thread():
            previous = signal.signal(signal.SIGTERM, interrupted)
        with inp.open('rb') as stdin, log.open('wb') as stdout:
            proc = subprocess.Popen(argv, cwd=case, stdin=stdin, stdout=stdout, stderr=subprocess.STDOUT,
                                    env=env, preexec_fn=_limits, start_new_session=True)
            while proc.poll() is None:
                if time.monotonic() - started > TIMEOUT or _available() < 10 * GIB:
                    raise RuntimeError('Surface timeout or 10 GiB memory reserve reached')
                time.sleep(0.1)
            if proc.returncode:
                raise RuntimeError(f'Surface CHARMM exit code {proc.returncode}')
    except BaseException:
        if proc is not None:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait(timeout=5)
        raise
    finally:
        if previous is not None:
            signal.signal(signal.SIGTERM, previous)
    usage = (case / 'time.txt').read_text()
    return {'wall_seconds': time.monotonic() - started, 'returncode': proc.returncode,
            'max_rss_kib': int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', usage)[1]),
            'address_limit_bytes': 20 * GIB, 'reserve_bytes': 10 * GIB, 'timeout_seconds': TIMEOUT}


def calculate_surface(workspace, stem, chain, residue, partner, executable, output_name):
    workspace, executable = map(_ascii, (workspace, executable))
    stem, chain, residue, partner, output_name = map(_name, (stem, chain, residue, partner, output_name))
    resources, out = workspace / 'sinputfiles', workspace / 'jobout'
    template = out / f'surface_{stem}.inp'
    paths = [template, out / f'{stem}_charmm.psf', out / f'{stem}_charmm.crd',
             resources / 'surface_residue.str', resources / 'top_all36_prot.rtf', resources / 'par_all36_prot.prm']
    # Bound the physical coordinate domain used in the conservative-radius proof.
    coordinate_lines = paths[2].read_text().splitlines()
    expected_atoms = int(next(line.split()[0] for line in coordinate_lines if line.strip() and not line.startswith('*')))
    coordinates = [list(map(float, line.split()[4:7])) for line in coordinate_lines
                   if len(line.split()) >= 10 and line.split()[0].isdigit()]
    if (not coordinates or len(coordinates) != expected_atoms
            or any(not math.isfinite(v) or abs(v) > 10000 for xyz in coordinates for v in xyz)):
        raise ValueError('Invalid or unsupported surface coordinates')
    context = {'generation': GENERATION, 'implementation': _hash(__file__),
               'inputs': {str(path): _hash(path) for path in paths}, 'runtime': _runtime(executable),
               'selection': [stem, chain, residue, partner]}
    case = out / 'optimized_surface' / output_name
    deck, mode = make_deck(template, resources)
    (case / 'result.json').unlink(missing_ok=True)
    timing = run_deck(executable, workspace, stem, chain, residue, partner, case, deck)
    values = audit((case / 'surface.log').read_text(), local=mode == 'local')
    if context['inputs'] != {str(path): _hash(path) for path in paths}:
        raise ValueError('Surface inputs changed during calculation')
    if context['runtime'] != _runtime(executable):
        raise ValueError('Surface runtime changed during calculation')
    result = {'context': context, 'mode': mode, 'audit': values, 'process': timing,
              'deck_sha256': _hash(case / 'surface.inp'), 'log_sha256': _hash(case / 'surface.log')}
    (case / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    # Publish only after the process and all required values pass validation.
    destination = out / f'{output_name}.out'
    temporary = destination.with_name(destination.name + f'.tmp.{os.getpid()}')
    try:
        shutil.copyfile(case / 'surface.log', temporary)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return result


def install(engine, workspace):
    """Override only engine.surface; caller still runs engine.surfacefile first."""
    workspace = _ascii(workspace)
    def surface():
        results = []
        with Path(str(engine.in_file) + '.cleaned').open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                mutation = row['Mutation_cleaned']
                stem = row['PDBid'].lower() + '_' + mutation.lower()
                chain, residue = mutation[1], mutation[2:-1]
                partners = row['NewPartner1'] + row['NewPartner2']
                if partners.count(chain) != 1:
                    raise ValueError('Mutation chain is not uniquely mapped to a segment')
                segment = 'ch' + str(partners.index(chain) + 1)
                partner = 'a' if chain in row['NewPartner1'] else 'b'
                results.append(calculate_surface(workspace, stem, segment, residue, partner,
                               Path(engine.pathcharmm) / 'charmm', f'surface_{stem}_{chain.lower()}_{residue}'))
        return {'generation': GENERATION, 'mutants': results, 'wt_surface_skipped': True}
    engine.surface = surface
