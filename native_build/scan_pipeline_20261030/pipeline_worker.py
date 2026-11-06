"""Run the original batch engine, waiting for admission before its energy stage."""
import json
from pathlib import Path
import sys
import time

INSTALL = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(INSTALL))
from mutabind2_tool.scan import task


def wait_for_energy(workspace, state, operation):
    workspace = Path(workspace)
    gate = workspace / 'pipeline_energy_ready.json'
    started = time.monotonic()
    if not gate.exists():
        state.update(stage='waiting_pipeline_energy', pb_operation=operation,
                     pipeline_prepared_utc=task.now(), pipeline_energy_released=False)
        task.atomic_json(workspace / 'status.json', state)
        task.atomic_json(workspace / 'pipeline_prepared.json', {
            'prepared_utc': state['pipeline_prepared_utc'],
            'protocol_sha256': state['protocol_sha256']})
    while not gate.exists():
        time.sleep(0.5)
    admission = json.loads(gate.read_text())
    if admission['protocol_sha256'] != state['protocol_sha256']:
        raise ValueError('Pipeline energy admission has a different protocol')
    state.update(pipeline_energy_released=True,
                 pipeline_energy_wait_seconds=round(state.get('pipeline_energy_wait_seconds', 0)
                                                    + time.monotonic() - started, 3))


def main():
    if len(sys.argv) < 2 or sys.argv[1] != 'mutabind2_tool.scan.task':
        raise ValueError('The pipeline wrapper only supports the original batch task')
    original = task.run_pb
    def gated(workspace, state, operation, function, *parameters):
        wait_for_energy(workspace, state, operation)
        return original(workspace, state, operation, function, *parameters)
    task.run_pb = gated
    sys.argv = sys.argv[1:]
    task.main()


if __name__ == '__main__':
    main()
