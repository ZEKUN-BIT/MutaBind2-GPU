"""Wait for the unchanged SASA memory admission check before native execution."""
import json
import time
from pathlib import Path

from mutabind2_tool.scan import task
from mutabind2_tool.scan.optimize import surface

ADMISSION_ERROR = 'Surface stage requires 30 GiB MemAvailable at admission'
MINIMUM = 30 * surface.GIB


def install():
    original = surface.run_deck
    if getattr(original, '_waits_for_memory', False):
        return

    def admitted(executable, workspace, *args, **kwargs):
        workspace = Path(workspace)
        started = time.monotonic()
        record = workspace / 'surface_admission.json'
        waiting = False
        saved_stage = None

        def publish(stage, available):
            nonlocal saved_stage
            info = dict(stage=stage, updated_utc=task.now(),
                        minimum_bytes=MINIMUM, available_bytes=available,
                        wait_seconds=round(time.monotonic() - started, 3))
            task.atomic_json(record, info)
            path = workspace / 'status.json'
            if path.exists():
                state = json.loads(path.read_text())
                if saved_stage is None:
                    saved_stage = state.get('stage', 'surface')
                state['stage'] = stage if stage == 'waiting_surface_memory' else saved_stage
                state['surface_memory_admission'] = info
                task.atomic_json(path, state)

        while True:
            available = surface._available()
            if available < MINIMUM:
                waiting = True
                publish('waiting_surface_memory', available)
                time.sleep(1)
                continue
            if waiting:
                publish('admitted', available)
            try:
                return original(executable, workspace, *args, **kwargs)
            except RuntimeError as exc:
                # The native function checks again before mkdir/Popen. Retry
                # only this pre-launch rejection if memory changed meanwhile.
                if str(exc) != ADMISSION_ERROR:
                    raise
                waiting = True
                publish('waiting_surface_memory', surface._available())
                time.sleep(1)

    admitted._waits_for_memory = True
    surface.run_deck = admitted
