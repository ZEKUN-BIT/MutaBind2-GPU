"""Read production logs and save a performance audit without running calculations."""

from collections import defaultdict
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "projects/7mkj_256_374_20261102"
HERE = Path(__file__).resolve().parent


def read_json(path):
    return json.loads(path.read_text())


def main():
    state = read_json(OUTPUT / "batch_status.json")
    completed = state["batches_complete"]
    rows = []
    for index in range(max(1, completed - 9), completed + 1):
        path = OUTPUT / f"batches/batch_{index:04d}/status.json"
        status = read_json(path)
        timings = status["timings_seconds"]
        groups = {}
        for prefix in ("topology_", "minimize_", "geometry_", "mutant_pb_energy_"):
            values = [v for k, v in timings.items() if k.startswith(prefix)]
            groups[prefix] = {
                "count": len(values),
                "sum_seconds": sum(values),
                "mean_seconds": sum(values) / len(values) if values else None,
            }
        rows.append({
            "batch": index,
            "source": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "timings_seconds": timings,
            "per_mutant_groups": groups,
        })

    sample = None
    workspace = Path(state["workspace"])
    for path in sorted((workspace / "jobout/optimized_energy").glob("*/pb.log")):
        try:
            raw = path.read_bytes()
        except FileNotFoundError:
            continue
        text = raw.decode(errors="replace")
        if "NORMAL TERMINATION BY NORMAL STOP" not in text:
            continue
        grids = re.findall(r"grid (\d+)x(\d+)x(\d+)", text)
        iterations = [int(v) for v in re.findall(
            r"PBEQ_CUDA_SUCCESS iterations=\s*(\d+)", text)]
        elapsed = re.findall(r"ELAPSED TIME:\s*([\d.]+)\s+SECONDS", text)
        if len(grids) != 3 or len(iterations) != 3 or not elapsed:
            continue
        stages = defaultdict(float)
        for name, seconds in re.findall(r"PBEQ_STAGE (\w+) seconds=\s*([\d.]+)", text):
            stages[name] += float(seconds)
        sample = {
            "source": str(path),
            "mutation_directory": path.parent.name,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "grids": [[int(v) for v in grid] for grid in grids],
            "iterations": iterations,
            "elapsed_seconds": float(elapsed[-1]),
            "stage_seconds": dict(stages),
            "stage_accounting": "Stage timers are nested; do not sum all stage values.",
        }
        (HERE / "sample_pb_profile.log").write_text("\n".join(
            line for line in text.splitlines() if any(key in line for key in (
                "PBEQ_STAGE", "PBEQ_CUDA_SUCCESS", "PBEQ1 CUDA:", "ELAPSED TIME:"
            ))
        ) + "\n")
        break
    if sample is None:
        raise RuntimeError("No completed native PB log in the active workspace; no benchmark was launched.")

    nx, ny, nz = sample["grids"][0]
    sample["first_grid_points"] = nx * ny * nz
    sample["speedup_required_for_100ms"] = sample["elapsed_seconds"] / 0.1
    sample["speedup_required_for_1ms"] = sample["elapsed_seconds"] / 0.001
    source = ROOT / "toolchain/external/charmm-pb-cuda/source/pbeq_cuda.cu"
    result = {
        "captured_beijing": datetime.now(timezone(timedelta(hours=8))).isoformat(),
        "constraints": {
            "retain_physical_methods": True,
            "retain_precision": True,
            "surrogate_model_authorized": False,
        },
        "method": "Read existing archived timings and one completed PB log; no new physical calculation.",
        "timing_scope": {
            "foldx_vmd_charmm_surface_features_predict_contributions": "One batch of 19 mutations.",
            "per_mutant_groups": "Measured individual tasks; counts, totals and means are separate.",
            "mutant_pb_energy_": "Combined VDW and PB process time, despite the historical field name.",
            "preparation": "Normally overlaps the preceding batch's energy work; do not add to its wall time.",
        },
        "scan_snapshot": state,
        "recent_completed_batches": rows,
        "sample_pb": sample,
        "solver_source": {
            "path": str(source),
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "observed": [
                "All six large grids remain on device throughout each solve.",
                "Six grids are allocated and uploaded at each solve entry.",
                "Residual reduction and synchronous scalar download occur every iteration.",
                "Potential grid download and host commit occur after convergence.",
                "Explicit float/double rounded operations preserve expression evaluation.",
            ],
        },
        "limitations": [
            "Observed times do not establish a hardware lower bound or an optimized speedup.",
            "Native stage timers are nested and cannot all be added together.",
            "Batch throughput and latency of a new single mutation are different metrics.",
            "Cached retrieval latency is not the latency of a new physical calculation.",
        ],
    }
    (HERE / "profile.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({
        "profile": str(HERE / "profile.json"),
        "completed_batches_sampled": len(rows),
        "sample_pb_seconds": sample["elapsed_seconds"],
        "grid_points": sample["first_grid_points"],
        "iterations": iterations,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
