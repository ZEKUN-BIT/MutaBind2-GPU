"""``mutabind2 scan`` command group."""

import argparse
import json
import fcntl
from pathlib import Path

from . import plan as plan_module


def add_plan_arguments(parser, output_required=True):
    parser.add_argument("--pdb", type=Path, help="input complex PDB")
    parser.add_argument("--partner1", help="comma-separated chain IDs of partner 1")
    parser.add_argument("--partner2", help="comma-separated chain IDs of partner 2")
    parser.add_argument("--mutations", type=Path,
                        help="TSV with columns chain<TAB>mutation[,id,note]")
    parser.add_argument("--sites", help="expand observed sites, e.g. A:10-20,A:30-35")
    parser.add_argument("--mutated-partner", choices=["1", "2"],
                        help="binding partner containing all mutations (inferred when omitted)")
    parser.add_argument("--fragment-mode", choices=["auto", "none", "explicit"], default="auto",
                        help="auto splits chains at physical peptide breaks (default)")
    parser.add_argument("--fragments", type=Path, help="explicit fragment JSON for --fragment-mode explicit")
    parser.add_argument("--gap-threshold", type=float, default=2.0,
                        help="C-N distance in Angstrom that defines a chain break (default 2.0)")
    parser.add_argument("--cap-mode", choices=["neutral", "charged"], default="charged",
                        help="original charged termini (default), or explicit neutral ACE/ACP + CT2 caps")
    parser.add_argument("--surface-backend", choices=["local", "original"], default="local",
                        help="optimized local CHARMM SASA (default) or original full surface")
    parser.add_argument("--provean-mode", choices=["fragment", "full"], default="fragment",
                        help="fragment PROVEAN (default) or a precomputed full-sequence cache")
    parser.add_argument("--provean-cache", type=Path, help="full-sequence PROVEAN cache directory")
    parser.add_argument("--max-per-batch", type=int, default=19,
                        help="maximum substitutions per batch (default 19)")
    parser.add_argument("--workers", type=int, default=2,
                        help="independent FoldX processes per batch (default 2)")
    parser.add_argument("--foldx-workers", type=int, choices=[1, 2], default=2)
    parser.add_argument("--pb-executable", type=Path,
                        help="optional validated split PB CHARMM executable")
    parser.add_argument("--backend", choices=["gpu", "cpu"], default="gpu",
                        help="GPU split PB (default) or seeded split CPU PB for validation")
    parser.add_argument("--vdw-executable", type=Path,
                        help="optional compact VDW CHARMM executable")
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--no-annotations", action="store_true",
                        help="skip partner/terminus distance annotations")
    if output_required:
        parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")


def _plan_summary(manifest):
    summary = plan_module.summary(manifest)
    return json.dumps(summary, ensure_ascii=False, indent=2)


def add_parser(subparsers):
    scan = subparsers.add_parser("scan", help="batch mutation scan over a complex structure")
    children = scan.add_subparsers(dest="scan_command", required=True)

    plan = children.add_parser("plan", help="validate input and write a resumable scan plan")
    add_plan_arguments(plan)

    run = children.add_parser("run", help="prepare baseline and run batches")
    add_plan_arguments(run)
    run.add_argument("--plan-only", action="store_true",
                     help="write the plan and print the batch schedule without calculations")
    run.add_argument("--baseline-only", action="store_true",
                     help="prepare and cache the WT baseline, then stop")
    run.add_argument("--only-batches", help="comma-separated batch indexes to run")
    run.add_argument("--max-batches", type=int, help="stop after this many new batches")
    run.add_argument("--stop-after-batch", action="store_true",
                     help="write STOP_AFTER_BATCH before starting; stop when it is seen")
    run.add_argument("--no-contributions", action="store_true")
    run.add_argument("--threads", type=int, default=8,
                     help="PROVEAN threads per sequence search (default 8)")
    run.add_argument("--workspace-root", type=Path,
                     help="lowercase ASCII scratch root; defaults to INSTALL_ROOT/workspaces")
    run.add_argument("--keep-workspace", action="store_true",
                     help="retain successful batch workspaces instead of archiving them")

    status = children.add_parser("status", help="show scan progress")
    status.add_argument("--output", type=Path, required=True)

    aggregate = children.add_parser("aggregate", help="rebuild predictions.tsv from checkpoints")
    aggregate.add_argument("--output", type=Path, required=True)

    provean = children.add_parser("provean",
                                  help="compute a full-sequence PROVEAN cache for --provean-mode full")
    provean.add_argument("--sequence", type=Path, required=True,
                         help="single-record FASTA of the complete mutated chain")
    provean.add_argument("--chain", required=True, help="original PDB chain for this sequence")
    provean.add_argument("--mutations", type=Path, required=True,
                         help="TSV with chain<TAB>mutation")
    provean.add_argument("--output", type=Path, required=True)
    provean.add_argument("--threads", type=int, default=8)
    provean.add_argument("--config", type=Path, default=None)
    return scan


def _fill_defaults(args):
    for name, value in (("fragment_mode", "auto"), ("gap_threshold", 2.0),
                        ("cap_mode", "charged"), ("surface_backend", "local"),
                        ("provean_mode", "fragment"), ("max_per_batch", 19),
                        ("workers", 2), ("foldx_workers", 2),
                        ("no_annotations", False), ("resume", False)):
        if not hasattr(args, name):
            setattr(args, name, value)


def dispatch(args):
    _fill_defaults(args)
    if args.scan_command == "plan":
        args.output = args.output.expanduser().resolve()
        args.output.mkdir(parents=True, exist_ok=True)
        with (args.output / "scan.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            manifest, _ = plan_module.build(args)
        print(_plan_summary(manifest))
        return 0
    if args.scan_command == "run":
        from .scheduler import run_scan
        return run_scan(args)
    if args.scan_command == "status":
        from .scheduler import print_status
        return print_status(args.output)
    if args.scan_command == "aggregate":
        from .scheduler import aggregate
        args.output = args.output.expanduser().resolve()
        manifest = plan_module.load(args.output)
        with (args.output / "scan.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            rows = aggregate(args.output, manifest)
        print(f"Aggregated {len(rows)} requests from checkpoints")
        return 0
    if args.scan_command == "provean":
        from .provean import build_cache
        result = build_cache(args.sequence, args.chain, args.mutations,
                             args.output, args.threads, args.config)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    raise SystemExit(f"Unknown scan command: {args.scan_command}")
