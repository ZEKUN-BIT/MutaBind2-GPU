# PB REEN spatial candidates, 2026-10-30

The new Fortran bridge builds ordered candidate atom lists for occupied spatial
cells before the existing REEN neighbor calculation. It keeps the original
distance and SEGPAR operations, sorting order, random stream, CUDA kernels,
physical parameters, precision, grids and solver tolerance.

Cell width is four times the largest expanded radius. Spheres whose cells are
not adjacent are separated by at least that width, whereas intersecting spheres
can be at most twice the largest expanded radius apart. This leaves a factor of
two margin on the validated coordinate/radius domain. Candidate lists contain
the original selected indices in increasing order. Small, sparse, unsupported
or oversized cases retain the existing scan. Allocation failure also falls back.

## Measured result

Both comparisons used the same frozen PSF/CRD and one CPU thread. They held
`runtime/energy.lock`, so benchmark energy calculations did not overlap production
energy calculations. The existing 30 GiB admission, 20 GiB address limit and
10 GiB memory reserve remained in force.

| Frozen structure | Previous neighbor preparation | New neighbor preparation | Previous complete PB check | New complete PB check |
| --- | ---: | ---: | ---: | ---: |
| Current scan WT | 6.369 s | 1.652 s | 57.020 s | 52.020 s |
| Archived V130A | 6.319 s | 1.638 s | 57.026 s | 52.019 s |

The complete checks include writing three PHI files. Their wall times are about
8.8% shorter. The isolated neighbor preparation is about 74% shorter. Production
previously spent about 66.7 seconds per mutant on VDW plus PB; saving about
4.7–5 seconds predicts roughly 7% less time for those stages together. This is a
two-structure estimate, not a measurement across the full mutation table.

The continuation activated at 13:17:26 Asia/Shanghai. Its first production energy
calculation (S133A) completed with 15.511 s VDW plus 48.023 s PB, or 63.534 s
combined; neighbor preparation totaled 1.643 s. The prior run's typical combined
time was about 66.7 s on other mutants. That historical comparison is not a
controlled comparison of identical inputs, and a new batch average is still
pending. The active native process's loaded library was checked directly.

All six full PHI files have identical SHA-256 hashes between implementations.
Logged energies, grids, residuals and iteration counts agree exactly. Current WT
coordinates and its previous cached PB audit also match this validation.

Evidence: `validation/comparison.json`, `native_validation.log`,
`cell_test_results.json`, `transition_test_results.json`, `transition_tests.log`.
Production evidence is in `validation/production_activation.json` and
`validation/production_first/`.
There are 46 compiled Fortran index cases and five continuation tests. The latter
include a real isolated continuation plan, all 80 frozen WT cache files, 19
prepared mutant geometries and rejection of an altered coordinate checkpoint.

## Deployment and continuation

The previous installed libraries and scientific Python package are untouched.
The new installed backend is `bin/charmm-pb-gpu-cells`, backed by
`toolchain/external/charmm-pb-cuda-cells/lib/libchmm.so`.

`continue_scan.py run` waits for the preceding scan's current batch to finish,
then creates `projects/7mkj_128_255_20261030_cells` with a new scientific protocol
fingerprint. It keeps the exact frozen WT cache, records its original protocol
and numerical equivalence evidence, and migrates the waiting next batch's frozen
geometries. It performs no new WT or donated mutant minimization. All previously
completed predictions keep their original scan/protocol provenance in
`projects/7mkj_128_255_20261030_combined`.

The existing batch pipeline is reused with its recorded source hashes. Energy
concurrency remains one and CPU worker settings are unchanged. Original scans
retain stop markers so resuming them cannot duplicate the continuation.

From `/home/zek/mutabind2-tool`, inspect combined progress with:

```bash
.env/bin/python native_build/pb_cells_20261030/continue_scan.py status
```

Resume this continuation controller after an interruption with:

```bash
nohup .env/bin/python -u native_build/pb_cells_20261030/continue_scan.py run >> native_build/pb_cells_20261030/continuation.log 2>&1 &
```

An already running controller holds `coordinator.lock` and rejects a duplicate.
Runtime libraries, immutable plans and recorded continuation scripts must retain
their recorded hashes. This command does not reuse energy results from an
interrupted native calculation; completed batch archives remain reusable.

During migration the old scan's stop marker belongs to the transition. To pause
the new scan cleanly after its current batch once it is active:

```bash
touch projects/7mkj_128_255_20261030_cells/STOP_AFTER_BATCH
```

Remove that new scan's stop marker before explicitly resuming a requested pause.
Keep the stop markers on the two older scans.
