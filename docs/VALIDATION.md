# Validation scope

## Historical native measurements

The original local deployment was tested on Linux x86-64 with an RTX 5070 Ti. The public example is PDB 1CSE, partners E/I, mutation I:L45P. The source distribution preserves the original model bytes and calculation settings. It does not include the local licensed executables used for those measurements.

`validation/1cse_energy_comparison.json` is a compact extract of the 2026-10-21 native comparison, containing the measured energy audits, grids, residuals, iteration counts, timings and CPU/GPU differences. Host paths and dependency inventories have been omitted. The underlying comparison used identical PSF/coordinates, force-field inputs, seeds and grids for each geometry.

| Geometry | CPU VDW + PB | GPU VDW + PB | PB grid |
|---|---:|---:|---|
| WT | 50.188 s | 4.509 s | 211 × 225 × 243 |
| L45P | 46.577 s | 4.509 s | 211 × 225 × 243 |

The GPU time is the sum of the recorded VDW and PB subprocess times. Lock waits are excluded. These figures describe this calculation stage on this machine, not complete prediction throughput or a universal speedup. CPU/GPU energies agree at the logs' printed precision; no claim is made about equality of all floating-point bits. Rebuilt GPU-library comparisons are also included in the JSON.

The protocol retains 0.4 Å spacing, three component cold starts, a residual threshold of 2×10⁻⁶, an iteration limit of 2000, and explicit failure on non-convergence. WT iteration counts are `[891, 843, 884]`; mutant counts are `[891, 843, 883]`. The energy-comparison script requires matching grids and an energy difference no greater than 10⁻⁴ kcal/mol.

The complete charged-termini 1CSE deployment produced +4.83 kcal/mol and PROVEAN −2.191; the compact prediction table is included in `validation/1cse_predictions.tsv`. Those are historical outputs from the documented tool/database configuration. Independent FoldX/NAMD runs, different database releases and other toolchains can change the result.

## Public source tests

The Python test suite covers observed-site expansion, peptide continuity versus numbering gaps, partner-one mutations, residue/chain/result-ID mapping, WT identities, unchanged model hashes, cache corruption, mismatched full sequences, pause markers, scheduler exclusion, recovery and archive checks. It does not require licensed native tools or a GPU. CI runs the suite on Python 3.11, 3.12 and 3.13 and builds the Python wheel.

The source-release preparation also checks that the distributed native patches reconstruct the checksum-verified PB/VDW sources from the supported CHARMM baseline. Complete CHARMM source is supplied by the installing user. Native compilation and a numerical comparison are required when changing the native toolchain or GPU architecture.

The PROVEAN audit on 2026-10-22 reproduced −2.191 for the 63-aa 1CSE I-chain fragment (local mutation L38P), with 84 support sequences in 30 clusters. The unmodified PROVEAN binary gave the same score when using the same support set and cluster assignments, with full-ID formatting for modern blastdbcmd. This is an implementation/compatibility check, not biological accuracy validation.

## Limits of these results

Only a small public complex has been validated through the general native workflow. Multichain/break mappings have Python test coverage, not exhaustive native validation. The general GPU interface supports independent single substitutions; the combination-mutation GPU workflow has not been validated. Confidence labels come from the original model and are not experimental confidence intervals. See [known limitations](KNOWN_LIMITATIONS.md).
