# Native acceleration sources

This directory is part of the unofficial MutaBind2-GPU improvement. It contains newly authored CUDA kernels, C ABI headers, source hash manifests and minimal integration patches. The complete CHARMM Fortran sources, build objects, executables and shared libraries must come from the user's authorized installation. Original CHARMM terms still apply to the underlying code referenced by the patches.

| Directory | Scope |
|---|---|
| `pb/` | PB SOR grid iteration and INTBP boundary-potential sums |
| `reen/` | REEN probe occlusion checks and negative-grid restoration |
| `vdw/` | Compact destination allocation for supported INTER3 partner selections |

PB defaults preserve mixed float/double expression behavior and disable CUDA fast math / fused multiply-add contraction. The supported boundary subset is zero salt, no membrane, no image atoms and a checked nonperiodic odd grid. The CPU retains the original interpolation, random-point generation, neighbor sorting and final full-grid cleanup. The REEN bridge uses the same generated random points on a fallback, preserving RNG progression.

Unsupported geometry/settings or exhausted REEN neighbor capacity use the checked CPU path. CUDA runtime errors abort and are not reported as successful fallback calculations. Native logs distinguish `PBEQ_CUDA_SUCCESS`, `PBEQ_BOUNDARY_CUDA_SUCCESS`, `PBEQ_REEN_CUDA_SUCCESS`, non-convergence and CPU fallback; inspect them when reporting GPU timings.

The VDW optimization changes allocation size for disjoint, exhaustive partner selections while retaining pair generation/filtering and energy kernels. Unsupported settings use the original allocation path. It is a CPU/memory optimization, separate from CUDA PB acceleration.

Use `scripts/build_native.py` from the repository root; see [installation instructions](../docs/INSTALL.md). It verifies original and patched file hashes and links against an existing CHARMM CMake build. All outputs go into a new directory. Source patching or compilation success alone does not establish numerical equivalence; run the comparison workflow with fixed geometries on your own installation.
