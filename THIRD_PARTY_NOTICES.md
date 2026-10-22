# Attribution and third-party components

## MutaBind2

The original MutaBind2 implementation, random-forest models and distributed resource files originate from [MinghuiLab/MutaBind2](https://github.com/minghuilab/MutaBind2). The original MIT copyright and permission notice is retained in `LICENSE`. The model files are unchanged; their SHA-256 checksums are recorded in `docs/upstream_models.sha256.json`. Existing notices and scientific references in force-field files are preserved.

This is an independently maintained GPU/workflow improvement by ZEKUN-BIT, not an official upstream release. No upstream endorsement is claimed.

## CHARMM integration

`native/pb/` and `native/reen/` contain CUDA implementations and integration patches; `native/vdw/` contains compact-allocation patches. Complete CHARMM Fortran source, object files, executables, and shared libraries are not distributed here. Users must obtain a compatible CHARMM source/build under its applicable license. Patch files do not grant rights to the underlying CHARMM code, which remains subject to its original terms. References to CHARMM routine names describe compatibility.

The repository MIT license covers original Python/workflow additions and newly authored CUDA code to the extent owned by the contributors; it does not relicense third-party content or underlying source referenced by integration patches.

## PROVEAN

`toolchain/vendor/provean-1.1.5/` is a separate GPL-licensed component, copyright 2012–2014 J. Craig Venter Institute. Its `COPYING`, `LICENSE`, source notices, and `MODIFICATIONS.md` are retained. The local compatibility changes affect sequence identifier handling and clustering integration; they do not create a new official PROVEAN version. PROVEAN is built separately and invoked as an external executable.

## External dependencies and data

VMD, CHARMM, FoldX, NAMD, BLAST+, CD-HIT, DSSP, CUDA, R packages and NR databases are external dependencies, each with its own terms. Their binaries, installers, runtime environments and databases are not bundled. Obtain and cite them from their respective distributors. The Python wheel packages the Python CLI and MutaBind2 resources only; native sources remain in the Git repository.

The 1CSE example structure originates from the Protein Data Bank and the upstream MutaBind2 example set. Public validation files are limited to this example and contain no private scan datasets.
