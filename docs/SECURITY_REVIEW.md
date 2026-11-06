# Source publication review

Date: 2026-10-22. Reviewer: Codex primary task. Security-Sensitive: YES. Scope: source publication, portable native-build changes, configuration-based comparison paths, the separate PROVEAN build helper, package metadata and CI.

The CLI is for trusted local users. Configuration, native executables, source/build directories, model files and checkpoints are trusted installation inputs. This is not a service sandbox for untrusted users' programs. The review is limited to publication changes and does not claim a complete security audit of CHARMM, FoldX, NAMD, CUDA or all inherited upstream code.

| OWASP category | Publication review |
|---|---|
| A01 Access control | Build outputs use a new empty directory; original licensed source/library inputs remain unchanged. Existing scan ownership and locking tests pass. |
| A02 Cryptographic failures / exposure | No credentials, personal task datasets or host tool configuration are included. Git author metadata uses the GitHub no-reply address. |
| A03 Injection | Native build/patch/comparison operations use argument arrays. Build-helper paths are quoted. Architecture and native output-path restrictions are retained. |
| A04 Insecure design | Resource limits and trusted-local-input scope are retained. Unsupported GPU cases and known PROVEAN limitations are documented. |
| A05 Misconfiguration | Only a placeholder example configuration is tracked. CI has read-only repository permissions and pinned action revisions. |
| A06 Vulnerable components | The clean-install Python runtime dependency set was checked with pip-audit; no known advisories were reported for those resolved versions. Native dependencies require separate maintenance. |
| A07 Authentication failures | No login/session/authentication code is introduced. GitHub authentication is handled by the user's installed gh CLI. |
| A08 Integrity failures | Original model hashes are retained. Native patches reconstruct expected PB/VDW source hashes before compilation. GPL notices remain with PROVEAN. |
| A09 Logging | Build manifests/logs and scan checkpoints remain available locally; private execution logs are not part of the public repository. |
| A10 SSRF | No server-side request handling is introduced. External tool installation is explicit; predictions use the configured local database. |

The two unused upstream `pdb2crd` ELF resource files were removed from the public source copy. No complete CHARMM source, native ELF binaries, shared libraries, NAMD containers, Conda environments or NR databases are tracked. Newly authored CUDA sources and minimal integration patches are included with third-party notices.

Release checks: 29 Python tests passed in a new virtual environment; Python wheel built and checked for expected model/resource files; PROVEAN compiled from the distributed source; the compact VDW library compiled and passed loader/symbol checks; all three native source reconstruction hashes matched. The historical CPU/GPU numerical comparison is separately documented in `VALIDATION.md`; no new biological accuracy claim is made.

Dependency audit scope: NumPy 2.5.4, pandas 3.0.6, python-dateutil 3.0.0, six 1.17.0 and tzdata 2026.5, as resolved in the clean Python 3.12 installation. This describes the audited versions, not a promise that every version allowed by the package constraints is free of advisories.

No critical/high findings were identified in the publication changes. Status: PASS for this scope. Known PROVEAN quality/cache issues remain documented in `KNOWN_LIMITATIONS.md`.

## 2.1.0 synchronization review — 2026-11-06

Security-Sensitive: YES. Reviewed by the primary Codex task. Scope: native-build argument handling, incremental patch integrity, operational pipeline entry point, scratch-directory selection, result publication and local relocation.

| OWASP category | Assessment |
| --- | --- |
| A01 Access control | Local trusted-user CLI; existing output ownership, locks and checkpoint checks retained. Relocation preserved original bytes before cleanup. |
| A02 Cryptographic failures | SHA-256 checks cover patch inputs/outputs, result members and release assets. No credentials found in staged source/results. |
| A03 Injection | New commands use argument lists. CUDA architecture and generated executable/library paths are validated. No user input is evaluated as shell code. |
| A04 Insecure design | No increase in native concurrency or reduction in memory limits; numerical failures remain fatal. |
| A05 Misconfiguration | Machine configuration, environments, native binaries and databases are excluded from Git; wheels explicitly exclude legacy pdb2crd binaries. |
| A06 Vulnerable components | Audit of NumPy, pandas and their installed runtime dependencies found no known vulnerabilities. The broader old installation audit reported GHSA-xmj3-9gf7-h52w in setuptools 84.0.0, a Windows-only legacy wheel-installation path; this Linux tool does not use that path and the environment is not distributed. |
| A07 Authentication | No authentication changes or credential handling. |
| A08 Integrity | Exact CHARMM base/candidate checksums enforced. Historical protocol records are preserved; new scans have new fingerprints. Frozen WT/mutant comparisons matched all reported energy values. |
| A09 Logging | Resource retries preserve failed logs; release/migration inventories and result member checksums are retained. |
| A10 SSRF | No server-side network request feature added. |

Review status: PASS for the supported Linux deployment and published code/results. This review does not audit the full proprietary native toolchain. Validation is recorded in `validation/release_2.1.0_checks.json` and the release asset manifests.

Publication follow-up: GitHub identified Pillow advisories in two copied historical plotting-environment manifests. The installable lists now pin Pillow 13.0.0; both complete pinned lists pass pip-audit. Their original bytes remain in `historical_environment.freeze` and the immutable result archives. No prediction dependency, result or figure was changed. Security Review Status: ISSUES_FIXED for the newly published installable analysis manifests.
