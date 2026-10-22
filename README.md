# MutaBind2-GPU · GPU 加速改进版

**An unofficial GPU-accelerated improvement of [MutaBind2](https://github.com/minghuilab/MutaBind2), maintained by [ZEKUN-BIT](https://github.com/ZEKUN-BIT).**

**本项目是基于原始 MutaBind2 的非官方 GPU 加速改进版，保留原模型与上游署名；不代表原作者发布的新版本。** [中文说明](README.zh-CN.md)

MutaBind2-GPU estimates mutation-induced changes in protein–protein binding affinity from a complex structure. It retains the original random-forest models and adds CUDA acceleration to selected CHARMM Poisson–Boltzmann stages, resource controls, resumable mutation scans, and cached calculations. It is an engineering extension, rather than a newly trained prediction model.

## Improvements

- CUDA kernels for supported PB SOR grid iterations, boundary-potential sums, and REEN probe evaluation/grid restoration; the Fortran bridge checks supported conditions and reports CPU fallback.
- Compact VDW destination allocation, local analytical SASA, cached wild-type calculations, and at most two independent FoldX processes.
- Python 3 CLI for individual substitutions and saturation scans, with original chain/residue mapping, explicit chain-break handling, checkpoints, and scientific intermediate archives.
- PROVEAN 1.1.5 compatibility changes for modern BLAST sequence identifiers, plus per-sequence score caching and an optional complete-sequence workflow.

**The general GPU scan workflow has been validated for independent single substitutions.** The legacy combination-mutation interface remains available; combination-mutation GPU support has not been validated. PROVEAN and the final random forest run on the CPU. Native tools, licensed binaries, environments, and NR databases must be installed separately.

## Measured validation

Linux x86-64, RTX 5070 Ti, CHARMM 51, fixed 1CSE WT/L45P geometries:

| Case | CPU VDW + PB | GPU VDW + PB | CPU/GPU comparison |
|---|---:|---:|---|
| WT | 50.19 s | 4.51 s | Energies agree at printed precision |
| I:L45P | 46.58 s | 4.51 s | Energies agree at printed precision |

These are timings of the **VDW + PB stage**, excluding lock waits, and are not an end-to-end prediction speedup. Grid spacing remains 0.4 Å; each component starts a cold solve, with a residual threshold of 2×10⁻⁶ and explicit convergence checks. The comparison uses identical coordinates, force-field inputs, seeds, and grids. See [validation details](docs/VALIDATION.md) and [machine-readable energy evidence](validation/1cse_energy_comparison.json). This engineering comparison is not an independent biological accuracy benchmark.

## Install

Use Linux, Python 3.11+, an NVIDIA GPU/CUDA toolchain, and separately installed native dependencies. The GPU bridge currently targets a hash-checked CHARMM 51 source baseline; obtain CHARMM from its authorized distributor.

```bash
git clone https://github.com/ZEKUN-BIT/MutaBind2-GPU.git mutabind2-gpu
cd mutabind2-gpu
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[test]'
cp config/tools.example.json config/tools.json
```

Edit `config/tools.json` with your own executable and NR database paths. Build the native GPU libraries using [the installation guide](docs/INSTALL.md), then run:

```bash
./mutabind2 doctor --scope gpu
python -m pytest -q
```

The Python tests run without proprietary tools or a GPU. A complete prediction requires VMD, CHARMM, NAMD, DSSP, FoldX, PROVEAN, BLAST+, CD-HIT, an NR database, and R with `randomForest` and `forestFloor`. Installing the Python package alone does not install these tools.

## Predict or scan

The included public example is PDB **1CSE**, partners E and I, substitution I:L45P:

```bash
./mutabind2 scan run \
  --pdb examples/1cse.pdb --partner1 E --partner2 I \
  --mutations examples/mutations.tsv \
  --output projects/1cse_scan --plan-only

./mutabind2 scan run --resume --output projects/1cse_scan
```

A mutation table is tab separated:

```text
chain	mutation	id
I	L45P	1
```

Each row is an independent single substitution. For a saturation scan, replace `--mutations` with, for example, `--sites 'I:45-48'`. Each observed site expands to 19 non-WT substitutions. Missing coordinates are not rebuilt. Inputs may cover multiple chains within the same binding partner; split requests that mutate opposite partners.

`completed_predictions.tsv` contains completed predictions. Positive ΔΔG means weaker binding; negative ΔΔG means stronger binding, in kcal/mol. `PROVEAN` is a sequence-variation feature, not a binding ΔΔG. Run status is available through:

```bash
./mutabind2 scan status --output projects/1cse_scan
touch projects/1cse_scan/STOP_AFTER_BATCH
```

Remove the stop marker before resuming. Use a new output/cache directory after changing code, tools, models, or databases. See [known limitations](docs/KNOWN_LIMITATIONS.md), especially the fragment-sequence and PROVEAN cache limitations.

## Source and attribution

- Original MutaBind2: [MinghuiLab/MutaBind2](https://github.com/minghuilab/MutaBind2), MIT license, copyright retained.
- Original paper: [MutaBind2: Predicting the Impacts of Single and Multiple Mutations on Protein–Protein Interactions](https://doi.org/10.1016/j.isci.2020.100939), iScience (2020). Please cite it when using the method.
- This repository provides GPU/workflow improvements by ZEKUN-BIT. It is not an official MinghuiLab release and does not imply upstream endorsement.
- [License](LICENSE), [third-party notices](THIRD_PARTY_NOTICES.md), and [citation metadata](CITATION.cff).

The `native/` directory contains CUDA kernels and minimal integration patches; it does not contain the complete CHARMM source or binaries. The separately built PROVEAN source in `toolchain/vendor/provean-1.1.5/` retains its GPL license and modification notice. External tools and force-field resources retain their own terms and attribution.
