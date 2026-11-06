#!/usr/bin/env python3
"""Audit and describe saved single-mutant predictions; no scientific reruns.

All thresholds operate on the saved two-decimal scores. Outputs are descriptive;
the 19 substitutions at one site are different variants, not replicates.
"""
from __future__ import annotations

import collections
import hashlib
import json
import platform
import re
import sys
import tarfile
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image

OUT = Path(__file__).resolve().parent
SCAN = OUT.parent
ROOT = SCAN.parent.parent
INPUTS = ROOT / "projects/7mkj_256_374_20261102_inputs"
SCORE = "DeltaDeltaG_kcal_mol"
AA = "ACDEFGHIKLMNPQRSTVWY"
CONTRIBUTIONS = ["DDE_vdw", "DDG_solv", "DDG_fold", "CS", "SA_com_wt", "SA_part_wt", "N_cont_wt"]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(name, data):
    (OUT / name).write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")


def write_tsv(name, frame):
    frame.to_csv(OUT / name, sep="\t", index=False, float_format="%.6f")


def md_table(frame, columns, headers):
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in frame[columns].itertuples(index=False, name=None):
        lines.append("| " + " | ".join(f"{x:.2f}" if isinstance(x, (float, np.floating)) else str(x) for x in row) + " |")
    return "\n".join(lines)


def validate():
    paths = [SCAN / x for x in ["completed_predictions.tsv", "predictions.tsv", "scan.json", "protocol.json",
             "pipeline_status.json", "plan/mutations.tsv", "plan/site_annotations.tsv", "plan/fragments.pdb", "plan/fragments.json"]]
    paths += [INPUTS / x for x in ["256_374.original.txt", "mutations.tsv", "7MKJ2610.pdb", "RpoD_full_sequence.fasta", "input_audit.json"]]
    batch_paths = sorted((SCAN / "batches").glob("batch_*/results.tsv"))
    assert len(batch_paths) == 119
    for path in batch_paths:
        paths.extend([path, path.with_name("features.tsv"), path.with_name("status.json")])
    hashes = {str(x): sha(x) for x in paths}
    assert hashes[str(SCAN / "completed_predictions.tsv")] == hashes[str(SCAN / "predictions.tsv")]
    df = pd.read_csv(SCAN / "completed_predictions.tsv", sep="\t", keep_default_na=False)
    original_columns = df.columns.tolist()
    parsed = df.Mutation.str.extract(r"^([A-Z])(\d+)([A-Z])$")
    assert parsed.notna().all().all()
    df["WT"], df["Position"], df["Mutant"] = parsed[0], parsed[1].astype(int), parsed[2]
    for col in [SCORE, "PROVEAN", "Partner_min_heavy_distance_A"]:
        df[col] = pd.to_numeric(df[col], errors="raise")
        assert np.isfinite(df[col]).all()
    for col in ["Truncated_terminus", "Near_truncated_terminus_10A"]:
        df[col] = df[col].astype(int)
        assert df[col].isin([0, 1]).all()
    assert len(df) == 2261 and df.Result_Id.is_unique
    assert not df.duplicated(["Chain", "Mutation"]).any()
    assert df.Status.eq("complete").all() and df.Chain.eq("L").all()
    assert df.Model_Confidence.eq(np.where(df.Model_Interface.eq("yes"), "High", "Low")).all()
    assert df.Model_Interface.isin(["yes", "no"]).all() and df.PROVEAN_mode.eq("full").all()
    assert df.Predicted_binding_effect.eq(np.select([df[SCORE].lt(0), df[SCORE].gt(0)], ["stronger", "weaker"], default="zero_rounded")).all()
    requests = []
    for line in (INPUTS / "256_374.original.txt").read_text().splitlines():
        if line.strip():
            chain, wtpos, mut = line.replace("\\t", "\t").split()
            requests.append((chain, wtpos + mut))
    assert len(requests) == len(set(requests)) == len(df)
    assert set(requests) == set(zip(df.Chain, df.Mutation))
    normalized = pd.read_csv(INPUTS / "mutations.tsv", sep="\t", dtype=str)
    plan = pd.read_csv(SCAN / "plan/mutations.tsv", sep="\t", dtype=str, keep_default_na=False)
    identity = {(str(x.Result_Id), x.Chain, x.Mutation) for x in df.itertuples()}
    assert identity == set(zip(normalized.id, normalized.chain, normalized.mutation))
    assert identity == set(zip(plan.result_id, plan.chain, plan.mutation))
    scan = json.loads((SCAN / "scan.json").read_text())
    protocol = json.loads((SCAN / "protocol.json").read_text())
    pipeline = json.loads((SCAN / "pipeline_status.json").read_text())
    assert identity == {tuple(x) for x in scan["mutation_identity"]}
    assert scan["partner1"] == ["I", "J"] and scan["partner2"] == ["L"] and scan["mutated_partner"] == 2
    assert pipeline["stage"] == "complete" and pipeline["protocol_sha256"] == protocol["sha256"]
    assert sha(INPUTS / "7MKJ2610.pdb") == scan["source_sha256"]
    assert sha(SCAN / "plan/fragments.pdb") == scan["plan"]["prepared_pdb_sha256"]
    audit = json.loads((INPUTS / "input_audit.json").read_text())
    assert sha(INPUTS / "256_374.original.txt") == audit["source_sha256"]
    original_user_file = Path(audit["source"])
    if original_user_file.exists():
        assert sha(original_user_file) == audit["source_sha256"]
        hashes[str(original_user_file)] = sha(original_user_file)
    seq = "".join(x.strip() for x in (INPUTS / "RpoD_full_sequence.fasta").read_text().splitlines() if not x.startswith(">"))
    annotations = pd.read_csv(SCAN / "plan/site_annotations.tsv", sep="\t").set_index("original_number")
    groups = list(df.groupby("Position", sort=True))
    assert [pos for pos, _ in groups] == list(range(256, 375))
    for pos, group in groups:
        assert len(group) == 19 and group.WT.nunique() == 1
        assert set(group.Mutant) == set(AA) - {group.WT.iloc[0]}
        assert seq[pos - 1] == group.WT.iloc[0]
        for column, annotation in [("Partner_min_heavy_distance_A", "partner_min_heavy_distance_A"),
                ("Truncated_terminus", "is_truncated_terminus"), ("Near_truncated_terminus_10A", "near_truncated_terminus_10A")]:
            assert group[column].nunique() == 1 and group[column].iloc[0] == annotations.loc[pos, annotation]

    # Independent coordinate check of native residue identities and distances.
    aa3 = dict(zip(["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"], AA))
    partner, residues, names = [], collections.defaultdict(list), {}
    for line in (INPUTS / "7MKJ2610.pdb").read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        if line[76:78].strip() in {"H", "D"} or line[12:16].strip().lstrip("0123456789").startswith(("H", "D")):
            continue
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        if line[21] in {"I", "J"}:
            partner.append(xyz)
        elif line[21] == "L":
            pos = int(line[22:26])
            residues[pos].append(xyz)
            names[pos] = aa3[line[17:20]]
    partner = np.asarray(partner)
    errors = []
    for pos, group in groups:
        assert names[pos] == group.WT.iloc[0]
        distance = np.sqrt(min(np.sum((partner - point) ** 2, axis=1).min() for point in residues[pos]))
        errors.append(abs(float(distance) - group.Partner_min_heavy_distance_A.iloc[0]))
    assert max(errors) <= 0.000051

    results, features = [], []
    for path in batch_paths:
        status = json.loads(path.with_name("status.json").read_text())
        assert status["stage"] == "complete" and status["predictions"] == 19
        assert status["protocol_sha256"] == protocol["sha256"]
        for filename in ["results.tsv", "features.tsv"]:
            assert status["saved_sha256"][filename] == sha(path.with_name(filename))
        r = pd.read_csv(path, sep="\t", keep_default_na=False)
        f = pd.read_csv(path.with_name("features.tsv"), sep="\t", keep_default_na=False)
        assert len(r) == len(f) == 19
        r["Source_batch"] = path.parent.name
        results.append(r)
        features.append(f)
    rr, ff = pd.concat(results).set_index("Result_Id"), pd.concat(features).set_index("Result_Id")
    assert rr.index.is_unique and ff.index.is_unique and set(rr.index) == set(ff.index) == set(df.Result_Id)
    saved = rr.loc[df.Result_Id]
    raw = ff.loc[df.Result_Id]
    for a, b in [(SCORE, "SunDDG"), ("Model_Interface", "Interface"), ("Model_Confidence", "Confidence")]:
        assert np.array_equal(df[a].to_numpy(), saved[b].to_numpy())
    assert np.array_equal(df.PROVEAN.to_numpy(), raw.provean.to_numpy())
    assert np.array_equal(df.Model_Interface.to_numpy(), raw.interface.to_numpy())
    assert np.array_equal(df.Mutation.to_numpy(), raw.Mutation_PDB.to_numpy())
    df["Source_batch"] = saved.Source_batch.to_numpy()
    # Raw model input features and model attribution columns are kept distinct.
    for col in ["ddVdw", "ddPb", "Foldx_fold", "ACC_wt", "ACCmp_wt", "len_mp_contact_op_wt10"]:
        df["Feature_" + col] = pd.to_numeric(raw[col]).to_numpy()
    for col in CONTRIBUTIONS:
        df["Contribution_" + col] = pd.to_numeric(saved[col]).to_numpy()
    df["Structure_caution"] = np.where(df.Near_truncated_terminus_10A.eq(1), "within_10A_of_truncated_terminus", "no_10A_truncation_flag")
    df["Protocol_sha256"] = protocol["sha256"]
    checks = {"all_2261_requested_identities_match": True, "all_119_sites_have_19_non_WT_substitutions": True,
              "all_scores_finite_and_status_complete": True, "all_WT_match_PDB_and_full_sequence": True,
              "all_predictions_and_labels_match_119_batches": True, "all_batch_result_feature_hashes_match_saved_status": True,
              "single_protocol_sha256": protocol["sha256"], "distance_max_absolute_error_A": float(max(errors)),
              "all_structural_flags_match_plan_annotations": True, "completed_predictions_identical_to_predictions": True}
    return df, original_columns, checks, hashes, pipeline


def audit_interface(df, hashes):
    records = []
    # Resolve the two sites with mutation-dependent interface labels from saved SASA logs.
    for batch in ["batch_0001", "batch_0004"]:
        path = SCAN / "batches" / batch / "scientific_intermediates.tar.gz"
        hashes[str(path)] = sha(path)
        status = json.loads(path.with_name("status.json").read_text())
        assert hashes[str(path)] == status["saved_sha256"][path.name]
        with tarfile.open(path) as archive:
            for member in archive.getmembers():
                match = re.search(r"surface_fragments_([ar])g(131|132)([a-z])_g_\d+\.out$", member.name)
                if not match or not member.isfile():
                    continue
                wt, number, mutant = match.groups()
                mutation = f"{wt.upper()}{int(number) + 241}{mutant.upper()}"
                text = archive.extractfile(member).read().decode()
                values = re.findall(r'Parameter:\s+DRASAP\s+<-\s+"([^"\n]+)"', text)
                assert len(values) == 1
                value = float(values[0])
                row = df[df.Mutation.eq(mutation)].iloc[0]
                assert row.Model_Interface == ("yes" if value > 0 else "no")
                records.append({"Mutation": mutation, "Source_batch": batch, "Archive_member": member.name,
                                "DRASAP": value, "Model_Interface": row.Model_Interface, "Model_Confidence": row.Model_Confidence})
    frame = pd.DataFrame(records).sort_values("Mutation")
    assert len(frame) == 38 and frame.Mutation.is_unique
    write_tsv("interface_label_audit.tsv", frame)
    return frame


def summarize_sites(df):
    records = []
    for pos, g in df.groupby("Position", sort=True):
        g = g.sort_values([SCORE, "Mutation"], kind="stable")
        v = g[SCORE]
        records.append({"Chain": "L", "Position": pos, "WT": g.WT.iloc[0], "Site": f"{g.WT.iloc[0]}{pos}",
            "N_substitutions": len(g), "Mean_DeltaDeltaG_kcal_mol": float(v.mean()),
            "Median_DeltaDeltaG_kcal_mol": float(v.median()), "Min_DeltaDeltaG_kcal_mol": float(v.min()),
            "Max_DeltaDeltaG_kcal_mol": float(v.max()), "N_negative": int(v.lt(0).sum()),
            "N_zero_rounded": int(v.eq(0).sum()), "N_positive": int(v.gt(0).sum()),
            "N_le_minus_0_5": int(v.le(-0.5).sum()), "N_ge_plus_1_5": int(v.ge(1.5).sum()),
            "Best_scoring_mutation": g.Mutation.iloc[0], "Worst_scoring_mutation": g.Mutation.iloc[-1],
            "N_interface_yes": int(g.Model_Interface.eq("yes").sum()), "N_confidence_High": int(g.Model_Confidence.eq("High").sum()),
            "Partner_min_heavy_distance_A": float(g.Partner_min_heavy_distance_A.iloc[0]),
            "Near_truncated_terminus_10A": int(g.Near_truncated_terminus_10A.iloc[0])})
    return pd.DataFrame(records)


def save_figure(fig, name):
    for ext in ["png", "svg"]:
        fig.savefig(OUT / f"{name}.{ext}", dpi=180, facecolor="white")
    plt.close(fig)


def plots(df, sites, bands):
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10, "svg.fonttype": "none",
                         "axes.spines.top": False, "axes.spines.right": False}):
        fig = plt.figure(figsize=(14, 9), layout="constrained")
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15])
        ax = fig.add_subplot(gs[0, 0])
        edges = np.round(np.arange(-1.0, 2.31, 0.1), 2)
        ax.hist(df[SCORE], bins=edges, color="#536C82", edgecolor="white", linewidth=0.4)
        for value, style in [(0, "-"), (1.5, "--")]:
            ax.axvline(value, c="#222222", ls=style, lw=1)
        ax.set(xlim=(-1, 2.3), xlabel=r"Predicted $\Delta\Delta G$ (kcal/mol)", ylabel="Number of substitutions",
               title="A   Distribution of all 2,261 predictions")
        ax.text(.97, .96, f"Median: {df[SCORE].median():.2f}\nRange: {df[SCORE].min():.2f} to {df[SCORE].max():.2f}",
                transform=ax.transAxes, ha="right", va="top")
        ax = fig.add_subplot(gs[0, 1])
        labels, counts = zip(*bands)
        ax.barh(range(5), counts, color=["#26609B", "#6E9EC9", "#D2D7DC", "#D99C7B", "#A64B35"], edgecolor="#333333", linewidth=.4)
        ax.set_yticks(range(5), labels, fontsize=9)
        ax.invert_yaxis()
        ax.set(xlabel="Number of substitutions", xlim=(0, max(counts) * 1.25), title="B   Descriptive score bands (kcal/mol)")
        for i, count in enumerate(counts):
            ax.text(count + 12, i, f"{count:,} ({count / len(df):.1%})", va="center", fontsize=10)
        ax = fig.add_subplot(gs[1, :])
        for flag, marker, color, label in [(0, "o", "#315C80", "No 10 A truncation flag"), (1, "s", "#9A571C", "Within 10 A of truncated terminus")]:
            s = sites[sites.Near_truncated_terminus_10A.eq(flag)]
            ax.vlines(s.Position, s.Min_DeltaDeltaG_kcal_mol, s.Max_DeltaDeltaG_kcal_mol, color=color, lw=.9)
            ax.scatter(s.Position, s.Median_DeltaDeltaG_kcal_mol, marker=marker, c=color, s=21, label=label, zorder=3)
        high = df[df.Model_Confidence.eq("High")]
        ax.scatter(high.Position, high[SCORE], s=20, marker="^", facecolors="none", edgecolors="#222222", label="Individual predictions labeled High", zorder=4)
        ax.axhline(0, c="#222222", lw=.8)
        ax.axhline(1.5, c="#444444", lw=.7, ls="--")
        lookup = sites.set_index("Position")
        for pos in [260, 275, 296, 304, 318, 323, 334, 373]:
            row = lookup.loc[pos]
            y = row.Min_DeltaDeltaG_kcal_mol if pos in [275, 304, 318, 334] else row.Max_DeltaDeltaG_kcal_mol
            ax.annotate(row.Site, (pos, y), xytext=(pos, y + (-.17 if y < 0 else .14)), ha="center", fontsize=9)
        ax.set(xlim=(254, 376), ylim=(-1.25, 2.95), xticks=list(range(260, 375, 10)), xlabel="Original chain L residue number",
               ylabel=r"Predicted $\Delta\Delta G$ (kcal/mol)", title="C   Site medians and min-max spans across 19 different substitutions (not uncertainty)")
        ax.legend(loc="upper center", ncol=3, fontsize=9, frameon=False)
        fig.suptitle("7MKJ: chain L positions 256-374 against partners I + J\nNegative = predicted stronger binding; positive = predicted weaker binding", fontsize=14)
        save_figure(fig, "overview")

        fig = plt.figure(figsize=(17, 14.5))
        gs = fig.add_gridspec(3, 1, left=.06, right=.9, bottom=.09, top=.91, hspace=.36)
        limit = float(df[SCORE].abs().max())
        norm = Normalize(vmin=-limit, vmax=limit)
        cmap = matplotlib.colormaps["RdBu_r"].with_extremes(bad="#D9D9D9")
        for i, (left, right) in enumerate([(256, 295), (296, 335), (336, 374)]):
            ax = fig.add_subplot(gs[i])
            positions = list(range(left, right + 1))
            subset = df[df.Position.between(left, right)]
            values = subset.pivot(index="Mutant", columns="Position", values=SCORE).reindex(index=list(AA), columns=positions)
            assert int(values.isna().sum().sum()) == len(positions)
            img = ax.imshow(values.to_numpy(), cmap=cmap, norm=norm, interpolation="nearest", aspect="auto")
            labels = []
            for j, pos in enumerate(positions):
                site = lookup.loc[pos]
                wtrow = AA.index(site.WT)
                ax.text(j, wtrow, "x", ha="center", va="center", color="#333333", fontsize=7)
                labels.append(site.Site + ("*" if site.Near_truncated_terminus_10A else ""))
            for row in subset.itertuples():
                x, y = row.Position - left, AA.index(row.Mutant)
                score = getattr(row, SCORE)
                if row.Model_Confidence == "High":
                    ax.plot(x + .3, y - .27, marker="o", color="black", ms=2)
                if score <= -.5 or score >= 1.5:
                    ax.text(x, y, f"{score:+.2f}", ha="center", va="center", fontsize=6.2, color="white" if abs(score) > 1.3 else "#111111")
            ax.set_xticks(range(len(positions)), labels, rotation=90, fontsize=8)
            ax.set_yticks(range(20), list(AA), fontsize=9)
            ax.set_xticks(np.arange(-.5, len(positions), 1), minor=True)
            ax.set_yticks(np.arange(-.5, 20, 1), minor=True)
            ax.grid(which="minor", color="#EEEEEE", lw=.3)
            ax.tick_params(which="minor", length=0)
            ax.set(ylabel="Mutant amino acid", title=f"{chr(65+i)}   Positions {left}-{right}")
        cax = fig.add_axes([.925, .27, .016, .45])
        bar = fig.colorbar(img, cax=cax)
        bar.set_label(r"Predicted $\Delta\Delta G$ (kcal/mol)")
        fig.suptitle("Single-mutant predictions: all 119 sites x 19 substitutions\nShared symmetric color scale centered at zero; no smoothing or clipping", fontsize=16, y=.98)
        fig.text(.06, .025, "Gray x = WT identity (not calculated).   * = site within 10 A of a truncated terminus.   Black dot = model confidence High.\nNumbers mark scores <= -0.50 or >= +1.50 only; these are screening bands, not uncertainty bounds. Full values: all_results_with_features.tsv.", fontsize=11)
        save_figure(fig, "mutation_heatmap")
    return edges.tolist(), [-limit, limit]


def report(df, sites, summary):
    enhancing = df[df[SCORE].le(-.5)].sort_values([SCORE, "Mutation"])
    weakening = df[df[SCORE].ge(1.5)].sort_values([SCORE, "Mutation"], ascending=[False, True])
    text = f"""# 7MKJ L 链 256–374 位点突变预测分析

分析对象为 L 链与 I+J 两条链的结合，全部为单点突变。119 个连续位点、每位点 19 种非野生型替换，共 2261 条结果，均完整且与原始请求逐一对应；逐批分数、标签、结果及特征表哈希、野生型序列核对通过。任务完成时间：2026-11-04 01:04:34（北京时间）。

## 总体结果

ΔΔG = 突变体结合自由能 − 野生型结合自由能，单位 kcal/mol。负值预测结合增强，正值预测结合减弱。这里的分数是将物理计算、FoldX、PROVEAN 及结构特征输入 MutaBind2 模型后的预测，不等于单独的 PB 能量，也不是所有能量项的直接加和。

- 预测增强：{summary['negative']} 条（{summary['negative']/len(df):.2%}）；预测减弱：{summary['positive']} 条（{summary['positive']/len(df):.2%}）；保存精度下为零：{summary['zero_rounded']} 条。
- 均值 {df[SCORE].mean():.3f}，中位数 {df[SCORE].median():.2f}，范围 {df[SCORE].min():.2f} 至 +{df[SCORE].max():.2f}。
- ΔΔG ≤ −0.50 有 {len(enhancing)} 条；≤ −1.50 为 0；≥ +1.50 有 {len(weakening)} 条。
- ±0.50 只是本报告整理候选所用的描述性阈值；±1.50 是项目代码用于幅度分类的界限。二者都不是本批结果的误差界、统计显著性或实验成功率。大量小正值不能据此全部称为确定有害。

## 增强方向：优先探索 V275，其次跨位点候选

{md_table(enhancing, ['Mutation', SCORE, 'Model_Confidence', 'Partner_min_heavy_distance_A'], ['突变', 'ΔΔG (kcal/mol)', '模型标签', '野生型至 I/J 最近重原子距离 (Å)'])}

V275 的 19 种替换中 11 种为负值，位点中位数 −0.40；其中 T、K、Q、R、S 五种替换 ≤ −0.50。相比仅有一个极端分数，它提供了多个同方向的候选，但这些替换互相关联，不能当作独立重复实验。按分数优先考虑 V275T、V275K；若希望覆盖不同位点，可增加 A318K、S334Q、T304K。所有这 11 条均为 Low/非界面且没有 10 Å 截断邻近标记；没有找到同时满足 High 且较大幅度增强的候选。

V275T 的原始模型输入 ddVdw 约 −0.011、ddPb 约 −0.022、FoldX_fold 约 −0.917；其最终 −0.89 分数不能解释为 PB 作用大幅增强。A318K/S334Q 的 ddPb 为负，也不能仅凭这一项推断机制，因为模型同时使用其他特征。PROVEAN 本身已参与评分，不能把它作为独立验证证据。

## 减弱方向：区分界面候选和结构影响

{md_table(weakening, ['Mutation', SCORE, 'Model_Confidence', 'Near_truncated_terminus_10A'], ['突变', 'ΔΔG (kcal/mol)', '模型标签', '距截断端点 ≤10 Å'])}

R373P（+1.80）是本批唯一同时为 High 且达到 +1.50 的突变，值得优先作为结合减弱候选。R373 所有 19 种替换的分数都为正，其中其余 18 种为 +0.21～+0.64。R373P 的 FoldX_fold 约 +2.96，提示同时存在折叠稳定性改变的可能；验证结合时应同时检查蛋白稳定性/可溶性，避免把结构受损全部归因于直接界面接触变化。

N323P（+2.23）、K296P（+2.04）、K296I（+1.94）幅度更大，但均为 Low。N323P/K296P 的 FoldX_fold 分别约 +4.40/+4.60；这些结果更适合作为可能影响结构与结合的探索候选，不能直接认定为界面热点。R260 的五个 ≥+1.50 结果均靠近截断端点，优先级应更谨慎。

如设计一个小规模验证组合，可选 WT、V275T、V275K、A318K、S334Q、T304K、R373P，并加入 R373K（+0.21）作为同位点较小效应的比较对象。它仍是计算筛选方案，具体选择取决于实验目的；单突变分数不能相加来预测组合突变。

## High / Low 和结构边界如何解释

全部 19 个 High 对应 19 个界面标记：R373 的 18 种替换加 A372W；R373G 为 Low。原始逐批表与存档中的表面积输出一致。界面判定采用每个突变结构的 DRASAP > 0：R373G 为 0，A372W 为 0.006917212，R373P 为 0.00381787；同一位点因此可出现不同标签。该规则对是否严格大于零进行分类，不能把 High 理解为校准后的概率、低误差保证或效应一定大。

171 条结果（9 个位点）距某个结构截断端点 ≤10 Å：F256、K257、Q258、F259、R260、L261、G352、T354、Q357；没有突变恰好位于截断端点。该标记考虑双方结构中的截断端点，不只考虑 L 链。没有标记也不能保证完全不受缺失结构影响。自动片段化和中性封端处理仍不能还原缺失片段；全序列 PROVEAN 也不会补回三维结构。

结果针对当前结构、I+J/L 伙伴划分和现有模型，不能直接推断整套复合物的功能、转录活性、表达量或生物学适应度。当前没有本体系的实验标定或独立重复计算，不能给出本批特异的置信区间，也不应按两位小数的微小差别过度排序。

## 文件与复现

- 原始完整结果：[completed_predictions.tsv](../completed_predictions.tsv)。
- 全部增强/减弱排序：[ranked_enhancing.tsv](ranked_enhancing.tsv)、[ranked_weakening.tsv](ranked_weakening.tsv)。
- 23 条幅度筛选候选：[screening_candidates.tsv](screening_candidates.tsv)；每位点的 19 种替换汇总：[site_summary.tsv](site_summary.tsv)。
- 汇总及原始模型输入：[all_results_with_features.tsv](all_results_with_features.tsv)。Feature_ 前缀为模型输入，Contribution_ 前缀为模型归因列，归因值不是独立物理能量。
- 图：[overview.png](overview.png)、[mutation_heatmap.png](mutation_heatmap.png)，另有 SVG；[图注与文字说明](figure_captions_zh.md)。
- [机器可读汇总](summary.json)、[审计及来源记录](analysis_manifest.json)、[界面标签原始输出核对](interface_label_audit.tsv)。
- [分析脚本](analyze_results.py)只读取原始数据、向本分析目录写衍生文件；本次未重新运行任何物理计算。分析依赖版本记录在 requirements-analysis.txt，解释器记录在 analysis_manifest.json。

本次分析仅针对 256–374 批次，未与此前不同来源/基线的 128–255 批次合并比较。
"""
    text = re.sub(r"\]\(([^)]+)\)", lambda m: "](" + str((OUT / m.group(1)).resolve()) + ")", text)
    (OUT / "report_zh.md").write_text(text)


def main():
    df, original_columns, checks, hashes, pipeline = validate()
    audit_interface(df, hashes)
    sites = summarize_sites(df)
    enhancing = df[df[SCORE].lt(0)].sort_values([SCORE, "Mutation"], kind="stable")
    weakening = df[df[SCORE].gt(0)].sort_values([SCORE, "Mutation"], ascending=[False, True], kind="stable")
    candidates = df[df[SCORE].le(-.5) | df[SCORE].ge(1.5)].copy()
    candidates["Screening_group"] = np.where(candidates[SCORE].lt(0), "exploratory_enhancement", "predicted_weakening")
    candidates = candidates.sort_values(["Screening_group", SCORE, "Mutation"])
    for name, frame in [("ranked_enhancing.tsv", enhancing), ("ranked_weakening.tsv", weakening),
                        ("screening_candidates.tsv", candidates), ("site_summary.tsv", sites), ("all_results_with_features.tsv", df)]:
        write_tsv(name, frame)
    v = df[SCORE]
    bands = [("ddG <= -1.5", int(v.le(-1.5).sum())), ("-1.5 < ddG <= -0.5", int((v.gt(-1.5) & v.le(-.5)).sum())),
             ("-0.5 < ddG < +0.5", int((v.gt(-.5) & v.lt(.5)).sum())), ("+0.5 <= ddG < +1.5", int((v.ge(.5) & v.lt(1.5)).sum())),
             ("ddG >= +1.5", int(v.ge(1.5).sum()))]
    assert sum(count for _, count in bands) == len(df)
    summary = {"n_predictions": len(df), "n_sites": len(sites), "substitutions_per_site": 19,
        "position_range": [256, 374], "positive": int(v.gt(0).sum()), "negative": int(v.lt(0).sum()),
        "zero_rounded": int(v.eq(0).sum()), "mean": float(v.mean()), "median": float(v.median()),
        "min": float(v.min()), "max": float(v.max()), "descriptive_score_bands": dict(bands),
        "confidence_counts": df.Model_Confidence.value_counts().to_dict(), "interface_counts": df.Model_Interface.value_counts().to_dict(),
        "near_truncated_terminus_predictions": int(df.Near_truncated_terminus_10A.sum()),
        "near_truncated_terminus_sites": sites.loc[sites.Near_truncated_terminus_10A.eq(1), "Site"].tolist(),
        "best_enhancement": enhancing.iloc[0][original_columns].to_dict(), "strongest_weakening": weakening.iloc[0][original_columns].to_dict(),
        "pipeline_completed_utc": pipeline["updated_utc"]}
    write_json("summary.json", summary)
    bins, colorlimits = plots(df, sites, bands)
    report(df, sites, summary)
    (OUT / "figure_captions_zh.md").write_text("""# 图注与文字说明

## overview

7MKJ 中 L 链 256–374 位点对 I+J 结合的 2261 个单突变预测。A：保存至两位小数的 ΔΔG 频数分布，直方图宽度 0.1 kcal/mol；全部数据均纳入。B：描述性分数分档，±0.5 仅用于整理候选，不代表误差范围；±1.5 为代码的幅度分类界限。C：每位点 19 种不同替换的中位数和最小值至最大值范围，线段不是误差条或置信区间；方块指距截断端点 ≤10 Å 的位点，空心三角标出单条 High 预测。负值预测增强，正值预测减弱。

文字替代：分数主要集中在小的正值；62 个为负，2195 个为正，4 个保存为零。V275 的 11/19 替换为负，最小 −0.89；最强减弱 N323P 为 +2.23。R373P 为 +1.80，是唯一同时达到 +1.5 和 High 的预测。不存在 ≤−1.5 的预测。

## mutation_heatmap

119 列表示原始 L 链连续位点，20 行表示突变后的氨基酸，拆成三幅分段热图。灰色 x 表示野生型自身，未计算，不能当作零。三幅共用以零为中心的 −2.23 至 +2.23 kcal/mol 对称线性色标；不平滑、不截断数值。星号表示位点距截断端点 ≤10 Å；黑点表示该突变被模型标为 High；只有 ≤−0.5 或 ≥+1.5 的格子显示数值。全部数值可查 all_results_with_features.tsv。

文字替代：大部分格子为正值；V275 有多种负分替换。N323P、K296P/K296I、R373P 及 R260 的若干替换是较大的正值。High 标签集中在 R373（G 替换除外）以及 A372W。色彩提供总体浏览，精确筛选应使用表格及模型/结构标签。

这两幅图是本项目的通用分析图，未以特定期刊规范设计或认证。数据为不同突变的模型预测，没有实验重复或可据此绘制的统计置信区间。
""")
    packages = ["numpy", "pandas", "matplotlib", "pillow", "contourpy", "cycler", "fonttools", "kiwisolver", "packaging", "pyparsing", "python-dateutil", "pytz", "tzdata"]
    versions = {name: version(name) for name in packages}
    (OUT / "requirements-analysis.txt").write_text("\n".join(f"{name}=={value}" for name, value in versions.items()) + "\n")
    assert all(sha(path) == before for path, before in hashes.items()), "An input file changed during analysis"
    checks["all_source_hashes_unchanged_after_analysis"] = True
    figure_metadata = {}
    for filename in ["overview.png", "mutation_heatmap.png"]:
        with Image.open(OUT / filename) as image:
            figure_metadata[filename] = {"size_pixels": list(image.size), "mode": image.mode, "dpi": image.info.get("dpi"), "bytes": (OUT / filename).stat().st_size}
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(), "analysis_script_sha256": sha(Path(__file__)),
        "python_executable": sys.executable, "python_version": platform.python_version(), "libraries": versions,
        "scope": "descriptive analysis of one saved single-mutant scan; no new physical calculations",
        "validation": checks, "source_sha256": hashes, "histogram_bin_edges_kcal_mol": bins,
        "heatmap_normalization": {"kind": "linear, symmetric, zero-centered", "limits_kcal_mol": colorlimits, "missing_WT_cells": 119, "interpolation": "nearest"},
        "uncertainty": "none estimated; 19 substitutions per site are different variants, not replicates",
        "figure_metadata": figure_metadata,
        "outputs_sha256": {x.name: sha(x) for x in sorted(OUT.iterdir()) if x.is_file() and x.name != "analysis_manifest.json"}}
    write_json("analysis_manifest.json", manifest)
    print(json.dumps({"summary": summary, "validation": checks, "output_directory": str(OUT)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
