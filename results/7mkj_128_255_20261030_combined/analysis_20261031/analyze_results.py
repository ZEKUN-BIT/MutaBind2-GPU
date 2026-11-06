#!/usr/bin/env python3
"""Descriptive analysis of the finished scan; never reruns or edits predictions.

Run with the bundled Python runtime (numpy, pandas, matplotlib).
All effect categories use the rounded, two-decimal prediction scores as saved.
"""
from __future__ import annotations

import collections
import hashlib
import json
import math
import platform
import re
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
COMBINED = OUT.parent
ROOT = COMBINED.parent.parent
INPUTS = ROOT / "projects/7mkj_128_255_20261030_inputs"
RESULTS = COMBINED / "completed_predictions.tsv"
INPUT = INPUTS / "128-255.original.txt"
PDB = INPUTS / "7MKJ2610.pdb"
AA = "ACDEFGHIKLMNPQRSTVWY"
RANGES = [(128, 166), (213, 236), (242, 255)]
SCORE = "DeltaDeltaG_kcal_mol"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(name, data):
    (OUT / name).write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")


def write_tsv(name, frame):
    frame = frame.copy()
    formats = {SCORE: ".2f", "PROVEAN": ".3f", "Partner_min_heavy_distance_A": ".4f",
               "Median_DeltaDeltaG_kcal_mol": ".2f", "Min_DeltaDeltaG_kcal_mol": ".2f",
               "Max_DeltaDeltaG_kcal_mol": ".2f"}
    for col, fmt in formats.items():
        if col in frame:
            frame[col] = frame[col].map(lambda value: format(value, fmt))
    frame.to_csv(OUT / name, sep="\t", index=False, float_format="%.6f")


def md_table(frame, columns, headers):
    text = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for _, row in frame.iterrows():
        vals = []
        for col in columns:
            v = row[col]
            if isinstance(v, (float, np.floating)):
                v = f"{v:.2f}"
            vals.append(str(v))
        text.append("| " + " | ".join(vals) + " |")
    return "\n".join(text)


def load_and_validate():
    original_hash = sha(RESULTS)
    df = pd.read_csv(RESULTS, sep="\t", keep_default_na=False)
    original_columns = df.columns.tolist()
    parsed = df.Mutation.str.extract(r"^([A-Z])(\d+)([A-Z])$")
    assert not parsed.isna().any().any(), "Unexpected mutation notation"
    df["WT"] = parsed[0]
    df["Position"] = parsed[1].astype(int)
    df["Mutant"] = parsed[2]
    for col in [SCORE, "PROVEAN", "Partner_min_heavy_distance_A"]:
        df[col] = pd.to_numeric(df[col], errors="raise")
        assert np.isfinite(df[col]).all(), col
    for col in ["Truncated_terminus", "Near_truncated_terminus_10A"]:
        df[col] = df[col].astype(int)
        assert df[col].isin([0, 1]).all()
    df["Source_scan_name"] = df.Source_scan.map(lambda x: Path(x).name)
    requested = []
    for line in INPUT.read_text().splitlines():
        if line.strip():
            chain, wt_position, mutant = line.split()
            requested.append((chain, wt_position + mutant))
    found = list(zip(df.Chain, df.Mutation))
    assert len(requested) == len(set(requested)) == len(found) == len(set(found))
    assert set(requested) == set(found)
    assert df.Result_Id.is_unique and df.Status.eq("complete").all()
    assert df.Chain.eq("L").all()
    assert df.Predicted_binding_effect.tolist() == [
        "weaker" if v > 0 else "stronger" if v < 0 else "zero_rounded"
        for v in df[SCORE]
    ]
    groups = list(df.groupby("Position", sort=True))
    for _, group in groups:
        assert len(group) == 19 and group.WT.nunique() == 1
        assert set(group.Mutant) == set(AA) - set(group.WT.iloc[0])
        assert group.Model_Interface.nunique() == group.Model_Confidence.nunique() == 1
        for col in ["Truncated_terminus", "Near_truncated_terminus_10A", "Partner_min_heavy_distance_A"]:
            assert group[col].nunique() == 1, col

    # Independently check distances and WT identity against original coordinates.
    aa3 = dict(zip(
        ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"], AA))
    partner = []
    chain_l = collections.defaultdict(list)
    residue_names = {}
    for line in PDB.read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        if line[76:78].strip() in {"H", "D"} or line[12:16].strip().lstrip("0123456789").startswith(("H", "D")):
            continue
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        if line[21] in {"I", "J"}:
            partner.append(xyz)
        elif line[21] == "L":
            pos = int(line[22:26])
            chain_l[pos].append(xyz)
            residue_names[pos] = aa3[line[17:20]]
    partner = np.asarray(partner)
    distance_errors = []
    for pos, group in groups:
        assert residue_names[pos] == group.WT.iloc[0]
        points = np.asarray(chain_l[pos])
        distance = np.sqrt(np.min(np.sum((points[:, None, :] - partner[None, :, :]) ** 2, axis=2)))
        distance_errors.append(abs(float(distance) - group.Partner_min_heavy_distance_A.iloc[0]))
    assert max(distance_errors) <= 0.000051, "Distance annotations differ from coordinates"
    missing = sorted(set(range(128, 256)) - set(df.Position))
    assert missing == [pos for pos in range(128, 256) if pos not in chain_l]

    # Verify each merged score and label against its actual source batch.
    source_counts = {}
    for source, group in df.groupby("Source_scan"):
        source_path = Path(source)
        source_counts[source_path.name] = len(group)
        manifest = json.loads((source_path / "scan.json").read_text())
        detail = pd.read_csv(source_path / "plan/mutations.tsv", sep="\t", dtype=str)
        assert detail.result_id.is_unique
        mapping = detail.set_index("result_id")
        records = []
        for batch_file in sorted((source_path / "batches").glob("batch_*/results.tsv")):
            batch_rows = pd.read_csv(batch_file, sep="\t", dtype=str)
            records.extend(batch_rows.to_dict("records"))
        batches = {}
        for row in records:
            identity = (mapping.loc[row["Result_Id"], "chain"], mapping.loc[row["Result_Id"], "mutation"])
            assert identity not in batches
            batches[identity] = row
        for row in group.to_dict("records"):
            saved = batches[(row["Chain"], row["Mutation"])]
            assert float(saved["SunDDG"]) == row[SCORE]
            assert saved["Interface"] == row["Model_Interface"]
            assert saved["Confidence"] == row["Model_Confidence"]
        assert manifest["partner1"] == ["I", "J"] and manifest["partner2"] == ["L"]
    df["Structure_caution"] = np.select(
        [df.Truncated_terminus.eq(1), df.Near_truncated_terminus_10A.eq(1)],
        ["at_truncated_terminus", "within_10A_of_truncated_terminus"],
        default="no_10A_truncation_flag",
    )
    return df, original_columns, {
        "input_count": len(requested), "result_count": len(df),
        "exact_input_identity_match": True, "duplicate_mutations": 0,
        "duplicate_result_ids": 0, "all_complete": True,
        "finite_scores": True, "all_WT_identities_match_coordinates": True,
        "each_site_has_all_19_non_WT_substitutions": True,
        "merged_scores_labels_match_source_batches": True,
        "distance_max_absolute_recheck_error_A": float(max(distance_errors)),
        "missing_positions_match_missing_coordinates": True,
        "missing_positions_within_128_255": missing,
        "source_sha256_before": original_hash,
        "source_counts": source_counts,
    }


def summarize_sites(df):
    records = []
    for pos, g in df.groupby("Position", sort=True):
        g = g.sort_values([SCORE, "Mutation"], kind="stable")
        v = g[SCORE]
        records.append({
            "Chain": g.Chain.iloc[0], "Position": pos, "WT": g.WT.iloc[0],
            "Site": f"{g.WT.iloc[0]}{pos}", "N_substitutions": len(g),
            "Mean_DeltaDeltaG_kcal_mol": float(v.mean()),
            "Median_DeltaDeltaG_kcal_mol": float(v.median()),
            "Min_DeltaDeltaG_kcal_mol": float(v.min()),
            "Max_DeltaDeltaG_kcal_mol": float(v.max()),
            "N_negative": int(v.lt(0).sum()), "N_zero_rounded": int(v.eq(0).sum()),
            "N_positive": int(v.gt(0).sum()),
            "N_le_minus_0_5": int(v.le(-0.5).sum()),
            "N_ge_plus_0_5": int(v.ge(0.5).sum()),
            "N_ge_plus_1_5": int(v.ge(1.5).sum()),
            "Best_scoring_mutation": g.Mutation.iloc[0],
            "Worst_scoring_mutation": g.Mutation.iloc[-1],
            "Partner_min_heavy_distance_A": float(g.Partner_min_heavy_distance_A.iloc[0]),
            "Truncated_terminus": int(g.Truncated_terminus.iloc[0]),
            "Near_truncated_terminus_10A": int(g.Near_truncated_terminus_10A.iloc[0]),
            "Structure_caution": g.Structure_caution.iloc[0],
            "Model_Interface": g.Model_Interface.iloc[0],
            "Model_Confidence": g.Model_Confidence.iloc[0],
            "Source_scan_name": g.Source_scan_name.iloc[0],
            "Protocol_sha256": g.Protocol_sha256.iloc[0],
        })
    return pd.DataFrame(records)


def band_counts(values):
    # +/-0.5 is a descriptive screening band, not a validated confidence limit.
    return [
        ("ddG <= -1.5", int(values.le(-1.5).sum())),
        ("-1.5 < ddG <= -0.5", int((values.gt(-1.5) & values.le(-0.5)).sum())),
        ("-0.5 < ddG < +0.5", int((values.gt(-0.5) & values.lt(0.5)).sum())),
        ("+0.5 <= ddG < +1.5", int((values.ge(0.5) & values.lt(1.5)).sum())),
        ("ddG >= +1.5", int(values.ge(1.5).sum())),
    ]


def save_figure(fig, stem):
    for ext in ["png", "svg"]:
        fig.savefig(OUT / f"{stem}.{ext}", dpi=180, facecolor="white")
    plt.close(fig)


def make_overview(df, sites, bands):
    fig = plt.figure(figsize=(14, 8.4), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1])
    ax = fig.add_subplot(gs[0, 0])
    edges = np.round(np.arange(-1.5, 2.31, 0.1), 2)
    ax.hist(df[SCORE], bins=edges, color="#5C7389", edgecolor="white", linewidth=0.5)
    for value, style in [(0, "-"), (-1.5, "--"), (1.5, "--")]:
        ax.axvline(value, color="#222222", ls=style, lw=1)
    ax.set(xlabel=r"Predicted $\Delta\Delta G$ (kcal/mol)", ylabel="Number of substitutions",
           title="A   Distribution of 1,463 single-mutant predictions", xlim=(-1.6, 2.3))
    ax.text(.98, .95, f"Median = {df[SCORE].median():.2f}\nRange = {df[SCORE].min():.2f} to {df[SCORE].max():.2f}",
            ha="right", va="top", transform=ax.transAxes, fontsize=10)
    ax = fig.add_subplot(gs[0, 1])
    labels, counts = zip(*bands)
    colors = ["#26609B", "#6E9EC9", "#D2D7DC", "#D99C7B", "#A64B35"]
    ax.barh(range(5), counts, color=colors, edgecolor="#333333", linewidth=.4)
    ax.set_yticks(range(5), labels, fontsize=9)
    ax.invert_yaxis()
    ax.set(xlabel="Number of substitutions", xlim=(0, max(counts) * 1.22),
           title="B   Descriptive score bands (kcal/mol)")
    for i, count in enumerate(counts):
        ax.text(count + 12, i, f"{count:,} ({count / len(df):.1%})", va="center", fontsize=10)
    ax = fig.add_subplot(gs[1, :])
    for left, right in [(167, 212), (237, 241)]:
        ax.axvspan(left - .5, right + .5, facecolor="#EFEFEF", hatch="///", edgecolor="#CCCCCC", lw=0)
    for flag, marker, color, label in [
        (0, "o", "#315C80", "No 10 A truncation flag"),
        (1, "s", "#BC6D25", "Within 10 A of a truncated terminus"),
    ]:
        s = sites[sites.Near_truncated_terminus_10A.eq(flag)]
        y = s.Median_DeltaDeltaG_kcal_mol
        ax.errorbar(s.Position, y,
                    yerr=np.array([y - s.Min_DeltaDeltaG_kcal_mol, s.Max_DeltaDeltaG_kcal_mol - y]),
                    fmt=marker, color=color, markersize=4.2, capsize=2, lw=.8, label=label)
    ax.axhline(0, c="#222222", lw=.8)
    ax.axhline(1.5, c="#444444", lw=.7, ls="--")
    lookup = sites.set_index("Position")
    for pos, dy in [(152, .25), (154, .26), (156, -.27), (218, .35), (226, -.3)]:
        row = lookup.loc[pos]
        ax.annotate(row.Site, (pos, row.Median_DeltaDeltaG_kcal_mol),
                    xytext=(pos, row.Median_DeltaDeltaG_kcal_mol + dy), ha="center", fontsize=9,
                    arrowprops={"arrowstyle": "-", "color": "#555555", "lw": .5})
    ax.set(xlim=(126, 257), ylim=(-1.5, 2.5), xticks=list(range(130, 256, 10)),
           xlabel="Original chain L residue number", ylabel=r"Predicted $\Delta\Delta G$ (kcal/mol)",
           title="C   Site medians and min-max spans across 19 substitutions per site")
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    fig.suptitle("7MKJ: chain L mutations against partners I + J\nAll predictions: non-interface / model confidence Low", fontsize=14)
    save_figure(fig, "overview")
    return edges.tolist()


def make_heatmap(df, sites):
    fig = plt.figure(figsize=(19, 7.1))
    gs = fig.add_gridspec(2, 3, left=.045, right=.935, bottom=.22, top=.81,
                          height_ratios=[.35, 5], width_ratios=[39, 24, 14], hspace=.05, wspace=.10)
    limit = float(df[SCORE].abs().max())
    norm = Normalize(vmin=-limit, vmax=limit)
    cmap = matplotlib.colormaps["RdBu_r"].with_extremes(bad="#E5E5E5")
    lookup = sites.set_index("Position")
    for j, (left, right) in enumerate(RANGES):
        positions = list(range(left, right + 1))
        matrix = np.full((20, len(positions)), np.nan)
        for col, pos in enumerate(positions):
            for row in df[df.Position.eq(pos)].itertuples(index=False):
                matrix[AA.index(row.Mutant), col] = getattr(row, SCORE)
        ax_top = fig.add_subplot(gs[0, j])
        for col, pos in enumerate(positions):
            s = lookup.loc[pos]
            color = "#333333" if s.Truncated_terminus else "#D9A268" if s.Near_truncated_terminus_10A else "#FFFFFF"
            ax_top.add_patch(Rectangle((col - .5, -.5), 1, 1, facecolor=color, edgecolor="#999999", lw=.3,
                                      hatch="///" if s.Near_truncated_terminus_10A and not s.Truncated_terminus else None))
        ax_top.set(xlim=(-.5, len(positions) - .5), ylim=(-.5, .5), xticks=[], yticks=[], title=f"Residues {left}-{right}")
        for spine in ax_top.spines.values(): spine.set_visible(False)
        ax = fig.add_subplot(gs[1, j])
        im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        for col, pos in enumerate(positions):
            ax.plot(col, AA.index(lookup.loc[pos, "WT"]), marker="x", ms=3.5, mew=.65, color="#444444")
        ax.set_xticks(range(len(positions)), [lookup.loc[pos, "Site"] for pos in positions], rotation=90, fontsize=8)
        ax.set_yticks(range(20), list(AA) if j == 0 else [], fontsize=10)
        ax.tick_params(length=0)
        if j == 0: ax.set_ylabel("Mutant amino acid")
    cax = fig.add_axes([.95, .24, .012, .50])
    fig.colorbar(im, cax=cax, label=r"Predicted $\Delta\Delta G$ (kcal/mol)")
    fig.suptitle("7MKJ chain L: all 1,463 single substitutions\nNegative = predicted stronger binding; positive = predicted weaker binding", y=.965, fontsize=15)
    handles = [Patch(facecolor="white", edgecolor="#888888", label="No 10 A truncation flag"),
               Patch(facecolor="#D9A268", edgecolor="#888888", hatch="///", label="Within 10 A of truncated terminus"),
               Patch(facecolor="#333333", label="At truncated terminus"),
               Line2D([], [], color="#444444", marker="x", ls="none", label="WT identity: not predicted")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .075), ncols=4, frameon=False, fontsize=10)
    fig.text(.5, .025, "All model confidence labels are Low. Positions 167-212 and 237-241 are absent from both coordinates and the mutation input.\nSymmetric linear color scale; no clipping or smoothing. Truncation flags describe structural context, not calibrated uncertainty.",
             ha="center", fontsize=10, color="#333333")
    save_figure(fig, "mutation_heatmap")
    return {"normalization": "symmetric linear around zero", "vmin": -limit, "vmax": limit,
            "WT_cells": "not computed; grey with x", "omitted_number_ranges": [[167, 212], [237, 241]],
            "interpolation": "nearest; no smoothing"}


def make_report(df, sites, summary):
    n = len(df)
    safe = df.Near_truncated_terminus_10A.eq(0)
    best_safe = df[safe & df[SCORE].le(-.5)].sort_values([SCORE, "Mutation"])
    best_all = df[df[SCORE].lt(0)].sort_values([SCORE, "Mutation"]).head(12)
    worst = df[df[SCORE].ge(1.5)].sort_values([SCORE, "Mutation"], ascending=[False, True])
    hot = sites.sort_values(["Mean_DeltaDeltaG_kcal_mol", "Position"], ascending=[False, True]).head(10)
    cols = ["Mutation", SCORE, "Partner_min_heavy_distance_A", "Truncated_terminus", "Near_truncated_terminus_10A"]
    headers = ["突变（链 L）", "ΔΔG (kcal/mol)", "距伙伴 (Å)", "位于截断端", "距截断端≤10 Å"]
    band_labels = ["≤ −1.5", "−1.5 < ΔΔG ≤ −0.5", "−0.5 < ΔΔG < 0.5", "0.5 ≤ ΔΔG < 1.5", "≥ 1.5"]
    band_lines = [f"| {label} | {count} | {count/n:.2%} |" for label, (_, count) in zip(band_labels, summary["score_bands"])]
    text = f"""# 7MKJ 链 L 的突变预测结果分析

本批次已完成 **{n}/{n}** 条预测。输入与结果逐条对应，无重复、无漏项、无无效分数，合并值与各来源批次的结果一致。完成时间：2026-10-31 10:03:56（北京时间）。

**核心结论：A156V、A156L、A226K、A226R 可作为增强结合的探索候选；E154P、E155P、E154I 可作为减弱结合的探索候选。全部预测均被模型标为非界面、Low，当前不能认定已找到可靠的亲和力增强突变。**

## 1. 计算对象与适用范围

- 结构：本地 7MKJ2610.pdb；伙伴：I+J 对 L；全部突变在 L 链。本地全长序列文件标注为 RpoD。结论针对该复合物构象和这一伙伴定义。
- 实际覆盖 77 个位点，每个位点 19 种非野生型替换：128–166（39 位点）、213–236（24 位点）、242–255（14 位点）。
- 167–212 和 237–241 共 51 个位置在当前 L 链坐标和输入表中均缺失。这次已完成全部已提交请求，并没有计算这些缺失位置。
- ΔΔG = 突变体结合自由能 − 野生型结合自由能，单位 kcal/mol。正值预测结合减弱，负值预测结合增强。这里是 MutaBind2 的最终模型预测，PB、范德华、FoldX、PROVEAN 等是输入特征，不能把某个输入能量直接当成最终 ΔΔG。
- PROVEAN_mode 全部为 full，但这只表示保守性特征使用了完整序列，**没有补全结构缺失段**。

## 2. 总体分布

| 指标 | 结果 |
| --- | --- |
| ΔΔG > 0，方向上预测减弱 | {summary['positive']}（{summary['positive']/n:.2%}） |
| ΔΔG < 0，方向上预测增强 | {summary['negative']}（{summary['negative']/n:.2%}） |
| 显示为 0.00 | {summary['zero_rounded']} |
| 均值 / 中位数 | {summary['mean']:.3f} / {summary['median']:.2f} kcal/mol |
| 最小 / 最大 | {summary['min']:.2f} / {summary['max']:.2f} kcal/mol |

按幅度分层：

| ΔΔG 区间 (kcal/mol) | 数量 | 占比 |
| --- | --- | --- |
{chr(10).join(band_lines)}

原程序以 |ΔΔG| ≥ 1.5 标记较大幅度变化；本批次达到该幅度的 23 条均为正值。±0.5 在本报告中只用于方便筛选和展示，**不代表经过验证的中性边界、误差范围或统计显著性阈值**。有 {summary['negative']} 条负值不等于它们都具有可靠增强作用，尤其是接近零的预测。

## 3. 增强结合的探索候选

优先查看未触发“距截断端≤10 Å”标记、且 ΔΔG ≤ −0.5 的 9 条。它们仍全部是模型 Low。此筛选只是减少一种结构建模疑虑，不能提升为高置信度。

{md_table(best_safe, cols, headers)}

A156 的 V、L、M 替换都出现负向信号（分别 −1.14、−0.92、−0.59），A226 的 K、R、E、Q 替换也出现负向信号（−0.87、−0.72、−0.54、−0.51）。这两个位点值得在初筛中覆盖不同替换类型，但 ΔΔG 排名不能证明机制，也不支持把单突变值相加来预测双突变。

若只做第一轮小规模实验，可优先考虑 **A156V、A156L、A226K、A226R、V134I**，同时保留 WT 对照并测量蛋白表达/折叠状态与结合，避免把蛋白量或结构变化误认为结合变化。

按分数直接排序的前 12 条如下，用于保留全部候选信息。T244K、T234A、A243K 等虽然负值较大，但邻近截断处，建议先核对/改善结构上下文再决定优先级。

{md_table(best_all, cols, headers)}

没有任何单突变达到 ΔΔG ≤ −1.5，不能报告本批次找到“高幅度、高置信度的结合增强突变”。

## 4. 减弱结合及敏感位点

全部 ΔΔG ≥ 1.5 的突变：

{md_table(worst, cols, headers)}

**E154P +2.13、E155P +1.77、E154I +1.71** 没有触发 10 Å 截断标记，可作为减弱方向的探索候选，但它们也不是当前结构中的直接界面残基。E154 的不同替换方向相反（P +2.13，G −0.62），说明应按具体替换讨论，不能笼统认定 E154 的任意突变都会减弱结合。

**R218 是分数上最敏感的位点**：19/19 个替换均为正值，范围 +1.26 至 +1.67，中位数 +1.55，其中 12/19 达到 +1.5。但其距 I/J 伙伴约 33.1 Å、邻近截断端，不能据此认定它是已验证的直接结合热点。保守替换 R218K 也给出 +1.55，可作为后续结构核查或实验对照点。

按该位点 19 种替换的平均值排序：

{md_table(hot, ['Site', 'Mean_DeltaDeltaG_kcal_mol', 'Median_DeltaDeltaG_kcal_mol', 'Min_DeltaDeltaG_kcal_mol', 'Max_DeltaDeltaG_kcal_mol', 'N_ge_plus_1_5', 'Near_truncated_terminus_10A'], ['位点', '平均 ΔΔG', '中位 ΔΔG', '最小值', '最大值', '≥1.5 数量/19', '距截断端≤10 Å'])}

E152 的 19 种替换全部为正，均值约 +0.895、中位数 +0.89，且没有 10 Å 截断标记，表现为较一致的温和减弱信号。D213 的 5 条较大正值全部位于截断端，应格外谨慎。

## 5. 为什么全部是 Low，以及结构限制

模型根据结合前后溶剂可及面积变化判定是否为界面残基。原代码中，非界面突变会被统一覆盖为 Low；这个标签不提供数值误差或概率。本批次 77 个位点全部为非界面。独立重算原始坐标中位点到 I/J 伙伴的最短重原子距离为 **13.2355–40.0084 Å**，与结果注释的最大差异不足 0.000051 Å，支持它们在该构象中并不直接接触这些伙伴。

因此，这些分数可作为间接影响或结构影响的筛选线索，不能直接当作接触界面的优化证据，也不能由距离推断这些位点在生物学上完全无作用。若研究目标是直接增强 L 与 I/J 的结合，下一步应先核对是否选中了所需伙伴及相关功能构象，并评估是否需要扫描实际界面或改善缺失结构。

有 **{summary['near_truncation_mutations']} 条（{summary['near_truncation_mutations']/n:.2%}），即 38/77 位点** 位于截断端 10 Å 范围内。其中 V166、D213、K236、H242 四个位点就在截断端。23 条较大正值中，20 条有邻近截断标记；19 条 ΔΔG ≤ −0.5 的负值中，10 条有该标记。中性封端减少端基电荷处理问题，但无法恢复缺失片段的原子、接触和构象约束。

这些替换是同一结构模型下的相关预测，不是生物学重复。本报告只做描述性汇总，不计算 p 值、置信区间或“统计显著”。输出分数只保留两位小数，0.00 是舍入后为零；较小的名次差异不支持精细的可靠性排序。不要将 ΔΔG 直接当作酶活、转录活性、表达量、耐热性或适应度，也不据此报告实际亲和力倍数。

## 6. 来源及复现

合并表保留了 4 个扫描来源和各自的协议指纹：

{md_table(pd.DataFrame([{'source': k, 'n': v} for k, v in summary['validation']['source_counts'].items()]), ['source', 'n'], ['来源', '条数'])}

不能把不同协议指纹解释为完全同一次几何准备。已有过渡审计记录：最早 38 条与后续批次使用了独立最小化的 WT 缓存，存在轻微坐标差别；同一几何上的旧/新后端复核日志结果一致。本分析保留来源，没有重算或校正这些差别，也没有凭此估算模型误差。

生成文件：`ranked_enhancing.tsv`（全部 83 条负值）、`ranked_weakening.tsv`（全部 1377 条正值）、`site_summary.tsv`（77 位点）、`screening_candidates.tsv`（ΔΔG≤−0.5 或≥1.5，共 42 条）、`summary.json`（统计与完整性核验）。排名表保留原结果字段和来源，并附加位点与截断上下文标记。

![总体分布与位点分数]({OUT / 'overview.png'})

图中竖线表示同一位点 19 种替换的最小至最大分数，不是误差条或置信区间。缺失编号段以斜线背景标出，不跨缺失段连线。

![全部突变热图]({OUT / 'mutation_heatmap.png'})

热图采用以零为中心的对称线性色标，保留全部数值，不裁剪、不平滑。灰色叉号表示未预测的野生型身份；三个面板之间保留原始编号跳跃。顶部标记显示截断上下文。底层数据见原合并结果及排名表。

本目录 `analyze_results.py` 可在 Python 3.12 环境复现，依赖版本见 `requirements_analysis.txt`。所有图形和表格由已完成结果生成，未重新运行预测，未修改原始结果。文件哈希、软件版本和绘图参数见 `analysis_manifest.json`。
"""
    (OUT / "report_zh.md").write_text(text)


def main():
    df, original_columns, validation = load_and_validate()
    sites = summarize_sites(df)
    bands = band_counts(df[SCORE])
    assert sum(n for _, n in bands) == len(df)
    extra = ["WT", "Position", "Mutant", "Structure_caution", "Source_scan_name"]
    columns = original_columns + extra
    enhancing = df[df[SCORE].lt(0)].sort_values([SCORE, "Position", "Mutation"], ascending=[True, True, True])
    weakening = df[df[SCORE].gt(0)].sort_values([SCORE, "Position", "Mutation"], ascending=[False, True, True])
    write_tsv("ranked_enhancing.tsv", enhancing[columns])
    write_tsv("ranked_weakening.tsv", weakening[columns])
    write_tsv("site_summary.tsv", sites)
    candidates = df[df[SCORE].le(-.5) | df[SCORE].ge(1.5)].copy()
    candidates["Screening_group"] = np.where(candidates[SCORE].lt(0), "enhancing_ddG_le_minus_0_5", "weakening_ddG_ge_plus_1_5")
    candidates = candidates.sort_values(["Screening_group", "Near_truncated_terminus_10A", SCORE, "Mutation"])
    write_tsv("screening_candidates.tsv", candidates[columns + ["Screening_group"]])
    summary = {
        "n_mutations": len(df), "n_sites": len(sites), "substitutions_per_site": 19,
        "position_ranges": RANGES, "score_units": "kcal/mol",
        "score_sign": "positive = weaker; negative = stronger binding",
        "min": float(df[SCORE].min()), "max": float(df[SCORE].max()),
        "mean": float(df[SCORE].mean()), "median": float(df[SCORE].median()),
        "positive": int(df[SCORE].gt(0).sum()), "negative": int(df[SCORE].lt(0).sum()),
        "zero_rounded": int(df[SCORE].eq(0).sum()),
        "score_bands": bands, "model_interface": df.Model_Interface.value_counts().to_dict(),
        "model_confidence": df.Model_Confidence.value_counts().to_dict(),
        "near_truncation_mutations": int(df.Near_truncated_terminus_10A.sum()),
        "near_truncation_sites": int(sites.Near_truncated_terminus_10A.sum()),
        "at_truncated_terminus_mutations": int(df.Truncated_terminus.sum()),
        "at_truncated_terminus_sites": sites.loc[sites.Truncated_terminus.eq(1), "Site"].tolist(),
        "partner_distance_min_A": float(sites.Partner_min_heavy_distance_A.min()),
        "partner_distance_max_A": float(sites.Partner_min_heavy_distance_A.max()),
        "protocol_counts": df.Protocol_sha256.value_counts().to_dict(),
        "validation": validation,
    }
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "path", "savefig.facecolor": "white"}):
        histogram_edges = make_overview(df, sites, bands)
        heatmap_settings = make_heatmap(df, sites)
    make_report(df, sites, summary)
    validation["source_sha256_after"] = sha(RESULTS)
    assert validation["source_sha256_after"] == validation["source_sha256_before"]
    validation["original_results_unchanged"] = True
    write_json("summary.json", summary)
    dependencies = {name: version(name) for name in ["numpy", "pandas", "matplotlib", "pillow", "packaging",
                    "python-dateutil", "six", "contourpy", "cycler", "fonttools", "kiwisolver", "pyparsing"]}
    (OUT / "requirements_analysis.txt").write_text("".join(f"{name}=={ver}\n" for name, ver in dependencies.items()))
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "dependencies": dependencies,
        "sources": {str(p): sha(p) for p in [RESULTS, INPUT, PDB,
            ROOT / "mutabind2_tool/engines/single.py",
            COMBINED / "status.json", COMBINED / "baseline_transition_validation.json"]},
        "script_sha256": sha(Path(__file__)),
        "statistics": "descriptive only; rounded saved scores; no hypothesis tests or confidence intervals",
        "unit": "single substitution; 19 correlated predictions per site; no experimental replicates",
        "selection_rules": {"enhancing_ranking": "ddG < 0", "weakening_ranking": "ddG > 0",
            "screening_candidates": "ddG <= -0.5 or ddG >= 1.5; flags retained",
            "preferred_enhancing_subset": "ddG <= -0.5 and Near_truncated_terminus_10A == 0"},
        "figures": {"overview": {"histogram_bin_edges": histogram_edges,
                    "site_center": "median", "site_spans": "min-max over 19 substitutions; not uncertainty"},
                    "mutation_heatmap": heatmap_settings},
        "transforms": ["parse WT/position/mutant", "rank saved scores", "summarize 19 substitutions per site"],
        "exclusions": "none in plots or site summary; sign and threshold subsets explicitly named",
        "files": {p.name: sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name != "analysis_manifest.json"},
    }
    write_json("analysis_manifest.json", manifest)
    print(json.dumps({k: summary[k] for k in ["n_mutations", "n_sites", "positive", "negative", "zero_rounded", "score_bands", "model_confidence"]}, ensure_ascii=False, indent=2))
    print("Validation passed; original predictions unchanged.")
    print(f"Saved analysis to {OUT}")


if __name__ == "__main__":
    main()
