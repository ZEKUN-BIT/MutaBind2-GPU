# 7MKJ L 链 256–374 位点突变预测分析

分析对象为 L 链与 I+J 两条链的结合，全部为单点突变。119 个连续位点、每位点 19 种非野生型替换，共 2261 条结果，均完整且与原始请求逐一对应；逐批分数、标签、结果及特征表哈希、野生型序列核对通过。任务完成时间：2026-11-04 01:04:34（北京时间）。

## 总体结果

ΔΔG = 突变体结合自由能 − 野生型结合自由能，单位 kcal/mol。负值预测结合增强，正值预测结合减弱。这里的分数是将物理计算、FoldX、PROVEAN 及结构特征输入 MutaBind2 模型后的预测，不等于单独的 PB 能量，也不是所有能量项的直接加和。

- 预测增强：62 条（2.74%）；预测减弱：2195 条（97.08%）；保存精度下为零：4 条。
- 均值 0.442，中位数 0.41，范围 -0.89 至 +2.23。
- ΔΔG ≤ −0.50 有 11 条；≤ −1.50 为 0；≥ +1.50 有 12 条。
- ±0.50 只是本报告整理候选所用的描述性阈值；±1.50 是项目代码用于幅度分类的界限。二者都不是本批结果的误差界、统计显著性或实验成功率。大量小正值不能据此全部称为确定有害。

## 增强方向：优先探索 V275，其次跨位点候选

| 突变 | ΔΔG (kcal/mol) | 模型标签 | 野生型至 I/J 最近重原子距离 (Å) |
| --- | --- | --- | --- |
| V275T | -0.89 | Low | 14.99 |
| A318K | -0.81 | Low | 30.99 |
| S334Q | -0.81 | Low | 28.42 |
| V275K | -0.76 | Low | 14.99 |
| T304K | -0.74 | Low | 20.63 |
| V275Q | -0.73 | Low | 14.99 |
| V275R | -0.70 | Low | 14.99 |
| V275S | -0.65 | Low | 14.99 |
| S312N | -0.55 | Low | 27.87 |
| S334E | -0.53 | Low | 28.42 |
| S334K | -0.51 | Low | 28.42 |

V275 的 19 种替换中 11 种为负值，位点中位数 −0.40；其中 T、K、Q、R、S 五种替换 ≤ −0.50。相比仅有一个极端分数，它提供了多个同方向的候选，但这些替换互相关联，不能当作独立重复实验。按分数优先考虑 V275T、V275K；若希望覆盖不同位点，可增加 A318K、S334Q、T304K。所有这 11 条均为 Low/非界面且没有 10 Å 截断邻近标记；没有找到同时满足 High 且较大幅度增强的候选。

V275T 的原始模型输入 ddVdw 约 −0.011、ddPb 约 −0.022、FoldX_fold 约 −0.917；其最终 −0.89 分数不能解释为 PB 作用大幅增强。A318K/S334Q 的 ddPb 为负，也不能仅凭这一项推断机制，因为模型同时使用其他特征。PROVEAN 本身已参与评分，不能把它作为独立验证证据。

## 减弱方向：区分界面候选和结构影响

| 突变 | ΔΔG (kcal/mol) | 模型标签 | 距截断端点 ≤10 Å |
| --- | --- | --- | --- |
| N323P | 2.23 | Low | 0 |
| K296P | 2.04 | Low | 0 |
| K296I | 1.94 | Low | 0 |
| R373P | 1.80 | High | 0 |
| R260G | 1.75 | Low | 1 |
| N323I | 1.70 | Low | 0 |
| R279P | 1.70 | Low | 0 |
| R260A | 1.69 | Low | 1 |
| R260T | 1.68 | Low | 1 |
| R260D | 1.66 | Low | 1 |
| R260S | 1.63 | Low | 1 |
| N323T | 1.52 | Low | 0 |

R373P（+1.80）是本批唯一同时为 High 且达到 +1.50 的突变，值得优先作为结合减弱候选。R373 所有 19 种替换的分数都为正，其中其余 18 种为 +0.21～+0.64。R373P 的 FoldX_fold 约 +2.96，提示同时存在折叠稳定性改变的可能；验证结合时应同时检查蛋白稳定性/可溶性，避免把结构受损全部归因于直接界面接触变化。

N323P（+2.23）、K296P（+2.04）、K296I（+1.94）幅度更大，但均为 Low。N323P/K296P 的 FoldX_fold 分别约 +4.40/+4.60；这些结果更适合作为可能影响结构与结合的探索候选，不能直接认定为界面热点。R260 的五个 ≥+1.50 结果均靠近截断端点，优先级应更谨慎。

如设计一个小规模验证组合，可选 WT、V275T、V275K、A318K、S334Q、T304K、R373P，并加入 R373K（+0.21）作为同位点较小效应的比较对象。它仍是计算筛选方案，具体选择取决于实验目的；单突变分数不能相加来预测组合突变。

## High / Low 和结构边界如何解释

全部 19 个 High 对应 19 个界面标记：R373 的 18 种替换加 A372W；R373G 为 Low。原始逐批表与存档中的表面积输出一致。界面判定采用每个突变结构的 DRASAP > 0：R373G 为 0，A372W 为 0.006917212，R373P 为 0.00381787；同一位点因此可出现不同标签。该规则对是否严格大于零进行分类，不能把 High 理解为校准后的概率、低误差保证或效应一定大。

171 条结果（9 个位点）距某个结构截断端点 ≤10 Å：F256、K257、Q258、F259、R260、L261、G352、T354、Q357；没有突变恰好位于截断端点。该标记考虑双方结构中的截断端点，不只考虑 L 链。没有标记也不能保证完全不受缺失结构影响。自动片段化和中性封端处理仍不能还原缺失片段；全序列 PROVEAN 也不会补回三维结构。

结果针对当前结构、I+J/L 伙伴划分和现有模型，不能直接推断整套复合物的功能、转录活性、表达量或生物学适应度。当前没有本体系的实验标定或独立重复计算，不能给出本批特异的置信区间，也不应按两位小数的微小差别过度排序。

## 文件与复现

- 原始完整结果：[completed_predictions.tsv](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/completed_predictions.tsv)。
- 全部增强/减弱排序：[ranked_enhancing.tsv](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/ranked_enhancing.tsv)、[ranked_weakening.tsv](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/ranked_weakening.tsv)。
- 23 条幅度筛选候选：[screening_candidates.tsv](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/screening_candidates.tsv)；每位点的 19 种替换汇总：[site_summary.tsv](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/site_summary.tsv)。
- 汇总及原始模型输入：[all_results_with_features.tsv](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/all_results_with_features.tsv)。Feature_ 前缀为模型输入，Contribution_ 前缀为模型归因列，归因值不是独立物理能量。
- 图：[overview.png](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/overview.png)、[mutation_heatmap.png](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/mutation_heatmap.png)，另有 SVG；[图注与文字说明](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/figure_captions_zh.md)。
- [机器可读汇总](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/summary.json)、[审计及来源记录](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/analysis_manifest.json)、[界面标签原始输出核对](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/interface_label_audit.tsv)。
- [分析脚本](/home/zek/mutabind2-tool/projects/7mkj_256_374_20261102/analysis_20261106/analyze_results.py)只读取原始数据、向本分析目录写衍生文件；本次未重新运行任何物理计算。分析依赖版本记录在 requirements-analysis.txt，解释器记录在 analysis_manifest.json。

本次分析仅针对 256–374 批次，未与此前不同来源/基线的 128–255 批次合并比较。
