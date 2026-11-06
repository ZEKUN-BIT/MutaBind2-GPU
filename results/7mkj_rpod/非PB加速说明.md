# MutaBind2：PB 之外的加速

局部解析 SASA、两个独立 FoldX 进程并行，以及 VDW 紧凑配对列表均已部署并通过验证。两组已有 pilot 的 7 项能量复核通过，两份最终预测表逐字节不变。本轮没有启动 8,949 项全量扫描；调度器已验证新配置后按 [STOP_AFTER_BATCH](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/STOP_AFTER_BATCH) 退出，仍为完成 3/8,949 项、无活动计算进程。

| 阶段 | 原方式实测 | 优化实测 | 当前结论 |
|---|---:|---:|---|
| CHARMM 解析 SASA，单结构 native 调用 | 9.4–10.8 秒 | 约 0.60 秒；含 Python 适配整体约 0.7 秒/突变 | 已验证并接入 |
| FoldX，两项独立突变 | 串行合计 36.845 秒 | 两个独立进程并行 19.242 秒 | 已验证，默认并发 2 |
| VDW，4 个大结构，无额外审计输出的调用 | 17.43–17.63 秒；约 14.09 GiB | 13.62–13.82 秒；约 6.01 GiB | 4 个结构合计耗时减少 22%，峰值内存减少约 57% |

这些数字只对应表中的阶段和测试。SASA 的速度提升不能作为全流程的提升倍数；也不能把各阶段的倍数相乘。每个突变仍需要独立建模、拓扑处理、100 步最小化、VDW 和三次 PB。

SASA 使用原 CHARMM 的解析算法、原 LJ 半径、1.4 Å 探针和原始 `*_charmm.crd`。它只保留目标残基及所有可能遮挡其表面的原子，随后仍用原公式计算目标残基面积。实际参数给出的保守邻域为 8.35 Å，包含额外 1 Å 数值裕量；不是任意裁剪距离，也没有改用 FreeSASA。原来不被消费的两次面积计算和每批一次 WT surface 调用被移除。RF 中的 WT SASA 特征继续使用原 DSSP 结果。

六个真实结构覆盖中性/带电封端、末端残基以及零/非零界面。MSr、ASAC、ASAMP 的目标原子 double 数组逐字节相同，DRASAP 和界面判断一致。五个已有突变的特征表、RF 预测也逐字节一致。证据：[SASA 原流程对照](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v4_surface/native_oracle.json)、[完整特征与预测复核](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v4/surface_integration.json)。

FoldX 保持“一项突变一个独立进程”，只让两个进程同时运行。E90A、E90D 的突变结构、配套 WT 结构及所有输出能量分量与串行独立运行一致；每个进程实测峰值约 2.2543 GiB。相邻两个任务同时占用内存时，实测峰值之和约 4.51 GiB，仍需给操作系统和其他阶段留余量。证据：[并行验证](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v4_foldx/parallel_report.json)。

没有采用“把多个突变列表交给同一个 FoldX BuildModel 进程”的方式。实测中，后续突变的坐标和能量随批内顺序改变，无法保持原独立运行结果，因此该方案被拒绝。证据：[FoldX 批处理对照](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v4_foldx/report.json)。

VDW 优化保留原 CHARMM 配对次序、排除项、1–4 相互作用和实际参与计算的原子对，减少临时筛选列表的多余容量和深复制；完整主配对列表仍会建立。最终库共通过 37 次 native 验证。四个大结构的全部 192 项能量分量和 10 类哈希对照通过，包括配对、排除、1–4 和力的字节一致性。它没有缩短原来的 999 Å 配对距离，也没有冻结伙伴坐标，也不是 GPU VDW。四个结构正式性能测试合计为 70.12 → 54.69 秒，单结构为表中范围。证据：[VDW 验证汇总](/home/zek/mutabind2_runs/7mkj_rpod/optimization_vdw_v4/SUMMARY.md)、[正式性能报告](/home/zek/mutabind2_runs/7mkj_rpod/optimization_vdw_v4/validation/large_performance_report.json)、[完整数值对照](/home/zek/mutabind2_runs/7mkj_rpod/optimization_vdw_v4/validation/large_audit_report.json)。部署库已冻结，SHA256 为 `7df3b9d75cbfc9112c295eaaf84030053810a8c00d5d53f72f1594218f2bf89e`。

正式 pilot 接入复核中，7 项 VDW 调用均约 14.009 秒，包含适配器 0.5 秒轮询粒度。已有 PB 日志经过身份与哈希核验后复用，逐字节未变；这轮迁移没有重新计算 PB，报告中保留的 PB 耗时属于此前验证运行。能量审计值完全相同，两份预测表逐字节一致。证据：[最终 pilot 验证](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v4/pilot_validation.json)。

已部署的运行入口与参数如下。

| 用途 | 入口或参数 |
|---|---|
| 原已验证的 GPU PB | [charmm-pb-gpu-reen](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/tools/charmm-pb-gpu-reen)，通过 `--pb-executable` 指定；`--pb-threads 1` |
| 紧凑 VDW 列表 | [charmm-vdw-compact](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/tools/charmm-vdw-compact)；通过 `--vdw-executable` 指定，并同时提供 `--pb-executable` |
| 局部解析 SASA | `--surface-backend local`，默认开启；需要原方式时显式选择 `--surface-backend original` |
| 批次调度中的 FoldX 并发 | [run_batches.py](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/scripts/run_batches.py) 使用 `--foldx-workers 2`，默认 2，可改为 1 |
| 单个工作区的 FoldX 并发 | [task_engine.py](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/scripts/task_engine.py) 使用 `--workers 2`，默认 2，可改为 1；该选项不会并行 VDW/PB |

更换 VDW 可执行文件时，已有 PB 结果只有在几何、参数资源、随机种子、PB deck、运行库、线程数和日志全部核验一致后才复用。已验证日志在复制前后仍核对哈希和审计值，不能把中途被改动的日志当成有效结果。失败或未知版本的记录不会捐献缓存。SASA 遇到未知科学模板、输入变化或计算错误会停止，不静默改变算法。批次归档保存实际优化 deck、运行日志、输入哈希及原始 SASA 坐标，再清理工作区。详见[适配器独立复审](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v4/adapter_review.json)。

资源保护没有放松：VDW/PB 全局同时只运行 1 个，启动要求至少 30 GiB 可用内存，地址空间上限 20 GiB，保留至少 10 GiB 可用内存，并有超时和进程组终止处理。FoldX 上限为两个独立进程，每个使用更严格的 8 GiB 地址空间上限，另有内存和超时检查。降低 VDW 内存峰值不会自动提高能量阶段并发数。

按现有测试，SASA 每突变省约 9–10 秒；FoldX 两项并行比这两项串行合计省约 17.6 秒，摊到这组测试约每突变 8.8 秒。后者是批次吞吐收益，单项 FoldX 本身并未缩短到一半。四个大结构的 VDW 阶段每项另省约 3.8–4.0 秒。

在已有 GPU PB 优化的基础上，新一轮优化预计使整流程再减少约 15%–20% 的时间。按两个独立 FoldX 进程、单个 VDW/PB 能量作业、既有 WT 与 PROVEAN 缓存估算，采用 **95–112 秒/突变** 的计划预算：

| 待测规模 | 计划耗时 |
|---|---:|
| 40 个位点 × 19 种替换＝760 项 | 约 20–24 小时 |
| 49 个位点 × 19 种替换＝931 项 | 约 25–29 小时 |
| 全部 8,949 项 | 约 10–12 天 |

这些是现有阶段实测的外推，未完整运行上述规模；暂停、重试和其他程序争用会增加时间。最终吞吐仍需用实际选定位点批次测量，不能据局部 SASA 数字宣布全流程“16 倍加速”。
