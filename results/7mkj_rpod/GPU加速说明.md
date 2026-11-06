2026-10-20 补充：GPU PB 入口保持不变；新的表面积、VDW 和 FoldX 优化见 [非PB加速说明.md](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/非PB加速说明.md)。本文的 PB 基准保留历史实测，完整流程的最新时间预算以 [耗时与内存评估.md](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/耗时与内存评估.md) 为准。

**第二轮 GPU 优化已接入：三次 PB 从约 87 秒降到 66 秒，耗时减少约 24%。** 2026-10-20 在本机 RTX 5070 Ti 上，真实 WT 与 E90A 均完成对照，六份完整电势文件与固定种子 CPU 结果逐字节相同。这是 PB 阶段的实测，完整预测还包括 FoldX、建拓扑、最小化、VDW、表面积和模型计算。

本轮把分子表面构建中的探针遮挡判断和网格恢复移到 GPU。整个 REEN 分子表面构建阶段合计由约 44 秒降至 22.5 秒。CPU 继续生成原来的随机点、构造并排序邻居、执行最后的网格清理；已验证的第一轮 GPU 网格迭代和边界电势求和继续使用。

**数值与性能验证**

| 真实结构 | 优化 CPU：三次 PB | 第一轮 GPU | 本轮 GPU | 本轮对照 |
|---|---:|---:|---:|---|
| WT，约 5 万原子 | 542.11 秒 | 87.02 秒 | 66.01 秒 | 三个完整电势文件逐字节相同 |
| E90A | 544.11 秒 | 87.02 秒 | 66.01 秒 | 三个完整电势文件逐字节相同 |

表中均包含保存三份电势文件的开销。本轮 PB 相比优化 CPU 约快 8.2 倍。正式流程不保存电势文件，具体 VDW+PB 耗时见下方的完整流程记录；不能把 66 秒当成一次完整突变预测的耗时。

- 116 个独立 GPU 用例：原 Fortran 与新实现的四个网格，在最终清理前后均逐字节一致。
- 实际 Fortran 接口在 CUDA 和 CPU 模式各通过 118 个用例；探针点、随机数最终状态、四个网格全部一致，包含 513 个邻居触发的容量回退。
- 41 个无效输入以及无 GPU 场景检查通过，错误时不修改调用者的四个网格；实际接口遇到 CUDA 错误明确中止。
- 已部署入口在项目自带 WT/突变小体系中通过能量、收敛轨迹及缓存复用检查。小体系 VDW+PB 从约 7 秒降到 4.5 秒。

这些检查证明已测实例和分支的实现等效性，不代表已完成全部 8949 项突变预测。

**完整预测的时间预算**

第二轮 GPU 部署时，WT 的能量阶段为 VDW 18.01 秒 + PB 62.52 秒 = 80.53 秒，E90A 约 81.04 秒。这些是使用已有几何的历史能量重算记录。新非 PB 优化进一步减少表面积、VDW 和 FoldX 的耗时；当前预算见上方的耗时评估文档。PB 单阶段快 8.2 倍不能解释为完整流程快 8.2 倍。

**保留的计算条件**

继续使用 449×489×505 网格、0.4 Å 间距、全复合物范围、原介电参数、原三次 PB 求解、原收敛标准和随机森林模型。四个 CLCG 种子仍为 `314159265 22222222 33333333 44444444`，没有减少随机探针数量。GPU 几何使用原双精度表达式和运算顺序，不启用 fast-math 或融合乘加。

多个探针可能命中同一个网格点。GPU 使用原子操作，仅将负值恢复为原值的相反数；同一点重复命中不会反复翻转。邻居筛选只提前排除轴向距离明显大于两倍半径和的远处原子，保留原筛选公式与排序。每原子最多传输 512 个邻居；超出容量或不适用时，整段使用原 CPU 流程，且不会再次生成随机点。CUDA 执行失败会中止。正式适配器继续要求三次 SOR 都实际使用 GPU；本轮性能验证还额外要求三次 REEN GPU 成功。

不能跨突变直接复用 B 伙伴的静电能：最小化也会移动 B，E90A 中 B 的能量与 WT 相差约 389.45 kcal/mol。只有输入、协议和运行库均匹配的完整缓存才可复用。测试过的 uint32 网格索引改动几乎没有真实提速，本轮没有采用。

**内存与运行方式**

本轮 PB 主机峰值约 **4.44 GiB**，比第一轮约 3.10 GiB 略高，增加的内存用于先完整接收四个网格、成功后再提交结果。REEN 的 GPU 分配约 1.75 GiB；SOR 的主数组约 2.48 GiB，二者按阶段释放，不同时保留。驱动及显示占用另计。

新 VDW 入口减少目标原子对列表的复制，峰值由约 **14.1 GiB 降到 6.0 GiB**；仍保持全局单路能量锁。每次启动要求至少 30 GiB 可用内存；原生进程地址空间上限 20 GiB，运行中至少保留 10 GiB 系统可用内存，超时 1800 秒。GPU 分配前保留至少 2 GiB 显存。

**当前入口与缓存**

本轮入口：[charmm-pb-gpu-reen](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/tools/charmm-pb-gpu-reen)。第一轮 `charmm-pb-gpu` 和原 `charmm` 入口仍保留，普通项目调用不会自动切换。任务适配器会将入口、实际可执行文件、加载库、PSF、坐标、输入和种子全部纳入缓存身份；不会把旧后端缓存冒充成本轮运行结果。

对本任务，通过 `task_engine.py` 或 `run_batches.py` 显式传入新入口。中性试算恢复命令如下；只用于已有工作区，匹配的完成缓存会复用：

```bash
'/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/.env/bin/python' -u \
  '/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/scripts/task_engine.py' \
  --input '/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/inputs/fragmented/baseline.input' \
  --workspace '/home/zek/mutabind2_runs/7mkj_rpod/baseline_01' \
  --cap-mode neutral --workers 1 --resume-existing \
  --pb-executable '/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/tools/charmm-pb-gpu-reen' --pb-threads 1 \
  --vdw-executable '/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/tools/charmm-vdw-compact'
```

带电对照将输入改为 `sensitivity.input`、工作区改为 `sensitivity_01`、封端模式改为 `--cap-mode charged`。不要同时对同一工作区启动多个任务。

完成试算后，调度器的使用参数是：

```bash
'/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/.env/bin/python' -u \
  '/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/scripts/run_batches.py' \
  --pb-executable '/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/tools/charmm-pb-gpu-reen' --pb-threads 1 \
  --vdw-executable '/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/tools/charmm-vdw-compact'
```

本轮保留 `STOP_AFTER_BATCH`，全量扫描尚未启动。上述调度器读取完整 **8949 项**输入；40–49 个位点只是讨论过的规模，没有自动替换输入。启动完整批次的操作见 [方法与使用说明.md](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/方法与使用说明.md)。完整序列 PROVEAN 已完成 8949 项，无需重新搜索 NR。

**证据与源码**

- [WT 实测](/home/zek/mutabind2_runs/7mkj_rpod/optimization_01/large_wt_reen_save/result.json)、[WT 电势对照](/home/zek/mutabind2_runs/7mkj_rpod/optimization_01/large_wt_reen_save/phi_validation.json)。
- [E90A 实测](/home/zek/mutabind2_runs/7mkj_rpod/optimization_01/large_mut_e90a_reen_save/result.json)、[E90A 电势对照](/home/zek/mutabind2_runs/7mkj_rpod/optimization_01/large_mut_e90a_reen_save/phi_validation.json)。
- [独立数值审查](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v3/reen_independent_review.json)、[部署适配器检查](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v3/reen_adapter_deployed.json)。
- [部署与重编译说明](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/external/charmm-pb-cuda-reen/README.md)、[部署清单](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/external/charmm-pb-cuda-reen/DEPLOYMENT.json)、[新增 CUDA 源码](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/mutabind2/external/charmm-pb-cuda-reen/source/pb_reen_cuda.cu)。

本轮用于完整电势对照的临时网格文件已在核对哈希后删除，保留输入、日志、测量及全部哈希；加上定位阶段，清理约 3.7 GiB。记录：[REEN 清理清单](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v3/reen_cleanup.json)、[定位阶段清理清单](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v3/profile_cleanup.json)。基准输入中的历史 PHI 路径不表示该文件仍存在。

**第二轮 GPU 部署时的试算与调度记录**

两组试算的 7 个能量任务全部完成，新后端的原生日志审计、收敛轨迹和有效缓存均通过；两份最终结果表与上一版逐字节一致。6 次无明显外部争用的能量阶段为 **79.53–81.04 秒**；D213P 与另两个绘图任务同时运行时为 **92.04 秒**，该观测单独保留，不用来宣称统一速度。因可用内存低于 30 GiB 而中止的一次 S506A 尝试已在内存恢复后续算完成，资源限制未降低。

中性结果为 E90A +1.17、D213P +0.58、S506A +0.62 kcal/mol；带电对照 E90A +0.75、S506A +0.48，仍通过既定端基敏感性检查。这些是基于缺失片段结构的近似筛选结果。调度器已识别匹配的新后端并在停止标记处退出，已完成数仍为 **3/8949**，没有启动新批次。

记录：[完整试算核验](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v3/reen_pilot_validation.json)、[调度器检查](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v3/reen_scheduler_validation.json)、[资源中断记录](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v3/resource_interruption.json)、[中断及重复文件清理](/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod/validation/optimization_v3/reen_interruption_cleanup.json)。
