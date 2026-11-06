# 2026-10-30：REEN 探针网格计算优化

新后端将一个可达探针覆盖的网格分给同一个 warp 中的 32 个 GPU 线程处理。原版本由单线程遍历这片网格。邻居判断、随机探针、双精度几何表达式、舍入顺序、网格、收敛条件和 MutaBind2 模型均保持原设置。

## 实测

在同一台 RTX 5070 Ti 上，用完全相同的冻结结构、输入参数和资源限制比较。测量时显卡也被另一个 ProteinMPNN 任务使用。PB 总耗时包含 AB/A/B 三次求解及保存完整 PHI 的验证开销。

| 结构 | 旧版 PB | 新版 PB | 耗时减少 |
|---|---:|---:|---:|
| WT | 75.019 秒 | 62.013 秒 | 17.34% |
| 当前扫描中的 V130A | 75.021 秒 | 62.014 秒 | 17.34% |

探针阶段合计从约 18.94 秒降至 5.48 秒。每个突变预计节省约 13 秒；按切换后 1,425 个待算突变估计，约节省 5.1 小时。该估计来自两个代表结构，实际收益受结构、系统负载和其它计算阶段影响。17.34% 是 PB 耗时下降，不能直接当作整个预测流程的下降比例。

## 数值验证

- 116 个独立 CUDA 案例：四个网格在清理前后均与原 Fortran 计算逐字节一致。
- 118 个实际 Fortran 桥案例：网格、随机点、RNG 最终状态一致，包括 513 邻居的 CPU 回退。
- 41 个错误输入案例和 1 个无 GPU 案例通过。
- WT 的三份完整 PHI 与已保存的 CPU 参考哈希一致。
- 当前扫描中 V130A 的 PSF、CRD、输入文件哈希与归档一致；新旧后端的三份完整 PHI 逐字节一致。
- 两种真实结构的静电能、网格大小、边界和迭代次数一致，均确认三个组件成功执行 CUDA REEN、边界及 SOR。

验证汇总：`validation/summary.json`。各测试原始报告、完整计算日志和 PHI 文件均保留。

## 独立安装及接续

新入口：`/home/zek/mutabind2-tool/bin/charmm-pb-gpu-warp`。

新动态库：`/home/zek/mutabind2-tool/toolchain/external/charmm-pb-cuda-reen-warp/lib/libchmm.so`。

SHA256：`492fbfadde7d41226d25a3656a7aef7d6fc94c0cd7f55b18ad4742bda504c71d`。

旧入口及动态库保持原样。新入口使用独立路径，原扫描在第二批完成后暂停，已完成的结果保留原协议记录。剩余突变在新扫描目录中运行，并重新计算 WT 基线。每个突变保留用户原始 Result_Id。

原扫描已完成 38 项，新扫描接续 1,425 项。两份计划的准备结构相同，保留的 1,425 个 PROVEAN 分数逐项相同；协议差别仅为 PB 后端和删除已完成突变后的 PROVEAN 子集。

需要区分内核等价性与重新生成结构的重复性：两次 WT 最小化的输入坐标完全相同，但最小化后的坐标均方根差为 0.000324 Å，最大差为 0.0383 Å。新 WT 的静电结合能为 222.226 kcal/mol，原 WT 为 222.264 kcal/mol，VDW 分别为 -413.901408 和 -413.891376 kcal/mol。这是两个独立生成的基线，不能声称整个接续流程逐字节重现了原基线。合并表保留每行对应的扫描和协议，两份 WT 基线也完整保留。GPU 等价性验证使用相同冻结坐标。

在新 WT 坐标上额外运行旧后端后，两版的静电能、网格和迭代次数完全一致，确认上述基线差异来自重新最小化。该次不保存 PHI 的检查测得旧后端 72.02 秒、新后端实际基线运行 57.53 秒；完整记录见合并目录中的 `baseline_transition_validation.json`。这些测量与上表保存 PHI 的测量条件不同，不合并计算平均加速率。

- 原扫描：`/home/zek/mutabind2-tool/projects/7mkj_128_255_20261030`
- 优化扫描：`/home/zek/mutabind2-tool/projects/7mkj_128_255_20261030_optimized`
- 合并结果：`/home/zek/mutabind2-tool/projects/7mkj_128_255_20261030_combined`

`continue_scan.py` 在独占接续锁后运行新扫描，每 10 秒刷新合并表。`completed_predictions.tsv` 只含完成项；`predictions.tsv` 含全部 1,463 项及当前状态。新增 `Source_scan` 和 `Protocol_sha256` 两列记录每项结果的来源。

查看总进度：

```bash
watch -n 10 '/home/zek/mutabind2-tool/.env/bin/python /home/zek/mutabind2-tool/native_build/pb_reen_warp_20261030/continue_scan.py status'
```

进程中断后，用同一个 `continue_scan.py run` 接续；保留原扫描的 `STOP_AFTER_BATCH`，避免两份任务重复计算。

## 构建记录

`build_candidate.py` 只重新编译 REEN CUDA 对象，复用已经验证的 Fortran、SOR 和边界对象，沿用原浮点编译选项。`build_manifest.json` 保存全部命令和输入哈希。`benchmark_pb.py` 是原测试工具的隔离副本，仅把非阻塞能量锁改为阻塞锁，并保留原路径、内存限制和超时保护。

不要在新后端被扫描使用期间重编或覆盖动态库。需要继续实验时，创建新的独立构建目录。
