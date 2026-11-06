# 单批预备流水线

调度器在当前批次执行 PB/VDW 能量计算时，准备下一批突变结构。
下一批完成 FoldX、拓扑、最小化、表面积和接触计算后，等待独立的能量放行文件。
当前批次结束、结果归档并校验后，调度器才放行下一批。

部署目录：`/home/zek/mutabind2-tool/native_build/scan_pipeline_20261030`。
扫描目录：`/home/zek/mutabind2-tool/projects/7mkj_128_255_20261030_optimized`。
合并结果：`/home/zek/mutabind2-tool/projects/7mkj_128_255_20261030_combined/completed_predictions.tsv`。

## 计算协议

- 调用现有 `mutabind2_tool.scan.task`，不修改科学计算包。
- 复用当前 WT 基线，不重新创建扫描计划或重新编号。
- 协议指纹保持 `1f60065e9cfc4796bc1c5ada830c8b68621d82ff60537b6dc958b72463dadb86`。
- 深度为当前批次加下一批；能量计算并发数仍为 1。
- 保留现有全局 `runtime/energy.lock`，以及 PB、VDW、FoldX 内存和进程限制。
- 下一批启动前要求至少 30 GiB 可用内存。
- 继续使用原归档、校验和与清理函数。归档成功后才删除工作目录。
- 流水线源文件哈希单独记录在扫描目录的 `pipeline_deployment.json`；恢复时校验。

## 查看状态

在任何目录执行：

```bash
watch -n 5 cat /home/zek/mutabind2-tool/projects/7mkj_128_255_20261030_optimized/pipeline_status.json
```

`energy` 表示当前批次，`preparation` 表示正在准备或等待的下一批。
下一批处于 `waiting_pipeline_energy` 时，结构准备已完成，尚未开始能量计算。
`lowest_available_gib` 是本次调度器启动后观测到的最低可用内存。
主日志为扫描目录的 `pipeline.log`，调度事件为 `pipeline_events.jsonl`。
原有 `continue_scan.py status` 仍可查看合并进度。

## 暂停与恢复

创建扫描目录下的 `STOP_AFTER_BATCH` 后，当前批次会算完并归档，调度器随后退出。
已经准备的下一批可能保留一个等待放行的进程；恢复时会直接接管它。
用户确认继续后，移除该暂停文件，并使用以下命令；不要另行启动旧的扫描调度器：

```bash
cd /home/zek/mutabind2-tool
.env/bin/python -u native_build/scan_pipeline_20261030/pipeline.py \
  --output /home/zek/mutabind2-tool/projects/7mkj_128_255_20261030_optimized \
  --combined /home/zek/mutabind2-tool/projects/7mkj_128_255_20261030_combined
```

调度器通过 `pipeline.lock`、`scan.lock`、`continuation.lock` 排除重复调度。
`--handoff` 用于旧调度器仍运行时的批次边界交接。
原扫描 `7mkj_128_255_20261030` 的暂停文件应保留，以免重算已有结果。

## 验证与本次重启恢复

```bash
.env/bin/python -m pytest -q tests/test_scan.py native_build/scan_pipeline_20261030/test_pipeline.py
```

39 项通过，包含 29 项原扫描测试和 10 项流水线/恢复测试。
流水线测试用真实子进程验证准备重叠、能量放行顺序、与串行结果一致、孤立进程接管、
准备失败清理、暂停、锁互斥、错误协议拒绝、内存不足延迟、旧调度器交接和归档恢复。

2026-10-30 10:38 左右，主机发生重启，原调度器和刚启动的预备进程退出。
日志没有留下能确定重启原因的记录。
已提交的 38 条预测保留。第 1 批的 19 个已准备结构通过原几何审计，
12 个完整能量结果通过原输入/程序/日志/输出校验后复用。
`recover_energy.py` 仅用于这次已完成准备的批次恢复：验证计划、协议、基线、结构、
最小化和 FoldX 输出，复用固定几何与校验通过的能量，重新计算便宜的表面积/接触输出，
然后继续调用原能量和预测函数。未完成或损坏的能量会由原缓存机制重新计算。
第 2 批重启前只完成了部分 FoldX，保留旧尝试，使用新的工作目录重做准备。

重启后恢复前的文件哈希、原启动记录、空文件和实时观测保存在 `validation/`。
恢复前状态同时保存在 `validation/energy_recovery.json` 和第 1 批归档内的
`jobout/energy_recovery.json`。归档后的批次状态保留缓存复用数量。

## 首次实际批次切换

2026-10-30 10:57:08，第 1 批完成归档，随后放行第 2 批，第 3 批开始准备。
第 2 批的准备于 10:48:34 开始，10:55:50 完成，共 436.04 秒，
整个准备阶段与上一批处理重叠。完成后等待了 77.206 秒，期间没有能量放行文件，
也没有 PB/VDW 计算。归档完成至下一批放行的间隔为 0.034 秒。

此次切换后合并结果共 57/1463 条，包含原有 38 条和新增 19 条。
206 个受保护文件未改变，12 个复用突变能量输出及 WT 能量输出在归档中均验证为原字节。
归档验证通过后才清理第 1 批工作目录。观测期间最低可用内存为 30.327 GiB，
没有观察到两批能量计算重叠。

完整验证结果见 `validation/end_to_end.json`；分阶段实时证据见
`validation/live_observations.jsonl`。第 2 批能量和第 3 批准备继续在后台运行，
调度器 PID 为 12031（仅代表本次启动，重启后应读取状态文件中的新 PID）。

流水线减少批次之间等待准备的时间。具体收益取决于准备与能量计算的重叠程度、
CPU 争用、内存等待及归档时间；应依据实时观测计算，不能把单个阶段的加速比当作总加速比。
