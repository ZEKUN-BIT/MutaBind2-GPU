# 流水线内存准入与恢复修复

2026-10-30 14:45，第 4 批预取在局部表面积计算入口遇到可用内存低于
30 GiB。原函数在启动原生进程之前立即报错，调度器又立即终止第 3 批，
使其第 8 项正在执行的能量计算中断。此前已保存 171 条预测，第 3 批
另有 7 项完整能量缓存。

## 修复

- `memory_admission.py` 在原表面积入口之前等待 30 GiB 准入条件。
  如果内存在两次检查之间下降，只重试原函数尚未启动进程的准入错误。
  原生错误、运行超时和内存余量错误继续上报。等待时状态显示
  `waiting_surface_memory`，并写入 `surface_admission.json`。
- `continue_scan.py` 在预取批次失败时等待当前批次完成并归档，再报告
  具体失败原因。预取和单路能量计算仍按原流水线顺序执行。
- `recover.py` 校验已准备批次的输入、科学协议、WT、FoldX、全部最小化
  与几何检查点，保留冻结坐标。从表面积阶段继续，复用通过原校验的
  能量缓存。恢复检查点记录冻结输入的 SHA256，后续恢复拒绝篡改。
- 原科学包、CHARMM 后端、模型、阈值和已部署脚本均未修改。科学协议仍为
  `91c31d483e948c0a961e9f58690f29591132467dc7cbd7f35271e2cb3201cece`。
  本扩展单独登记在扫描目录 `memory_admission_deployment.json`。

## 验证

20 项测试通过，含原流水线回归、低内存等待、准入检查竞争、原生错误
传播、中断、预取失败后当前批次提交、停止标记和冻结坐标篡改拒绝。
测试使用隔离的临时 scratch 和模拟工作进程，不启动额外 PB/GPU 任务。
第一次测试运行遇到旧测试 scratch 中遗留的空 owner.json；隔离 scratch
后完整重跑通过。没有修改这些无关 scratch 的归属标记。

`before_recovery.json` 记录恢复前 171 条预测、两批冻结坐标、80 个 WT
缓存文件、科学包和 7 项能量结果的哈希，以及根因错误。

## 使用

```bash
/home/zek/mutabind2-tool/.env/bin/python /home/zek/mutabind2-tool/native_build/pipeline_memory_20261030/continue_scan.py run
```

监控：

```bash
watch -n 5 '/home/zek/mutabind2-tool/.env/bin/python /home/zek/mutabind2-tool/native_build/pipeline_memory_20261030/continue_scan.py status'
```

原监控命令仍然读取同一进度。后续恢复请使用本目录的 `run` 入口。
暂停继续使用当前 GPU 扫描目录的 `STOP_AFTER_BATCH`；历史三个扫描的
停止标记必须保留。
