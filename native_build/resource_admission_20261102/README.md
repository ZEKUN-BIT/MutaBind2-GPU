# 256–374 扫描的资源等待与恢复入口

该独立运行扩展用于 `projects/7mkj_256_374_20261102`。科学包、计算输入、部署二进制和协议保持原设置。

原入口在 VDW/PB 每个原生进程启动前检查可用系统内存；低于 30 GiB 时直接失败。本扩展会先等待资源恢复，并处理等待结束到实际启动之间再次触发的同一准入拒绝。等待不占用新的能量计算进程。

GPU 侧在启动前要求至少 5 GiB 空闲显存；仅当原始 PB 日志同时明确包含显存准入不足及 `ERROR STOP 84` 时，保存失败日志并最多重试三次。其它原生失败和数值错误继续正常报错。

保留原有 30 GiB 系统内存门槛、2 GiB 原生显存余量和能量并发数 1。已有批次通过原校验后复用；准备完成的失败批次复用冻结的 FoldX/最小化坐标。准备尚未完成的批次保留原尝试目录，在新尝试目录重做准备。

## 恢复

今后恢复这个扫描请使用此入口：

```bash
nohup bash /home/zek/mutabind2-tool/native_build/resource_admission_20261102/resume.sh \
  >> /home/zek/mutabind2-tool/projects/7mkj_256_374_20261102_inputs/resource_admission_launcher.log \
  2>&1 < /dev/null &
```

原 `projects/7mkj_256_374_20261102_inputs/run_prediction.sh` 保留原文件，未加载本扩展。

## 查看状态

```bash
/home/zek/mutabind2-tool/mutabind2 scan status \
  --output /home/zek/mutabind2-tool/projects/7mkj_256_374_20261102
```

等待时阶段显示 `waiting_native_resources`，当前工作目录中的 `native_resource_admission.json` 记录内存、显存门槛及实时值。

## 验证

2026-11-03 部署前，15 项模拟测试通过，包含本次 VDW 完成后内存不足导致 PB 启动失败的回归案例。没有为测试运行额外的 GPU/物理计算。

```bash
cd /home/zek/mutabind2-tool
.env/bin/python -m unittest discover \
  -s native_build/resource_admission_20261102 -p 'test_*.py' -v
```

恢复前快照及运行证据保存在 `incident_20261103/`。该快照记录已完成的 1805 条预测和第 96 批冻结结构的校验值。
