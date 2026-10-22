# MutaBind2-GPU：GPU 加速改进版

这是 [ZEKUN-BIT](https://github.com/ZEKUN-BIT) 基于 [MinghuiLab/MutaBind2](https://github.com/minghuilab/MutaBind2) 开发的**非官方 GPU 加速改进版**。保留原随机森林模型、原方法引用和上游版权声明；不是原作者发布的新版本，也不代表原作者背书。

本版主要加速受支持的 CHARMM PB 网格迭代、边界电势求和及 REEN 探针计算，并加入紧凑 VDW 内存分配、局部 SASA、WT/PROVEAN 缓存、断点续算和结构断链映射。PROVEAN 与最终随机森林仍在 CPU 上执行。

## 实测范围

已验证的通用 GPU 入口适用于**独立单点替换与饱和扫描**。原多点组合接口保留，但尚未验证组合多点 GPU 流程。

在 RTX 5070 Ti、固定 1CSE WT/L45P 几何下，VDW＋PB 阶段从 CPU 的约 46.6–50.2 秒缩短到约 4.5 秒。相同网格、种子、收敛条件和输入下，能量在日志打印精度内一致。**这是能量阶段耗时，不是整次预测的加速比，也不是生物学准确率评估。** [验证说明](docs/VALIDATION.md)

## 安装与运行

完整步骤见 [安装指南](docs/INSTALL.md) 和 [英文 README](README.md)。Python 3.11+；完整预测还需要自行安装 CUDA、VMD、合法来源的 CHARMM/FoldX/NAMD、DSSP、PROVEAN、BLAST+、CD-HIT、本地 NR，以及 R 包 `randomForest` 和 `forestFloor`。仓库不附带这些工具的二进制、运行环境或 NR 数据。

```bash
git clone https://github.com/ZEKUN-BIT/MutaBind2-GPU.git mutabind2-gpu
cd mutabind2-gpu
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[test]'
cp config/tools.example.json config/tools.json
# 修改 config/tools.json，并按安装指南构建 GPU 库
./mutabind2 doctor --scope gpu

./mutabind2 scan run \
  --pdb examples/1cse.pdb --partner1 E --partner2 I \
  --mutations examples/mutations.tsv \
  --output projects/1cse_scan --plan-only
./mutabind2 scan run --resume --output projects/1cse_scan
```

结果 `completed_predictions.tsv` 中，正 ΔΔG 表示结合减弱，负值表示增强，单位 kcal/mol。PROVEAN 只是序列变异影响特征，不是结合 ΔΔG。默认对连续结构片段计算 PROVEAN；完整序列和片段分数可能不同。无同源序列退回自身计分、数据库变化后缓存复用等边界问题已写入 [已知限制](docs/KNOWN_LIMITATIONS.md)。

请引用原始 [MutaBind2 论文](https://doi.org/10.1016/j.isci.2020.100939)，同时注明使用本 GPU 加速改进版及其提交版本。[许可和第三方声明](THIRD_PARTY_NOTICES.md)
