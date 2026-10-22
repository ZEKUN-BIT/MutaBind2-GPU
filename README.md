# MutaBind2-GPU：GPU 加速改进版

MutaBind2-GPU 用蛋白质复合物结构预测氨基酸突变引起的结合亲和力变化（ΔΔG）。这是基于 [MinghuiLab/MutaBind2](https://github.com/minghuilab/MutaBind2) 的**非官方 GPU 加速改进版**，由 [ZEKUN-BIT](https://github.com/ZEKUN-BIT) 维护。

这个版本沿用原版的随机森林模型，主要改进计算速度和批量任务的运行方式，没有重新训练预测模型。原方法、模型和上游版权声明均予以保留；本项目不属于原作者的官方发布。

## 这个版本改了什么

- 为 CHARMM 中受支持的 PB SOR 网格迭代、边界电势求和，以及 REEN 探针计算和网格恢复加入 CUDA 加速。Fortran 接口会检查运行条件，不支持时回退到 CPU，并在日志中说明。
- 调整 VDW 内存分配，加入局部解析 SASA 计算和野生型（WT）结果缓存；FoldX 最多同时运行两个独立进程。
- 提供 Python 3 命令行入口，支持独立单点替换、饱和扫描和断点续算。保留原始链名、残基编号及结构断链的映射关系，并保存计算中间结果。
- 调整 PROVEAN 1.1.5 对新版 BLAST 序列标识符的处理，加入按序列缓存分数的功能，并提供可选的完整序列计算流程。

目前验证过的 GPU 扫描流程适用于**独立单点替换和饱和扫描**。原有的多点组合突变接口仍然保留，但组合突变的 GPU 流程尚未验证。PROVEAN 和最终的随机森林预测仍在 CPU 上运行。

## 实测能快多少

下面是在 Linux x86-64、RTX 5070 Ti 和 CHARMM 51 环境中，使用固定的 1CSE 野生型与 L45P 突变体结构得到的结果：

| 测试对象 | CPU：VDW + PB | GPU：VDW + PB | 能量对照 |
|---|---:|---:|---|
| 野生型（WT） | 50.19 秒 | 4.51 秒 | 日志打印精度内一致 |
| I:L45P | 46.58 秒 | 4.51 秒 | 日志打印精度内一致 |

这些时间只对应 **VDW + PB 计算阶段**，不含锁等待时间，不能直接当作整次预测的加速比。完整任务还包括结构处理、PROVEAN 等步骤，实际耗时取决于输入和运行环境。

CPU/GPU 对照使用相同的坐标、力场输入、随机种子和网格。网格间距为 0.4 Å，各组分均从冷启动开始求解，保留 2×10⁻⁶ 的残差阈值和收敛检查。详细设置见 [验证说明](docs/VALIDATION.md)，原始能量对照见 [JSON 记录](validation/1cse_energy_comparison.json)。这些结果用于验证计算实现的一致性，不代表预测的生物学准确率有所提高。

## 安装

运行环境为 Linux、Python 3.11 及以上版本。使用 GPU 加速还需要 NVIDIA GPU 和 CUDA 编译工具链。当前 GPU 接口针对经过哈希校验的 CHARMM 51 源码基线构建，CHARMM 需从授权渠道获取。

先下载代码并安装 Python 部分：

```bash
git clone https://github.com/ZEKUN-BIT/MutaBind2-GPU.git mutabind2-gpu
cd mutabind2-gpu
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[test]'
cp config/tools.example.json config/tools.json
```

在 `config/tools.json` 中填写本机的程序和 NR 数据库路径，再按 [安装指南](docs/INSTALL.md) 编译 GPU 库。配置完成后，检查环境并运行测试：

```bash
./mutabind2 doctor --scope gpu
python -m pytest -q
```

Python 测试本身不需要 GPU 或专有软件。要完成一次实际预测，还需要另外安装 VMD、CHARMM、NAMD、DSSP、FoldX、PROVEAN、BLAST+、CD-HIT、NR 数据库，以及 R 和 `randomForest`、`forestFloor` 包。安装 Python 包不会自动安装这些依赖。

## 跑一个示例

仓库附带 PDB **1CSE** 示例，两个结合伙伴分别为 E 链和 I 链，示例突变为 I:L45P。

先生成任务计划，再开始计算：

```bash
./mutabind2 scan run \
  --pdb examples/1cse.pdb --partner1 E --partner2 I \
  --mutations examples/mutations.tsv \
  --output projects/1cse_scan --plan-only

./mutabind2 scan run --resume --output projects/1cse_scan
```

突变表用制表符分隔，格式如下：

```text
chain	mutation	id
I	L45P	1
```

每一行代表一次独立的单点替换。如果要做饱和扫描，可以把 `--mutations examples/mutations.tsv` 换成 `--sites 'I:45-48'`。程序会为结构中有坐标的每个位点生成 19 种非野生型替换，不会补建缺失坐标。同一次扫描可以包含同一结合伙伴内的多条链；如果要分别突变两个结合伙伴，请拆成两个任务。

### 查看结果和暂停任务

完成的预测保存在输出目录下的 `completed_predictions.tsv`。ΔΔG 的单位是 kcal/mol：正值表示结合减弱，负值表示结合增强。结果中的 PROVEAN 分数反映序列变异影响，不是结合 ΔΔG。

查看进度：

```bash
./mutabind2 scan status --output projects/1cse_scan
```

如果需要暂停，可创建下面的标记文件，程序会在当前批次结束后停止：

```bash
touch projects/1cse_scan/STOP_AFTER_BATCH
```

继续运行前，先删除 `STOP_AFTER_BATCH`，再执行上面的 `--resume` 命令。修改代码、外部工具、模型或数据库后，请换一个新的输出和缓存目录，避免混用旧结果。

## 使用前需要了解的限制

默认 PROVEAN 流程按连续结构片段计算，分数可能与完整蛋白序列的结果不同。如果没有找到同源序列，程序会退回到仅用查询序列自身计分，目前最终结果表还没有为这种情况单独标注质量状态。PROVEAN 缓存也不会校验数据库和工具版本；更换它们后，需要使用新的缓存目录。

其他适用范围和边界情况见 [已知限制](docs/KNOWN_LIMITATIONS.md)。

## 来源、引用与许可

原始项目为 [MinghuiLab/MutaBind2](https://github.com/minghuilab/MutaBind2)，采用 MIT 许可，本仓库保留其版权声明。GPU 加速和运行流程的改进由 ZEKUN-BIT 维护。

使用本工具开展研究时，请引用原始论文：

> MutaBind2: Predicting the Impacts of Single and Multiple Mutations on Protein–Protein Interactions. iScience (2020). [DOI: 10.1016/j.isci.2020.100939](https://doi.org/10.1016/j.isci.2020.100939)

同时注明使用了本 GPU 加速改进版及对应的提交版本。引用信息见 [CITATION.cff](CITATION.cff)，许可说明见 [LICENSE](LICENSE) 和 [第三方声明](THIRD_PARTY_NOTICES.md)。

`native/` 目录提供 CUDA 源码和最小化的 CHARMM 集成补丁，不包含完整 CHARMM 源码或二进制程序。`toolchain/vendor/provean-1.1.5/` 中单独构建的 PROVEAN 源码保留 GPL 许可和修改说明。仓库保留了原模型与必要资源，但不附带专有软件二进制、NR 数据库或完整运行环境；外部工具和力场资源仍遵循各自的使用条款。
