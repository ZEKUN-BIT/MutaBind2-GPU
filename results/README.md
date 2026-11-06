# 计算结果与完整归档

本目录保存 2026-11-06 同步时的预测表、输入记录、分析图表及精简的验证证据。完整原始结果、批次科学中间文件、工作区和验证电势文件见 [v2.1.0 Release](https://github.com/ZEKUN-BIT/MutaBind2-GPU/releases/tag/v2.1.0)。

| 结果 | 完成预测 | 说明 |
| --- | ---: | --- |
| [128–255 合并结果](7mkj_128_255_20261030_combined/completed_predictions.tsv) | 1,463 | 合并原始、warp、cells 与 GPU 四个计算协议；每行保留来源 |
| [256–374 结果](7mkj_256_374_20261102/completed_predictions.tsv) | 2,261 | 119 批，全部完成 |
| [早期 7MKJ 结果](7mkj_rpod/completed_predictions.tsv) | 668 | 早期按优先位点完成的预测；不代表整个计划全部计算 |

这些任务之间可能有重复位点/突变，不应简单相加作为独立样本量。128–255 的四个来源目录也不应再次与合并表累计。

- [128–255 分析](7mkj_128_255_20261030_combined/analysis_20261031/)
- [256–374 分析](7mkj_256_374_20261102/analysis_20261106/report_zh.md)
- [早期分析](7mkj_rpod/analysis_20261022/)
- [索引](index.json)

## 下载和校验完整结果

安装 `gh` 和 `zstd` 后，在一个新的目录下载附件：

```bash
gh release download v2.1.0 --repo ZEKUN-BIT/MutaBind2-GPU
sha256sum --check SHA256SUMS
mkdir restored
cat projects.tar.zst.part* | zstd -d | tar -xpf - -C restored
cat native-validation.tar.zst.part* | zstd -d | tar -xpf - -C restored
```

解压所有八组结果后，再解压 `result-links.tar.zst`：

```bash
tar --zstd -xpf result-links.tar.zst -C restored
```

这个补充包将 579 个结果文件/目录别名恢复为归档内部的相对符号链接，避免依赖原机器绝对路径；映射见 `result_links.json`。

其余附件中，单个 `.tar.zst` 可用 `tar --zstd -xpf 文件名 -C restored` 解压；有 `.part000` 等后缀的组按名称顺序拼接后解压。具体分卷列表见 `bundle_summary.json`。恢复目录内 `_manifest/` 记录每个文件的 SHA-256 与原始大小，同样的清单也作为独立附件提供。部分内容相同的验证文件在包内以硬链接存储，解压后各路径的文件内容与原始记录一致。

结果包包括 `projects`、`workspaces`、`native-validation`、`validation`、`legacy-predictions`、`legacy-workspaces`、`toolchain-results` 和 `toolchain-validation` 八组。`excluded_non_results.json` 列出未公开的原生二进制、构建对象、完整 CHARMM 派生源码和原始绝对符号链接；它们均保留在本机完整安装中。结果之间的链接由上述补充包恢复，外部软件资源链接仍需用户自己的安装。软件安装资源不属于预测结果。

结果中的绝对路径、时间和协议哈希是历史计算证据，保留原字节。它们不是当前机器的路径配置。旧目录 `/home/zek/mutabind2-tool` 的文件已迁入 `/home/zek/MutaBind2-GPU`；早期外部结果在本机 `projects/7mkj_rpod` 与 `projects/legacy_workspaces` 下有完整副本。使用迁移后的软件创建新任务；不要把历史停止的过渡扫描直接重启，否则会重复计算已并入合并表的突变。

两份分析目录的可安装依赖清单已将 Pillow 更新至 13.0.0。原始分析环境版本单独保存在各目录的 `historical_environment.freeze` 和原始归档中；相关说明见各目录的 `ENVIRONMENT.md`。
