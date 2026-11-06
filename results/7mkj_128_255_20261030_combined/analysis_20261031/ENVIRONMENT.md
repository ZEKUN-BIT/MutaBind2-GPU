原始分析使用的软件版本逐字节保存在 `historical_environment.freeze`，也保留在完整结果归档中。

当前 requirements 清单将 Pillow 更新为 13.0.0，以修复 GitHub 报告的旧版本安全问题。它用于新建分析环境，不应解释为原图生成时的版本。原始结果、图表、分析脚本和校验记录没有重算或改写。
