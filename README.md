# BDS Dielectric Analysis

用于聚合物和离子材料宽频介电谱分析的 Codex skill，当前版本 **2.5**。

直接读取文件中的相对介电实部 ε′ 和损耗 ε″，核对电导率换算，生成导数介电谱，筛查候选弛豫峰与近似导电平台，并保留可追溯的处理记录。

## 功能

- 读取已支持的 FSM 多温度 TXT，或规范化长表 CSV。
- 计算 σ′ = 2πfε₀ε″，提供 S/m 与 S/cm 单位。
- 在实际频率网格上按频率跨度求导；支持可配置的点数窗口和统一 13 点主窗口。
- 保留负导数、端点空缺和原始采样点，检查频率缺口。
- 筛查并跟踪跨温度候选峰；在 ε″ 与 σ′ 图中标记同一候选频段。
- 输出数值表、PNG/SVG/PDF、配置、脚本、输入哈希及结果版本清单。
- 按已确认取值表绘制 σTτ–τ 主图或两栏对比图，保留温度标注和几何平均参考线。
- 输出乘积、离散程度、来源版本和 Origin 独立 XY 绘图数据；不重新选择 σ 窗口或峰。

本仓库仅包含通用代码、文档和合成数据测试，不包含个人路径、实验原始文件或具体项目的数据摘要。

## 安装

需要 Python 3.10 或更新版本。

```bash
git clone https://github.com/klu6163-droid/bds-dielectric-analysis.git
cd bds-dielectric-analysis
python -m pip install -r requirements.txt
```

作为 Codex skill 使用时，将仓库放入 Codex 的 `skills/bds-dielectric-analysis` 目录，并调用 `$bds-dielectric-analysis`。也可独立运行下列 Python 命令。

## 使用

按照 [配置与输入格式](references/input-and-run.md) 准备配置文件。配置中的相对路径按配置文件所在目录解析。

```bash
# 只检查文件，不生成图或结果目录
python -X utf8 -B scripts/run_bds.py --config bds_config.json --inspect

# 新建结果目录；已有输出目录会被拒绝
python -X utf8 -B scripts/run_bds.py --config bds_config.json --out results_new

# 在复用结果前检查来源和版本
python -X utf8 -B scripts/result_lock.py --run results_new --files all_data.csv candidate_peaks.csv candidate_IC_regions.csv
```

新分析仅使用文件中的 εr′、εr″，不从 Cp–D 重建。归档结果的历史复现需要明确标记原路线，不能混用不同版本的 σ 与 τ。

## σTτ–τ 主图与两栏对比

先按[选值和版本规范](references/result-version-and-sigma-tau.md)准备 manifest 绑定的逐温度选值表，再运行：

```bash
python -X utf8 -B scripts/plot_sigma_tau.py --run selected_run --table selected_values.csv --out product_plot_new --selection-note "Reuse the confirmed candidate sequence and conductivity windows; assignments remain unconfirmed." --layout main
```

`--layout both` 同时保存主图与两栏图；右栏默认相对本组 25 ℃ 归一化，可用 `--reference-temperature` 调整。缺少该温度时不替换或插值。完整参数和字段见[绘图运行说明](references/sigma-tau-plotting.md)。

输出 PNG/SVG/PDF、计算表、统计表和来源记录。Origin 交接文件包含各组独立 X/Y 列及参考线坐标；它们不是 OPJU，也不是可直接执行的 render plan。需要可编辑 Origin 项目时，按运行说明衔接可用的 `originpro-cli-plot` skill。

## 测试

```bash
python -X utf8 -B scripts/verify_numerics.py
python -X utf8 -B scripts/verify_sigma_tau.py
```

两组测试使用临时目录中的合成数据；乘积绘图测试检查公式、状态分组、缺失温度、混合版本拒绝、导出和原文件保护。既有谱测试覆盖单位、求导、13 点主窗口、缺口分段、平台/非平台筛查、候选峰跟踪、源文件保护和结果版本核验。

## 方法边界

候选峰不自动等于链段 α 弛豫，候选平台也不自动等于 σdc。未检出候选不代表相应过程不存在。窗口敏感性不是置信区间，σTτ近似恒定不能单独证明耦合机制。

详细说明：

- [Skill 入口](SKILL.md)
- [方法与科学边界](references/method.md)
- [结果版本与 σTτ–τ 图](references/result-version-and-sigma-tau.md)
- [σTτ 绘图运行说明](references/sigma-tau-plotting.md)
- [文献方法说明](references/literature-method.md)
