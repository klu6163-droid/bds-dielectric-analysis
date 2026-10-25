# 输入与运行

依赖Python 3.10+、NumPy、SciPy、Matplotlib。下列命令中的python应指向安装了这些依赖的解释器；脚本路径按实际安装位置调整。

## 配置

```json
{
  "datasets": [
    {"sample": "SampleA", "path": "raw/sample_a.csv", "format": "long_csv"},
    {"sample": "SampleB", "path": "raw/sample_b.csv", "format": "long_csv"}
  ],
  "formats": ["png", "svg", "pdf"],
  "protect_paths": ["raw"],
  "screening": {"max_frequency_gap_decades": 0.15},
  "derivative": {
    "mode": "frequency_span",
    "spans_decades": [0.25, 0.35, 0.50],
    "primary_span_decades": 0.35,
    "min_points": 5,
    "min_points_each_side": 2,
    "max_gap_decades": 0.15
  },
  "tracking": {
    "enabled": true,
    "max_shift_decades_per_C": 0.12,
    "max_temperature_gap_C": 10,
    "ambiguity_margin_decades": 0.10,
    "use_prediction": true
  }
}
```

相对路径按配置文件目录解析。仅支持原始介电列；可省略conversion，或明确设为exported_epsilon。旧Cp–D模式及两批配置字段会被拒绝，不能静默改用另一种输入。无需用户再选路线。

用户指定统一13点主窗口时，derivative部分改为：

```json
{
  "mode": "legacy_points",
  "windows_points": [7, 9, 13],
  "primary_window_points": 13,
  "max_gap_decades": 0.15
}
```

7和9点仅用于敏感性比较，主导数、候选峰和τpeak均取13点结果；“稳定候选”仍要求跨敏感性窗口通过匹配。窗口列表须含至少两个互不重复的奇整数且均≥5，主窗口必须在列表中。未指定的旧legacy_points配置仍以9点为主。点数与频宽模式不能因为默认值不同而在同批比较中混用。

screening.max_frequency_gap_decades控制平台扫描；derivative.max_gap_decades控制求导。两者分别记录。改变任一参数都产生新分析版本。

## 输入格式

- `fsm_txt`：UTF-8-SIG，已知FSM分区格式。使用εr′、εr″分区；温度列、索引和频率须一致。其他分区原样保留在表内，仅作来源记录。
- `long_csv`：每行一个温度下的频率点，必需temperature_C、frequency_Hz、epsilon_prime_export、epsilon_doubleprime_export。ε列须为无量纲相对介电值。其他列作为原始记录保留。
- 每个文件对应一个样品；CSV若有sample列，须与配置一致。历史样品映射不能套用到新实验。
- 面积与厚度不是必需的计算输入。FSM头部中可识别的几何值保留为元数据；如显式覆盖已有值，需geometry_note记录依据。
- 频率应为严格递增的正值，温度为℃。不自动排序、删点或去重。缺少必需介电列，或有非正/非有限值时停止自动分析并报告；先由用户决定是否另做清理。

## 命令

```powershell
python -X utf8 -B scripts/run_bds.py --config bds_config.json --inspect
python -X utf8 -B scripts/run_bds.py --config bds_config.json --out results_new
```

inspect只读文件结构，不建目录、不绘图；正式运行再验证全部介电值。输出目录必须不存在。analyze_bds.py也支持相同参数；run_bds.py是简短入口。

## 输出与核验

输出为一套直接值分析：各样品与温度的ε′、ε″、σ′、导数谱PNG/SVG/PDF，候选峰与IC表，频宽支持、IC频率缺口与阈值敏感性表，跨温度轨迹及链接审计，原始列与派生量总表，配置、脚本、输入SHA-256、result_manifest.json和ZIP。

conversion_branch=exported仅标识原始介电列来源。不会生成第二批、Cp–D派生列、两路对照或user_selection.json。源文件里的Cp和D若存在，仍作为未参与计算的原始记录保存；删除的是处理分支，不是原始数据。

重用结果前执行：

```powershell
python -X utf8 -B scripts/result_lock.py --run results_new --files all_data.csv candidate_peaks.csv candidate_IC_regions.csv
```

正常检查核验原文件和配置，以及已保存结果表、脚本等的哈希；程序不重算谱。仅在原文件不可用、确需离线核查已归档结果时显式使用--artifacts-only，并说明原文件未在当前环境复核。该选项仍核验结果表哈希。没有manifest的历史结果不能伪装为当前路线，按[版本规则](result-version-and-sigma-tau.md)在新目录另建历史快照。

运行scripts/verify_numerics.py检查单位、原值保真、主窗口选择、频率缺口、筛查、跟踪、缺列拒绝、旧配置拒绝、结果锁定与原文件保护；元数据用skill-creator的quick_validate.py验证。实测回归核对ε值逐点不变、派生结果与此前同路线一致。
