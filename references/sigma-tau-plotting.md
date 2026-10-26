# σTτ–τ 图：运行与 Origin 衔接

用于已经确认取值、需要主图或两栏对比的请求。先读[版本与科学边界](result-version-and-sigma-tau.md)。此入口不重新寻峰、不用乘积平直程度选择频段，也不负责给峰作物理归属。

## 输入

依赖 Python 3.10+、NumPy、Matplotlib。输入是 `result_lock.py` 验证通过的结果目录及其中一个已绑定的逐温度 CSV。原始候选峰表含多个峰时不能直接当成最终选值表；先按用户指定和谱证据选定候选，再另存为有来源记录的结果版本。

必需字段：

|字段|规则|
|sample、temperature_C|每个样品每个测得温度一行；重复行停止，不悄悄合并。|
|sigma_S_cm、f_start_Hz、f_end_Hz|已选 σ 和支持频段；不重新平均或移动窗口。|
|f_peak_Hz、tau_peak_s|核验 τ=1/(2πfpeak)；不从最高峰重新取值。|
|sigma_status|confirmed_dc / candidate_plateau / ac_reference / missing。|
|tau_status|assigned_segmental / unassigned / missing。|
|selection_reason|该行实际选择理由，不能写“让乘积更平”。|

原始来源行、求导窗口、边界警示、证据说明等其余列全部保留。confirmed_dc 必须有 sigma_evidence，assigned_segmental 必须有 tau_evidence；这些文本是记录，不是脚本自动验证的物理证明。缺失或非正 σ/τ 的测得温度仍留表，新增排除原因，不能为了上对数轴补正值。

已有 analysis_id、sigma_analysis_id、tau_analysis_id 时三个必须齐全且一致，匹配父 manifest 的 analysis_id；当父结果明确声明 science.selected_input_analysis_id 时按该绑定快照身份核验。若表未含任何身份列，三个字段从该已验证 manifest 的身份补齐，来源表必须已经绑定在此 manifest 中。不能给未核验的任意 CSV 填入 ID 绕过版本检查。

可选 comparison_group（优先）或 track_id 用于区分同一样品的候选序列。没有这两列时，以此次明确选择的序列记录为 selected，`--selection-note` 要写明选值版本及每个样品内的分组依据；这不等于确认同一物理机制。统计始终按 sample、sigma_status、tau_status、候选组分开，不能将不同身份混成一条参考线。此表只容纳每样品每温度一个选择；多条候选序列的比较分别制表运行。

## 命令

```powershell
python -X utf8 -B '<skill目录>/scripts/plot_sigma_tau.py' --run '<已锁定结果目录>' --table 'selected_values.csv' --out '<新输出目录>' --selection-note '沿用已确认取值表；每个样品内使用所选候选序列，物理归属未确认。'
```

- 默认 `--layout main`：只画左侧 σTτ–τ 主图。
- `--layout two-panel`：左侧主图，右侧随温度变化的相对乘积；`--layout both` 同时保存两种布局。
- `--reference-temperature 25` 指定右图归一化温度。某组没有该温度的有效值时，保留组记录并留空归一化序列，不改用其他温度、不插值。
- `--sample-order SampleA SampleB` 可指定图例/配色顺序，必须完整列出本表样品。不将具体样品名称或固定颜色硬编码为其他项目的要求。
- `--dpi 300` 设置 PNG 分辨率；SVG 保留文字，PDF 嵌入字体。
- `--artifacts-only` 仅用于明确要求的离线核查，结果注明未复核父版本原始来源。通常不使用。

T 使用 K，σ 使用 S/cm，τ 使用 s。主图默认双对数、空心散点、温度标注；至少两个有效点的组才画几何平均水平虚线。带 warning 标识的点加星号。均值线是描述性参考，不是拟合。单点/空组仍保留统计但注明无法评价恒定性。右图连线只表示测得温度顺序，已记录的缺失行用 NaN 断开。

已保存的 temperature_K、product_S_cm_inv_K_s 若与公式不一致会停止；符合容差的已有乘积原值保持，不为排版重新取点。原输入 CSV 另存快照，以便核对计算表中的新增列。

## 输出和检查

- `sigmaTtau_main.*` / `sigmaTtau_two_panel.*`：PNG、SVG、PDF。
- `sigmaTtau_values.csv`：保留全部输入行/来源列，补计算、有效性、分组及归一化状态；`descriptive_statistics.csv` 含几何均值、算术均值、样本标准差、CV、极差比、温度数及 τ 跨度。
- `origin_points_wide.csv`：每组独立一对 X/Y 列；`origin_reference_lines.csv`：每条参考线两端坐标；`origin_handoff.json`：真实列配对、颜色/符号、坐标类型、单位和标签来源。
- 输入快照、脚本、配置、父 manifest 和新 `result_manifest.json`。输出目录必须新建且不在父结果目录内，原始和之前结果不改。

运行后核验新 manifest，打开实际导出 PNG，检查标签、坐标、图例、星号和近邻温度标注。脚本成功不代表视觉检查已完成。窗口敏感性不是重复测量的置信区间，不制造误差条。

## 用户要求 Origin 可编辑图时

读取本机可用的 `originpro-cli-plot/SKILL.md`，沿用已明确的图型、版本和样式选择。`origin_handoff.json` 是数据交接说明，**不是 render plan，不直接交给渲染器，也不代表 OPJU 已生成**。

1. 从 wide CSV 读取每组独立 X/Y 对，画原生散点，横轴 τ、纵轴 σTτ；不要让不同样品共用第一组的 τ 列。用已导出的两端坐标添加独立细点线参考。温度标签来自计算表，不靠 X 排序猜测。
2. 新图用该 skill 支持的 scatter 流程生成真实 understand/plan；若该路由只接受共享 X，分别配置独立数据或实现经检查的原生模板，不能错误套用。原生 executor 只替换第一工作表/图页/图层；多页项目先核对保存后的顺序，不能靠激活顺序推断。
3. 图上保持 σ 的实际身份和 τ 的归属状态；确认链段归属前用 τpeak。以 Origin 原生上/下标处理 cm⁻¹、τpeak，生成可编辑数据和对象。程序化设置线宽若用 LabTalk，`-wp` 表示点值，不把 `-w` 当作 pt。
4. 依 Origin skill 的来源绑定交付流程保存、独立重开核对数据、导出所有图页并逐张查看最终 PNG。科学坐标不得因样式修改而改变。

若还要求实部、虚部、导数证据图：读取与取值表同一版本的 ε′、ε″和导数数组，按样品分别绘制各温度；虚部标记已选 σ 支持点并注明频段，导数标记当前采用的峰。保留负导数、边缘缺失和未选响应；导数纵轴可用线性以完整展示。该乘积脚本不重算或导出证据谱，需从已绑定谱表完成，不能从 σTτ 点反推证据曲线。
