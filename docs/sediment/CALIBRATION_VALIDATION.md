# CoLM202X-SEDIMENT：调参与验证工作规程

> **定位**：服务于亚马孙河道泥沙模型的「CoLM 运行 → OpenBench 评价 → 物理过程诊断 → 设计下一轮参数实验」闭环。模型原理、物理公式和完整参数字典请先阅读 [`SEDIMENT_MODEL.md`](SEDIMENT_MODEL.md)。  
> **版本**：v0.1，2026-10-15。  
> **文档核对基线**：CoLM202X_SYSU [`32ea6ce72b0a8ced75f5859d47c3d6cb3632cd07`](https://github.com/jellywzx/CoLM202X_SYSU/commit/32ea6ce72b0a8ced75f5859d47c3d6cb3632cd07)；zhongwangwei/OpenBench [`f2d841609d3c214a8726e1abde5f4ca177cefc60`](https://github.com/zhongwangwei/OpenBench/commit/f2d841609d3c214a8726e1abde5f4ca177cefc60)（3.0.8）。这**不是**历史运行时已证明使用的两个版本。  
> **证据标记**：**[S]** 已在上面固定源码中核实；**[O]** 历史对话/实验报告的观察，尚未重新读取原始产物；**[H]** 待测试物理机制假说；**[T]** 需补充实际配置/数据/运行结果验证。  
> **不可混淆**：调参结果不等于独立验证；模拟全时段统计不等于有观测的配对日统计；数值保护上限不等于物理输沙能力参数。

## 0. 一页式操作流程

1. **锁定基线**：保存 CoLM/OpenBench SHA 与差异、编译选项、可执行文件哈希、main/SEDIMENT/history namelist、forcing、routing 文件、初始状态与 reference 版本。
2. **验证数据合同**：核对 Q/SSC/SSL 单位和逐粒级计算、参考 QC、站点—河网对应、时间标签及有效样本数。
3. **基准复现**：先使用完全不变的物理参数重跑；不同实验一律从同一或物理一致的初态和相同 spin-up 方案开始。
4. **运行与评价**：OpenBench 对每个 `station × variable` 计算仅基于其有效配对日的指标；另输出全模拟期诊断。
5. **物理归因**：同步分析 `sedinp`、`Es_raw/eff`、`D_raw/eff`、四类 `CapDep`、`SSC_precap_max`、悬沙与河床储量及水动力。
6. **提出假说**：预写干预参数、期待现象和反证条件；优先一次只改一个参数，或者一组明确耦合的参数。
7. **选择下一轮**：综合守恒、Q、SSC、SSL 的误差、粒级和空间一致性及独立验证表现；接受/拒绝/重测都有记录。

**停止门禁**：若质量守恒异常、未收敛/崩溃、单位/时间错误、站点严重错配或关键指标无有效配对值，暂停基于性能分数的物理调参。不得因报错而随意关闭守恒检查。

## 1. 研究范围、阶段及验证与率定分离

- **研究区**：亚马孙流域；模式为 CoLM `GridRiverLakeFlow` 中的 SEDIMENT 粒子模块，默认按 clay/silt/sand 三粒级模拟。
- **主观测变量**：Q = `Discharge_For_Sediment`，SSC = `Suspended_Sediment_Concentration`，SSL = `Suspended_Sediment_Load`。
- **近期问题**：已有部分站点 Q 尚可，但 SSC 和 SSL 高出观测几十至上千倍，并且 SSC 在若干格网反复接近约 26,500 mg/L 上限 [O]。
- **阶段 A（baseline）**：追溯历史 2001–2002 运行，固定所有输入、输出和配对规则。
- **阶段 B（归因）**：排除 Q、站点匹配、单位、时间偏移和数值守恒问题，查清强迫/产沙/起悬/沉积/限幅主导关系。
- **阶段 C（筛选）**：进行物理可解释的单参数或分组敏感性实验。
- **阶段 D（率定）**：在事先声明的训练集上确定候选参数，报告多站多指标折中。
- **阶段 E（独立验证）**：在不参与选参的年份和/或站点上检验，并检验不同季节、流量范围与极端事件。

**原则**：不能连续使用同一 4 站 2001–2002 的全部观测日优化参数，然后把这批同样的日期称为独立验证。若单站观测日过少，应扩大可用时段和站点，或先明确称为“诊断/预率定”，不夸大独立性。

## 2. 历史 baseline 候选：必须保留的不确定性

### 2.1 2026-10-14 历史报告 [O]

之前提交的 OpenBench 报告描述 `Amazon_sed_0p25` 相对 `OpenBench_Sediment_Daily` 的 2001–2002 日尺度、4 站评估：三项变量成功，`Overall_Score` 均值 Q ≈ **0.492**、SSC ≈ **0.097**、SSL ≈ **0.092**。这些是当时项目的综合评分，**不是** KGE/NSE，也不是已与当前代码 SHA 对齐的结果。

| 站点 ID | 历史报告有效配对天数 | SSC 模拟均值 / 观测均值 | SSL 模拟均值 / 观测均值 | 名称核实状态 |
|---|---:|---:|---:|---|
| `SED006136` | 53 | 13.1× | 6.7× | [T] 待复核 |
| `SED006139` | 42 | 1104.2× | 2190.4× | [T] 待复核 |
| `SED006140` | 68 | 49.7× | 36.6× | [O] Porto Velho |
| `SED006141` | 58 | 113.6× | 73.7× | [O] Fazenda Vista Alegre |

**使用边界**：天数和倍数来自历史会话总结，未在编写本文件时重新读取配对 CSV。应分别重新统计 `n_Q/n_SSC/n_SSL/n_joint`；上表的 `n` 暂不能保证三变量完全相同。倍数指同一配对日期集合上的 `mean(sim)/mean(obs)`，不可替换为 `mean(sim/obs)`。

### 2.2 浓度平台：全时段与配对日不能混用 [O]

此前更进一步的测试报告称：

| 站点 | 全部模拟日接近浓度上限 | SSC 配对日接近浓度上限 |
|---|---:|---:|
| Porto Velho / `SED006140` | 665/730（91.10%） | 62/68（91.18%） |
| Fazenda Vista Alegre / `SED006141` | 390/730（53.42%） | 26/58（44.83%） |

“接近上限”的**数值容差和原始字段仍待确认**；这些数字仅为历史高估现象提供线索，不能视为已经核实的实际 CapDep 事件次数。后续应通过 `f_ssc_precap_max`（限幅前最大值）与 `f_capdep_*`（实际限幅转移）交叉验证。

### 2.3 源码 Debug 四站和 OpenBench 四站不是一回事

CoLM `CoLMDEBUG` 的固定诊断站是 **Serrinha、Porto Velho、Manacapuru、Obidos**；模型按照最近 unit-catchment 网格中心定位，并不检验水文连通或上游面积 [S]。[源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L111-L127)

OpenBench 历史 4 站集合包括 `SED006141 = Fazenda Vista Alegre` [O]，**不得**因为同为“四个站”就把两份名单按顺序对应。每站必须保存独立映射 `station_uid → CaMA 坐标 → model grid/unitcatchment ID`。

### 2.4 目前尚未达到“baseline 已复现” [S/T]

- 当前仓库 [`run/Amazon_0p25.nml`](../../run/Amazon_0p25.nml) 的示例模拟仅结束在 **2001-12-31**，而历史 OpenBench 评价为 **2001–2002 两年**；实际提交运行的 namelist 未在本文件中确定。
- 当前 [`run/sediment_parameter.nml`](../../run/sediment_parameter.nml) 和 [`run/standard_sediment_parameter.nml`](../../run/standard_sediment_parameter.nml) 粒径、浅水阈值和活动层厚度不一致，不能互相代替。
- 历史结果尚缺精确的 CoLM/OpenBench SHA、forcing/routing/restart 版本、有效站点映射、原始评估 CSV 与变量时间对齐策略。
- 曾有模型和 OpenBench 更新改变浓度限制、泥沙 history、hydro handoff 或时间偏移；旧结果不能自动当作现行 master 的性能。

先称为 `baseline_candidate_20261014`，直到以上内容补齐且可复算，才建立 `BASE_2001_2002` 正式基线文件。

## 3. 模拟与观测质量门禁（正式评估前）

### 3.1 模型和输入文件核对

| 检查项 | 必须记录/检查 | 发生问题时 |
|---|---|---|
| CoLM 实际代码 | Git HEAD、dirty diff、`TRACER`/`GridRiverLakeFlow`/`CoLMDEBUG` 宏、编译日志、可执行文件哈希 | 避免“修改代码后仍运行旧 exe” |
| 主 namelist | 实际绝对路径、SHA256、流域范围、起止年月、spin-up、陆面步长、routing dt | 对照真正运行的任务脚本而非仓库示例 |
| 泥沙参数 | `DEF_TRACER_PARAM_FILES` 的真实路径、`SEDIMENT EFFECTIVE PARAMETERS` 日志 | 区分实验文件和标准文件 |
| 静态河网 | `DEF_UnitCatchment_file` 的 NetCDF 版本、`sed_n`、`sed_frc`、`sed_slope`、`numucat` | 检查是否把全球 unitcatchment 当作区域路由 |
| forcing | 强迫名称、版本、日内时间标签、降水缺测/负值/单位 | 产沙模型对降雨非线性敏感 |
| history/restart | 频率、`groupby`、`one/block`、`seddiag`、重启动参数元数据与守恒状态 | 不复用与新参数物理不兼容的 restart |
| 观测 | 站点、QC、参考版本、Q/SSC/SSL 单位、覆盖期和缺测说明 | 缺测不是零 |

**参考 QC 推荐合同**：当源数据使用 0=good、1=derived、2=suspect、3=bad、9=missing 时，核心基准优先保留 0–1；更宽松/严格的标记做独立敏感性分析。不要假设 OpenBench 的 `reference_catalog.yaml` 自动读取 Q/SSC/SSL flag；需先验证 `OpenBench_Sediment_Daily_full.nc` 是否已按质量标志过滤，并保存处理依据。

### 3.2 当前 OpenBench 参考注册 [S]

[`reference_catalog.yaml`](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/data/registry/reference_catalog.yaml#L4296-L4388) 中包含：

| Dataset | 频率 | 字段 | 注册表年份（不等于各站实际观测覆盖） |
|---|---|---|---|
| `OpenBench_Sediment_Daily` | Day | `discharge`, `ssc`, `ssl` | 1948–2025 |
| `OpenBench_Sediment_Monthly` | Month | 同上 | 1938–2000 |
| `OpenBench_Sediment_Annual` | Year | 同上 | 1912–2020 |

三套数据都是 `stn`，并用 `Discharge_For_Sediment`、`Suspended_Sediment_Concentration` 和 `Suspended_Sediment_Load` 评价。Daily 参考目录注册为 `OPENBENCH_REF_ROOT/Station/Water/Sediment/Daily`。若要拿 2001–2002 测试 Monthly，必须先确认 Monthly 的实际数据和注册年份是否足够，而不是因为 Daily 覆盖这两年就默认 Monthly 也覆盖。

### 3.3 站点—河段的水文一致性门禁 [S/T]

当前 [`station_matcher.py`](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/data/station_matcher.py#L57-L108) 对 `cama_allocation` 的重要规则：

- `grid_res=0.25°` 对应 CaMA `15min`；
- 最小 `upstream_area=3000 km²`；
- 最大 `abs(cama_alloc_err_15min)≤0.2`；缺测的 allocation error 会使站点被剔除；
- 网格提取 [`_processing_station_extract.py`](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/data/_processing_station_extract.py#L76-L146) 使用已分配位置抽取最近网格，默认容差为模拟网格分辨率的一半（0.25° 时为 0.125°）。

匹配时仍要保存并人工审核：`station_uid`、站名、观测坐标/流域面积、`cama_lon/lat_15min`、`cama_alloc_err_15min`、CoLM grid/`ucat`、模拟上游面积、两者相对误差、是否同河段/连通、支流汇合处歧义、排除理由。**CaMA allocation 通过过滤不等于 CoLM 的河道位置一定水文正确**。

### 3.4 时间标签必须审计 [S]

当前 3.0.8 OpenBench 根据**实际选中的原生变量优先、评价项名称次之**解析模型时间偏移：[`_processing_time_adjustments.py`](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/data/_processing_time_adjustments.py#L72-L148)。

- `f_discharge`/`v_discharge` 在当前 CoLM2024 profile 中有原生 `0` 偏移特例；
- CaMa 风格 `outflw` 保留其周期末标记处理（常为 −1 day/月）；
- 泥沙 `f_sedcon_*`/`f_sedout_*` 在模型 profile 中属于计算项，可能按评价项规则产生偏移，**不能简单假定和 `f_discharge` 完全一样**。

正式运行前必须检查 `raw history time → OpenBench 处理后 time → 观测配对 time` 三列，至少抽查三到五个已知日期；查看本地 YAML 的 `time_alignment`、选中变量名与处理日志。一个日偏移足以改变稀疏 SSC 观测的 KGE/NSE/相关。

### 3.5 Q、SSC、SSL 单位及重复转换 [S]

[`model_catalog.yaml`](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/data/registry/model_catalog.yaml#L1408-L1419) 对 CoLM2024 的泥沙变量先做：

```yaml
Suspended_Sediment_Concentration:
  compute: ds.sum_prefix('f_sedcon_', 3) * 2650
  varunit: kg m-3
Suspended_Sediment_Load:
  compute: ds.sum_prefix('f_sedout_', 3) * 2650
  varunit: kg s-1
```

然后 [`unit.py`](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/data/unit.py#L291-L311) 把 `kg m-3` ×1000 换为 `mg L-1`，`kg s-1` ×86.4 换为 `t d-1`。独立校验公式（`rho_s=2650 kg/m³`，`nsed=3`）：

```text
SSC [mg/L]  = (f_sedcon_1+f_sedcon_2+f_sedcon_3) × 2,650,000
SSL [t/day] = (f_sedout_1+f_sedout_2+f_sedout_3) × 228,960
```

`f_sedcon_k` 为固体体积分数，`f_sedout_k` 为**有向悬浮泥沙**体积输运率；`f_bedout_k` 是推移质，**不要**在评价悬浮输沙量时混加。平均 `Q ×` 平均 `SSC` 一般不等于期间平均泥沙通量，因为输运、存储和时间协方差有差异。检查 CoLM 实际运行的 `grain_density` 与 `nsed`，若与默认 2650/3 不同，必须相应修改映射。

**验收**：在一个站点的正确网格及日期，直接从原始 NetCDF 复算 Q/SSC/SSL，并与 OpenBench *preprocessed* 结果核对，不允许已经换为 mg/L 后再次乘 1000。

## 4. 统一对齐样本、指标和报告方式

### 4.1 三种不得混用的日期集合

| 日期集 | 定义 | 用途 |
|---|---|---|
| `full_simulation_period` | 全部成功模拟日期；2001–2002 日序列共 730 日 | 全时期 SSC 平台比例、水深、Es/D、床沙库存等 |
| `paired_variable_dates` | 对每一变量，`sim/ref` 同日期且双方有效的交集 | 此变量的所有配对误差、相关、KGE/NSE |
| `paired_joint_dates` | 同日 Q、SSC、SSL 的 `sim/ref` 均有效 | Q 对 SSC/SSL 的共同驱动归因和联合分析 |

每站必须保存 `n_sim/n_obs_valid/n_Q/n_SSC/n_SSL/n_joint`；不能把 `full_period` 的模拟均值除以 `paired_days` 的观测均值。报告 full 和 paired 诊断可揭示“稀疏观测是否碰巧抽中了浓度平台期”。

### 4.2 推荐的评价指标

| 变量 | 优先指标 | 必须补充 |
|---|---|---|
| Q | `KGE`、`NSE`、`percent_bias` | 配对 n、Pearson r、均值比、丰枯季水文量级 |
| SSC | `percent_bias`、`rSpearman`、`RMSE` | 配对 n、模拟/观测均值比、近限幅频率、正值样本上的 `KGEln`（可选） |
| SSL | `percent_bias`、`KGE`、`NSE` | 配对 n、均值比、日峰值/时滞、Spearman 和季节性 |

OpenBench 当前 [`metrics.py`](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/core/metrics.py) 的重要约定：

- `KGE` 最优为 1，**负值没有固定下界**；同时受相关 `r`、标准差比、均值比控制。`KGESS=(KGE+0.41)/1.41` 是基准转换，不是更“准确”的独立物理量。
- `NSE=1−sum((sim−obs)^2)/sum((obs−mean_obs)^2)`，可以很负；观测常数时无法定义。
- `percent_bias=100×sum(sim−obs)/sum(obs)`，正为高估，不应简单裁剪到 ±100%。在 sim、obs 均非负、观测总和正的限定下，最低为 −100%，向上无界。
- `rSpearman` 检查排名一致性，不能代表绝对量级；常数序列、缺少有效样本时可能为 NaN。
- `KGEln` 按当前实现对 KGE 中的标准差比与均值比取对数，需严格正的比值；并非任意“先 log(SSC) 再调用 KGE”。
- 指标 NaN 要保留并写明原因；不能把 NaN 改成 0 来纳入总评分。

### 4.3 统一偏差定义

对于变量各自的 `paired_variable_dates`：

```text
paired_mean_ratio = mean(sim_paired) / mean(obs_paired)
paired_PBIAS[%]   = 100 × sum(sim_paired - obs_paired) / sum(obs_paired)
paired_log_bias   = log10(paired_mean_ratio), if ratio>0
```

`mean(sim/obs)` 与 `mean(sim)/mean(obs)` 不相同；对于接近 0 的 SSC/SSL 观测逐日比值可能爆炸。对正值且符合 QC 的日期可再报告 `median(sim/obs)`、分位数和对数偏差。不要只报告一个所有站点等权的 `Overall_Score`；至少同时报告**逐站全表、中位数、最差站点及有效样本数**。

### 4.4 运行 OpenBench：锁定本地安装与配置

下面假设已进入实际要运行的 OpenBench checkout，且有 `.venv`；这是**命令模板**，不是已经执行的服务器作业：

```bash
cd /absolute/path/to/OpenBench
git rev-parse HEAD
git status --short
.venv/bin/python -c "import openbench; print(openbench.__file__)"
.venv/bin/openbench check /absolute/path/to/sediment.yaml
.venv/bin/openbench run /absolute/path/to/sediment.yaml --dry-run
.venv/bin/openbench run /absolute/path/to/sediment.yaml --cores 1
```

**[S]** 当前 `run` 命令还有 `--variable`、`--resume`、`--comparison-only`、`--force`；`--resume` 需要兼容预处理产物，并不能与 `--comparison-only` 或 `--force` 同用。[CLI 源码](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/cli/run.py#L18-L147)

运行前检查：本地配置选中的 model profile、`grid_res=0.25`、station_matching、`time_alignment`、`unified_mask`、`sim/ref root`、三个 sediment evaluation item、比较图以及缺测处理。若运行不产生有限的 KGE/相关/NSE，先排查缺测、常数/异常序列和配对样本，不要直接调物理参数。

## 5. 必须输出的过程诊断及其解释

### 5.1 history 变量

在 `run/history_sediment.nml` 或实际被 `DEF_HIST_vars_namelist` 引用的文件中：

```fortran
DEF_hist_vars%discharge = .true.
DEF_hist_vars%riv_height = .true.
DEF_hist_vars%riv_veloct = .true.
DEF_hist_vars%floodarea = .true.
DEF_hist_vars%floodfrc = .true.
DEF_hist_vars%sedcon = .true.
DEF_hist_vars%sedout = .true.
DEF_hist_vars%bedout = .true.
DEF_hist_vars%sedinp = .true.
DEF_hist_vars%netflw = .true.
DEF_hist_vars%sedlayer = .true.
DEF_hist_vars%shearvel = .true.
DEF_hist_vars%seddiag = .true.
```

在当前 CoLM 源码中 `seddiag` 提供三粒级的：

- `f_es_raw_k`、`f_es_eff_k`：原始起悬与实际有效起悬；
- `f_d_raw_k`、`f_d_eff_k`：物理沉降与所有有效入床通量；
- `f_capdep_adv_k`、`f_capdep_exch_k`、`f_capdep_hill_k`、`f_capdep_other_k`：不同阶段的浓度上限导致的实际转床；
- `f_sedsto_k`：平均可见水体悬沙**固体体积**，`f_seddep_solid_k`：平均深层床沙**固体体积**；
- `f_ssc_precap_max`：history 窗口限幅之前的**最大总 SSC（mg/L）**；是窗口极值，不是日平均。

另外使用 `f_layer_k`（平均活动床层**含孔隙体积**）、`f_sedinp_k`、`f_netflw_k`、`f_shearvel`、`f_sedcon_k` 和 `f_sedout_k`。详细单位及算子关系见 [模型手册第 6 节](SEDIMENT_MODEL.md#6-诊断变量与-history-完整字典)。

**不能做的事**：将 `CapDep` 直接加到 `D_eff` 作为“总沉积”。在当前实现中 CapDep 是 `D_eff` 的一部分；`D_eff` 还可包含干载体入床和浅水外源直接入床，因此 `D_eff` 不等同 `D_raw`。

### 5.2 派生诊断量

| 指标 | 建议定义 | 物理用途 |
|---|---|---|
| `SSC_cap_ratio` | `sum(f_sedcon_k)/max_concentration` | 判断日均总 SSC 与体积分数上限的距离 |
| `CapDep_total` | `sum_k(capdep_adv + capdep_exch + capdep_hill + capdep_other)` | 限幅真实入床通量，单位 m³ solid/s |
| `CapDep_fraction` | `CapDep_total / sum_k(f_d_eff_k)`（分母>0） | 评价有效入床中限幅机制占比；不重复累计 |
| `Es_realization` | `sum_k(f_es_eff_k)/sum_k(f_es_raw_k)`（分母>0） | 起悬是否主要受库存限制 |
| `settling_fraction` | `sum_k(f_d_raw_k)/sum_k(f_d_eff_k)`（分母>0） | 物理沉降与有效入床的相对贡献，不预设 ≤1 |
| `susp_visible` | `sum_k(f_sedsto_k)` | 水体悬沙平均固体库存 |
| `bed_active_solid` | `(1−porosity)×sum_k(f_layer_k)` | 活动层平均固体库存 |
| `bed_deep_solid` | `sum_k(f_seddep_solid_k)` | 深层床沙平均固体库存 |
| `exchange_reconciliation` | `sum_k(f_es_eff_k−f_d_eff_k)` 对照 `sum_k(f_netflw_k)` | 检查有效交换的 history 口径 |

所有分式遇到零/缺测/无物理意义的分母应输出 NaN 而非 0。`f_ssc_precap_max` 是时间窗口最大值，`f_sedsto`、`f_layer` 是形态更新后时间加权平均；**严格质量守恒检查需用瞬时库存及完整边界积分通量**，不能把日均库存简单做逐日差分来判定守恒失败。

### 5.3 六类诊断图（每轮至少覆盖前四类）

1. **逐站 Q–SSC–SSL**：模拟全时段折线 + 观测配对日散点；标配对 n、均值比、KGE/PBIAS。对正 SSC/SSL 可用对数轴。
2. **Es/D 四分量 × 粒级**：clay、silt、sand 分面，分别绘制 `Es_raw`、`Es_eff`、`D_raw`、`D_eff`。
3. **CapDep 分阶段**：`adv/exch/hill/other` 按粒级与总量绘制；同步展示 `f_ssc_precap_max`、日均 SSC 和浓度上限线。
4. **悬沙—河床库存**：`f_sedsto`、`(1−porosity)f_layer`、`f_seddep_solid` 的长期趋势；可识别反复起悬/回床。
5. **驱动响应**：Q、`shearvel`、降水、`sedinp`、洪泛比例与 SSC/SSL，在丰枯季查看滞后。
6. **空间模式**：上限接近比例、`Es−D`、`CapDep/D_eff`、坡面输入空间分布，叠加已验证的观测河网匹配位置。

同时报告两套统计：**完整模拟时段**和**同一变量的观测配对日**；不得把平台比例、平均浓度或过程平均通量跨口径混用。

## 6. 诊断决策树与下一轮调参

| 已观察到的表现 | 优先确认的诊断 | 如果成立，下一步可以测试的假说 |
|---|---|---|
| Q 严重偏大；SSC/SSL 同时偏大 | `f_discharge`、水深、储水、`shearvel`、汇水面积、时间对齐 | [H] 水动力或站点匹配是主因；先不调泥沙 |
| Q 尚可；`sedinp` 很大，`capdep_hill` 显著 | 降雨、`sed_slope` 的值/单位、`yield_coefficient`、坡度指数 | [H] 坡面供沙过高；优先做关闭源项与系数扰动 |
| `Es_raw` 大、`Es_eff` 高；`capdep_exch` 很大 | 粒径、剪切临界值、`shearvel`、活动床库存 | [H] 河床起悬过强；先核对床层可用量与临界剪切 |
| `D_raw` 弱且细粒 SSC 持续高 | `setvel`、`settling_multiplier`、clay/silt/sand 配置 | [H] 沉降能力不足或粒级方案偏细 |
| `capdep_adv` 很大 | 河道水量端点、正负向流量、干湿切换和 donor | [H] 水体载体/平流导致强限幅回床 |
| `capdep_other` 很大 | BIF、levee、protected-water 交换 | [H] 分汊/堤防路径贡献不可忽略 |
| SSC 频繁平台，但床层库存不减 | Es/D、`CapDep`、`sedsto`、`layer`、`seddep` | [H] “起悬→限幅回床→再起悬”循环 |
| SSC 尚可、SSL 高估 | Q、`f_sedout` 单位/时间、流向、配对日期 | [H] 输运/流量而非局部浓度过程主导 |
| 修改 `max_timestep_s` 后结果剧变 | morph steps、数值收敛、守恒误差 | [H] 数值时间离散敏感，不能当成物理率定成功 |

### 6.1 建议的第一批受控实验（**尚未执行**）

| 实验 ID | 唯一计划改变 | 主要目的和评估条件 |
|---|---|---|
| `BASE_2001_2002` | 所有参数不改，先确认实际输入 | 重现 baseline，锁定各站各变量配对日 |
| `P1-01_yield_zero` | `yield_coefficient = 0` | **归因试验，不是最终率定参数**。隔离当期降水坡面源；库存仍可继续输出，所以 SSC 不一定立即为零 |
| `P1-02_yield_0p1x` | `yield_coefficient × 0.1` | 观察产沙通量变化与 SSC/SSL 响应比例、床层滞后 |
| `P1-03_settling_2x` | `settling_multiplier × 2` | 检查沉速、`D_raw`、库存与日均 SSC 的敏感性，不预设该值物理合理 |
| `P1-04_dt_half` | `max_timestep_s × 0.5` | **数值收敛/算子时序诊断**，非物理优化 |
| `P1-05_match_time_check` | 不更改 CoLM 输出，只重新审核 OpenBench 单位/时间/匹配 | 分离评估链造成的指标变化与真实模型改变 |

**不要第一轮直接修改 `max_concentration`。** 目前 0.01 m³/m³ 对应默认颗粒密度下约 26,500 mg/L，属于质量守恒式回床限幅；压低上限可能掩盖供沙/起悬/沉降过强。可在明确科学问题的独立敏感性试验中研究它，但不能用来冒充机制改进。

**不可把 `cfl_adv` 当作当前有效率定参数**：在固定的 CoLM commit 中它为兼容字段，悬沙平流已使用 donor-inventory 限制；`max_timestep_s` 属于数值步进，`unit_conversion` 属于公式单位约定。真正可解释的物理参数优先考虑产沙系数/坡度/雨强指数、粒径/沉降系数及床层结构，但要结合第 6 节证据再选择。

### 6.2 如何选下一轮的“更优参数”

先通过结构/守恒、站点、单位和样本门禁，再依序比较：

1. **Q**：若 Q 极差，优先修水文/匹配，不让 SSC、SSL 参数补偿 Q。
2. **绝对量级**：减少 SSC/SSL 的配对均值比、PBIAS 和极端异常日，不只看综合分。
3. **物理预算**：改善有可追溯的 Es/D、CapDep、库存变化；不允许仅靠数值上限造成“漂亮”的低浓度。
4. **时序表现**：同站 KGE/NSE、Spearman、洪峰时滞、丰枯季一致性。
5. **稳健性**：最差站点不应被牺牲；记录逐站全表和中位数；最后查看独立留出年份/站点。

对于 SSC/SSL 有千倍差异的情况，可先用 `log10(mean_ratio)` 表示数量级，再检查 PBIAS、KGE 和样本数。任何指标阈值若用于参数筛选，应在实验前写在 manifest 中，而不是根据结果事后修改。

### 6.3 Restart 与 spin-up 约束

CoLM 当前 sediment restart 校验许多参数及静态属性，包括粒径、沉速、床层、`sed_frc`/`sed_slope` 等。不同实验改变这些物理配置后，旧 restart 可能不兼容，且即使技术上能载入也未必物理等价。

- 为比较参数效果，保持相同 forcing、路由网格、参考数据、匹配和评估期。
- 发生参数更改时使用**一致的初始化和自旋方案**，并写明是 cold start 还是兼容 restart；不要手工篡改 restart 元数据来绕过校验。
- 如果实验只改 `max_timestep_s`，仍需检查历史结果是否数值收敛。
- 用一致的预热期平衡床层库存，报告关键 `sedsto`、`layer`、`seddep` 对初态的敏感性。

## 7. 每轮实验的文件与机器可追溯记录

### 7.1 建议目录结构

```text
CoLM202X_SYSU/
└── docs/
    └── sediment/
        ├── SEDIMENT_MODEL.md
        ├── CALIBRATION_VALIDATION.md
        └── experiments/
            ├── BASE_2001_2002.md
            ├── P1-01_yield_zero.md
            └── ...

# 大体积 NetCDF、日志、OpenBench 输出留在 HPC 或归档存储；
# Git 保留路径、SHA、评估结果摘要、绘图和可复算脚本。
```

### 7.2 新实验 manifest 模板

```yaml
experiment_id: P1-01_yield_zero
status: planned  # planned / running / complete / failed / rejected / accepted
parent_experiment: BASE_2001_2002
date: YYYY-MM-DD
hypothesis: "An excessive hillslope source dominates SSC high bias"
expected_change: "sedinp becomes zero while preexisting bed stores may still export"
falsification: "High SSC persists with small incoming source and bed exchange dominates"

code:
  colm_repo: jellywzx/CoLM202X_SYSU
  colm_sha: "<exact SHA>"
  colm_dirty_diff: "<file or clean>"
  executable_sha256: "<sha256>"
  openbench_repo: zhongwangwei/OpenBench
  openbench_sha: "<exact SHA>"
  openbench_runtime_path: "<python: openbench.__file__>"

inputs:
  main_namelist: "<path + sha256>"
  sediment_namelist: "<actual resolved path + sha256>"
  history_namelist: "<path + sha256>"
  forcing: "<name/version/path/date/sha256>"
  unitcatchment_netcdf: "<path/hash>"
  initial_conditions: "<cold or validated restart path/hash>"
  spinup_protocol: "<years/cycles>"
  simulation_period: "2001-01-01/2002-12-31"
  target_grid_res_deg: 0.25
  parameter_changes:
    yield_coefficient: {from: "<baseline>", to: 0}

reference:
  dataset: OpenBench_Sediment_Daily
  version: "<version/DOI/hash>"
  qc_rule: "0-1 or documented alternative"
  station_map: "<path + sha256>"
  time_alignment: "<actual settings and native offset>"
  unit_conversions: "<verified steps>"

outputs:
  colm_log: "<path>"
  history_manifest: "<files/hash>"
  restart: "<path/hash>"
  openbench_config: "<path + sha256>"
  openbench_results: "<directory>"
  paired_daily_csv: "<path>"
  station_metrics_csv: "<path>"
  sediment_process_plots: "<directory>"
  process_budget: "<path>"

result:
  exit_code: "<value>"
  conservation_pass: "<true/false>"
  paired_n_by_station_variable: "<table/path>"
  full_period_diagnostics: "<summary>"
  paired_period_diagnostics: "<summary>"
  process_explanation: "<evidence>"
  hypothesis_result: inconclusive  # supported/refuted/inconclusive
  decision: pending                 # accept/reject/retest
  next_experiment: "<id + exact proposed change>"
```

**新实验不能覆盖旧实验目录**。如果 OpenBench 更新了比较图代码，需保留原始 run、评估数值和配置，明确“仅重新绘图”还是“重新计算指标”。

### 7.3 标准的站点输出表

```text
experiment_id,station_uid,station_name,variable,tim_res,
n_sim,n_obs_valid,n_paired,n_joint,
mean_sim_paired,mean_obs_paired,mean_ratio_paired,
pbias_percent,kge,nse,pearson_r,spearman_r,
cap_near_count_full,cap_near_count_paired,
model_grid_id,matching_method,matching_quality,
time_offset_applied,qc_rule
```

缺测/未定义的指标填 `NaN`，不能填 0。另保存 `paired_daily.csv`（时间、各变量 sim/ref、QC），方便从原始配对值重算指标而不受后续 OpenBench 版本变化影响。

## 8. 每次向 ChatGPT / Codex 请求下一轮调参时，提供什么

推荐按以下顺序提供：

1. 当前实验 ID、父实验、CoLM/OpenBench commit、唯一改动的参数和是否保持相同 restart/spin-up。
2. OpenBench 逐站逐变量的 `n`、模拟/观测配对日均值、KGE/NSE/PBIAS/相关和历史 baseline 的同口径差值。
3. 模型日志中质量守恒、异常浓度、数值稳定性、极值位置、`SED_PERF` 或 morph step 状况。
4. 站点或匹配河段的 `f_sedinp`、`Es_raw/eff`、`D_raw/eff`、四类 CapDep、`f_ssc_precap_max`、`f_sedsto`、`f_layer`、`f_seddep_solid` 图/统计。
5. 如观测、站点、配对日期或 time_alignment 有变化，先用**不变的 CoLM 原输出**重评旧配置，隔离数据处理误差。

期望回答固定包含：**①代码/数据事实，②最可能的物理解释，③尚不能排除的机制，④下一轮参数及预期响应，⑤反证条件与验收指标**。当信息不足时优先索取针对性诊断，不凭“SSC 高了 100 倍”直接要求把某个系数缩小 100 倍。

## 9. 未完成验证项目（持续维护）

- [ ] 找到真实 2001–2002 两年 CoLM 运行的 main/sediment/history namelist、forcing、routing 文件版本及 restart 来源。
- [ ] 重读四站 `SED006136/139/140/141` 的 OpenBench 原始配对文件，核对各 `n`、倍数、`Overall_Score` 评分规则和平台计算口径。
- [ ] 站名/河段/上游面积匹配逐站复核，尤其 **Porto Velho、Fazenda Vista Alegre** 与硬编码 **Serrinha/Manacapuru/Obidos** 的名单差异。
- [ ] 核实模型运行实际 `grain_density`、`grain_diameter`、`max_concentration`、`active_layer_depth`、`ignore_depth_m`。
- [ ] 检查 `sed_slope` 原始单位、极值分布和洪泛坡面暴露面积；核对降雨非线性矩及有效缺测掩膜。
- [ ] 核验 OpenBench 3.0.8 对 `f_discharge` 与 `f_sedcon/f_sedout` 的实际时间偏移及单位转换。
- [ ] 用当前 `seddiag` 输出开展完整时段/观测配对日双口径六类过程分析，并记录模块物理守恒和数值收敛。
- [ ] 固化 `experiments/BASE_2001_2002.md`，再开始 P1-01/02/03/04。
- [ ] 增加真正独立的年份/站点验证；没有留出数据时不宣称泛化技能。

## 10. 源码参考（固定于审计 SHA）

- [模型原理与完整参数/变量](SEDIMENT_MODEL.md)
- [CoLM 主沉积模块](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90)
- [CoLM 亚马孙示例配置](../../run/Amazon_0p25.nml)、[实验泥沙参数](../../run/sediment_parameter.nml)、[history 选择](../../run/history_sediment.nml)
- [OpenBench reference_catalog.yaml](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/data/registry/reference_catalog.yaml#L4296-L4388)
- [OpenBench CoLM model_catalog.yaml](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/data/registry/model_catalog.yaml#L1158-L1419)
- [OpenBench station_matcher.py](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/data/station_matcher.py)
- [OpenBench 时间偏移实现](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/data/_processing_time_adjustments.py)
- [OpenBench 指标源码](https://github.com/zhongwangwei/OpenBench/blob/f2d841609d3c214a8726e1abde5f4ca177cefc60/src/openbench/core/metrics.py)
- [KGE benchmark：Knoben et al. (2019)](https://doi.org/10.5194/hess-23-4323-2019)

*以后每轮实验新增独立记录，更新本文的结论和未完成项目，但不覆盖过去的观测、运行与参数决策历史。*
