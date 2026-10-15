# CoLM202X-SEDIMENT：泥沙模块原理、参数与诊断手册

> **文档定位**：面向 CoLM202X_SYSU 中的河道泥沙模块开发、调参、运行诊断和验证。不是最终发表的方法描述，也不代替逐次实验记录。  
> **代码基线**：`jellywzx/CoLM202X_SYSU`，`master` @ [`32ea6ce72b0a8ced75f5859d47c3d6cb3632cd07`](https://github.com/jellywzx/CoLM202X_SYSU/commit/32ea6ce72b0a8ced75f5859d47c3d6cb3632cd07)。核对日期：2026-10-15。  
> **证据标记**：**[S]** 已核对上述源码/仓库配置；**[O]** 项目此前实验或对话报告（不保证对应此 commit）；**[H]** 物理机制假设或调参建议；**[T]** 尚需实际数据、运行或文献检验。  
> **维护原则**：每次修改模块物理、参数读取、history 单位或 restart 结构时，更新本文件的版本基线和有关条目；不要把历史结论当作最新版源码事实。

## 1. 快速索引与适用范围

| 内容 | 核心文件 | 用途 |
|---|---|---|
| 具体泥沙物理与数值算子 | [`main/TRACER/MOD_Tracer_Particle_Sediment.F90`](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90) | 参数读取、产沙、推移质、悬沙、沉降、床层、history、restart |
| tracer 描述符/生命周期 | [`main/TRACER/README.md`](../../main/TRACER/README.md)、[`MOD_Tracer_Lifecycle.F90`](../../main/TRACER/MOD_Tracer_Lifecycle.F90) | SEDIMENT 作为 `particle` provider 的注册和路由回调 |
| 河湖水动力耦合入口 | [`main/HYDRO/MOD_Grid_RiverLakeFlow.F90`](../../main/HYDRO/MOD_Grid_RiverLakeFlow.F90) | 供给逐水动力子步的水深、水量、流速、流量和降水映射 |
| 实验用泥沙参数 | [`run/sediment_parameter.nml`](../../run/sediment_parameter.nml) | 当前仓库亚马孙示例所指向的独立调参文件 |
| 另一套标准参数 | [`run/standard_sediment_parameter.nml`](../../run/standard_sediment_parameter.nml) | **与上一个文件不同，不能视为同一套参数** |
| 亚马孙主配置示例 | [`run/Amazon_0p25.nml`](../../run/Amazon_0p25.nml) | tracer、routing、输出频率和参数文件映射 |
| 泥沙 history 选择 | [`run/history_sediment.nml`](../../run/history_sediment.nml) | 基础泥沙变量及 `seddiag` 开关 |
| 历史实现说明 | [`main/TRACER/README.md#sediment-particle-provider`](../../main/TRACER/README.md#sediment-particle-provider) | 降水、restart、过程诊断的补充约定 |

**[S] 模块适用条件**：源码预处理条件是 `TRACER` **且** `GridRiverLakeFlow`；仅编译 tracer 不会自动启用 SEDIMENT。SEDIMENT 被注册为 `PARTICLE / PROVIDER / NONE`，专有状态由泥沙 provider 维护，无单独顶层 `SEDIMENT_ON` 开关。[源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L1-L40)

**[S] 基本范围**：该模块对 CoLM 的 **GridRiverLakeFlow 河网/单元汇水区**进行多粒级悬浮泥沙、推移质、河床—水体交换和降水驱动的坡面泥沙输入。历史站点序列验证还需要独立的站点—河网匹配；不能把模型控制台打印的“最近格网”视为经过水文一致性检验的站点匹配。

## 2. 代码结构与执行顺序

### 2.1 生命周期 hook

`register_sediment_tracer_provider` 将以下函数绑定到 tracer 生命周期（[源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L278-L301)）：

| Hook | 实现函数 | 作用 |
|---|---|---|
| `route_init` | `grid_sediment_init` | 读参数和河网静态数据、分配泥沙状态 |
| `route_read_restart` | `read_sediment_restart` | 读取并校验沉积物 restart |
| `route_forcing_put` | `sediment_forcing_put` | 按有效面积与时间积累降水及非线性产沙矩 |
| `route_diag_accumulate` | `sediment_diag_accumulate` | 从水动力子步积累流速、水深、水量、正负向流量等 |
| `route_calc` | `grid_sediment_calc` | 每个 routing period 执行泥沙过程 |
| `route_sediment_bif` | `sediment_bif_accumulate` | 分汊通量交接（启用时） |
| `route_sediment_levee` | `sediment_levee_repartition` | 堤防可见/保护水体的泥沙交接（启用时） |
| `route_history` | `write_sediment_history` | 通过 unit-catchment dispatcher 写 history |
| `route_flush_history` | `flush_sediment_history` | 冲刷/清零统计窗口累积量 |
| `route_write_restart` | `write_sediment_restart` | 写状态、参数元数据和诊断积累量 |
| `route_final` | `grid_sediment_final` | 释放状态 |

### 2.2 逐步流程

```text
主 namelist: DEF_TRACER_NAMES='SEDIMENT', TYPES='particle'
  └─ tracer_lifecycle_init / register_sediment_tracer_provider
      └─ grid_riverlake_flow_init → grid_sediment_init
          ├─ 读 unitcatchment NetCDF: sed_n、slope_layers、sed_frc、sed_slope
          ├─ 读 DEF_TRACER_PARAM_FILES 指向的泥沙 namelist
          └─ 初始化 grain/settling velocity/bed inventory
水文/陆面步
  ├─ land precipitation → unit catchment → sediment_forcing_put
  └─ GridRiverLakeFlow 的 HYDRO 子步
       └─ sediment_diag_accumulate（水深、速度²、储水量端点、正向/反向流）
Routing period 结束 → grid_sediment_calc
  ├─ 形成本期平均水动力量、摩阻速度与临界剪切
  ├─ 每个 dt_morph = min(剩余 routing 时长, max_timestep_s)
  │    ├─ 悬沙/推移质平流（取 donor 快照，限制供给库存）
  │    ├─ 坡面输入量计算（降水阈值、坡度、暴露面积）
  │    ├─ 河床起悬与沉降 Es/D（并实施浓度限幅）
  │    ├─ 坡面泥沙入河（浅/干河段直接入床）
  │    ├─ 活动层/沉积层重分配
  │    └─ 过程 history 的 dt 加权积分与窗口最大值
  └─ 更新 sedsto / sedcon，后续写 history、restart
```

**[S] 重要数值变化**：当前版本每个 `dt_morph` 进行 **一次**基于 donor 库存的平流计算。`cfl_adv` 仅为历史兼容而保留，并**不再控制平流子步长**；真正的形态子时段长度受 `max_timestep_s` 限制。水动力模块仍有独立时间步/稳定性约束。不能把旧日志中的“sediment CFL 控制大量 substeps”直接套用到此版本。[源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L62-L84)、[主计算](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L1050-L1057)

## 3. 状态定义、单位与质量守恒

下表是理解 history、限幅及守恒的基础。`k=1..nsed` 表示粒级，`i` 表示 unit catchment。

| 内部变量 | 内部单位 | 物理含义与注意事项 |
|---|---|---|
| `sedsto(k,i)` | m³ **固体** | 可见水体中悬浮泥沙的**预报主状态**，不随水体积改变而凭空增减 |
| `sedcon(k,i)` | m³ solid / m³ water | 从 `sedsto` / HYDRO 水体积推导的悬沙体积分数，不是独立的质量状态 |
| `layer(k,i)` | m³ **床层总体积** | 活动床层中各粒级；与 m³ 固体比较须乘 `(1-porosity)` |
| `seddep(k,layer,i)` | m³ **床层总体积** | 活动层以下的命名沉积层 |
| `sedsto_protected(k,i)` | m³ **固体** | 受保护水体的悬沙库存，与 `f_sedsto_k` 的“可见”库存不同 |
| `sedbed_protected(k,i)` | m³ **固体** | 保护区内不移动的床沙库存 |
| `sedout(k,i)` | m³ solid/s | 悬浮泥沙的有向出流；可有负号（逆流） |
| `bedout(k,i)` | m³ solid/s | 推移质有向通量；**不要直接加入悬浮输沙量 SSL** |
| `sedinp(k,i)` | m³ solid/s | 外部坡面泥沙输入 |
| `netflw(k,i)` | m³ solid/s | 水体→床层为负；床层→水体为正；包含浅水输入直接入床的负贡献 |
| `exch_es_raw` / `exch_d_raw` | m³ solid/s | 限制前起悬、可用悬沙库存限制后的物理沉降 |
| `exch_es_eff` / `exch_d_eff` | m³ solid/s | 真正执行后的**有效**起悬、总有效沉积 |
| `shearvel` / `critshearvel` | m/s | 水动力摩阻速度与粒级临界摩阻速度 |
| `susvel` / `setvel` | m/s | 起悬速度参数化结果与粒级沉速 |

**[S] 核心约定**：`sedsto` 是悬沙固体体积主状态，`sedcon = sedsto / visible_water_volume`；降水或水动力改变了水体积时可以改变 SSC，但不应直接改变固体质量。[源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L2763-L2814)

对于局部床层的固体体积，可按 `(1-λ) × (sum(layer)+sum(seddep))` 计算，其中 `λ=porosity`。完整河网收支还必须加入可见/保护悬沙、保护床库存、外部坡面输入、边界和分汊/堤防通量；不能仅用 `f_sedsto` 与 `f_layer` 代替整个质量守恒账本。

## 4. 物理过程及本版本的实际计算

### 4.1 降水驱动的坡面泥沙源项 [S]

原理来自代码注释所述 Sunada & Hasegawa (1993) 型经验产沙方案；具体实现以源码为准：

1. 陆面每一 forcing step 把降水映射到 unit catchment。降水是**强度量**，按“本步有效 patch 面积”加权平均；缺失 patch 不当作零降水，真实干 patch 仍计入有效面积。
2. 阈值为编译常数 `SED_PRECIP_THRESHOLD_MM_DAY = 10 mm/day`；只有本步 `P×86400 > 10` 时，非线性产沙矩 `(P_mm/h)^precipitation_exponent` 才计入分子；低于阈值的有效步长依然计入分母。
3. `sed_precip_yield / sed_precip_time` 形成逐步产沙矩的加权平均，**不是**先求平均降水再取幂（避免 Jensen 偏差）。
4. `calc_sediment_yield` 按 `yield_coefficient × 雨强矩 × sed_slope^slope_exponent × 面积 × 暴露层份额 × unit_conversion × sed_frc / 3600` 累加到 `sedinp(k)`；洪泛比例 `fldfrc` 调节可贡献的坡面层。

主要源码：[`sediment_forcing_put`](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L2375-L2423)、[`calc_sediment_yield`](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L3605-L3645)。

**[T] 必须核查**：`sed_slope` 来自 routing NetCDF，当前代码未在 namelist 中给出统一的“度/百分数/无量纲坡度”换算说明。修改 `slope_exponent` 前，应查看输入 NetCDF 的属性、分布与极端值，不能只凭日志里的 `max_slope` 猜测单位。此方案按 unit catchment 均一暴露面积计算，并非空间显式的坡面侵蚀—输送连通性模型。

### 4.2 水动力摩阻速度与起动条件 [S]

源码使用时间平均的速度平方 `<v²>` 与 Manning 粗糙度、水深估计摩阻速度：

```text
u* = sqrt( g × n_Manning² × <v²> × h^(-1/3) )
```

临界条件先由 `calc_critical_shear_vel_sq(diam)` 的分段经验式计算 **(cm/s)²**，再通过粒径混合的 Egiazaroff 修正得到 `critshearvel(k)`（取平方根并乘 0.01 转 m/s）；不要把分段式的返回值误当作 m/s。`calc_suspend_velocity` 在 `u* <= u*crit` 时令起悬速度为零，并使用沉速、`von_karman` 和 `porosity` 形成粒级起悬速度。[源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L2676-L2760)

### 4.3 沉速、起悬与沉积 [S]

逐粒级沉速（使用密度比值，`ν` 为运动黏度，`d_k` 为粒径）：

```text
w_s,k = settling_multiplier × [
  sqrt( (2/3) × ((rho_s-rho_w)/rho_w) × g × d_k + (6ν/d_k)² )
  - 6ν/d_k
]
```

`calc_sediment_exchange` 对有效水体积与 `bed_area × ignore_depth_m` 判断是否计算床水交换。原始起悬项近似为：

```text
Es_raw,k = susvel_k × (1 - porosity) × bed_area × layer_k / sum(layer)
```

沉积由沉速、底面积、悬沙浓度及近床剖面因子控制；低剪切接近静水时仍可沉降。源码把物理沉积量限制在现有悬沙 `sedsto_k/dt` 内，并对活动床库存限制实际可起悬量：

```text
D_raw,k = min(w_s,k × bed_area × sedcon_k × profile_factor,
              sedsto_k / dt)
netflw_k (初步) = Es_raw,k - D_raw,k
```

有效起悬和沉积 `Es_eff/D_eff` 会进一步反映库存限制、限幅沉积及浅水输入，**不能把 `D_eff` 当作纯粹的自然沉降**。[源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L3259-L3388)

### 4.4 悬沙平流与推移质 [S]

- `calc_sediment_advection` 在一个 `dt_morph` 内分别计算正向和反向的 donor-side 通量，并通过共同库存约束限制同一供给格网向多个出口输出的泥沙，减少负库存与“一个时段重复输出新输入泥沙”的风险。
- 悬沙通量基本形式为 `Q_signed × donor_concentration`，负流选择反向 donor；不能用净流量为零推断没有正反向总输运。
- `bedout` 采用源码标注的 Ashida–Michiue 型摩阻速度关系，含固定系数 `SED_BEDLOAD_COEFF=17`，且按正反向水流占比拆分贡献。`bedout` 与 `sedout` 在代码中分别记账。
- `DEF_USE_BIFURCATION`、`DEF_USE_LEVEE` 启用时分别有分汊/受保护水体通量与库存协调路径；是否足以表达自然洪泛平原泥沙滞留需另行研究。

[源码：平流](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L2817-L2965)、[正反向通量和推移质](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L2968-L3095)。

### 4.5 浓度上限不是泥沙质量损失 [S]

`sed_max_conc` 是**各粒级体积分数之和**的上限，单位 `m³ solid/m³ water`。超过允许体积 `water_volume × sed_max_conc` 的部分按粒级比例转移至床层（`layer`），同步修正 `netflw`、`D_eff` 等诊断，而非丢弃质量。主要阶段：

| 限幅阶段 | 主要入口 | history 诊断 |
|---|---|---|
| 正/反向平流前后（有水载体） | `calc_sediment_advection` | `f_capdep_adv_k` |
| 床水起悬—沉积交换后 | `calc_sediment_exchange` | `f_capdep_exch_k` |
| 湿河段坡面输入后 | `apply_sediment_input` | `f_capdep_hill_k` |
| 分汊、可见/保护堤防等其他路径 | 相关 BIF/levee 函数 | `f_capdep_other_k` |

浅/干单元接收到 `sedinp` 时可能**直接进入床层**，计入 `D_eff`、`netflw`，但**不是** CapDep；水体彻底变干后的库存回床也要单独理解。[源码：坡面入河](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L3206-L3256)、[床水交换](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L3363-L3385)。

当 `grain_density=2650 kg/m³` 时，`sed_max_conc=0.01` 对应 **26,500 mg/L** 的总 SSC 上限。若很多站点出现接近上限的 SSC，优先诊断供沙、起悬、沉降与反复“起悬→超限回床→再起悬”的可能性，而不是通过继续压低上限来掩盖物理偏差 **[H]**。

### 4.6 床层重分配与局部守恒 [S]

`calc_layer_redistribution` 按 `active_layer_depth × bed_area` 维持活动层厚度，多余物质下移至 `ndeposit_layers` 个床层；当活动层不足时从下方补充。`layer` 和 `seddep` 存的是**含孔隙床层总体积**，检查守恒时使用 `(1-porosity)` 转换为固体体积。[源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L3391-L3501)

局部守恒检查为：

```text
residual  = after_solid_volume - before_solid_volume - expected_change
tolerance = SED_BALANCE_ABS_TOL + SED_BALANCE_REL_TOL
            × max(abs(before),abs(after),abs(expected_change))
```

默认绝对和相对容差分别为 `1e-10`；一旦异常，会打印 worker、河网单元、粒级的储量与通量并终止，不能直接关闭质量守恒检查以“通过”实验。[源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L3505-L3601)

## 5. 完整可配置参数字典（当前源码 19 项）

### 5.1 文件选择与两套参数的区别

主配置通过 `DEF_TRACER_PARAM_FILES` 指定 **实际读取的绝对路径/相对路径**。建议复制实验参数文件，而不是每次覆盖公共 `standard_sediment_parameter.nml`：

```fortran
DEF_TRACER_NUM         = 1
DEF_TRACER_NAMES       = 'SEDIMENT'
DEF_TRACER_TYPES       = 'particle'
DEF_TRACER_PARAM_FILES = 'SEDIMENT:/absolute/path/to/my_experiment_sediment_parameter.nml'
```

参数文件内，`&nl_colm_sediment_parameter` 的 `DEF_SEDIMENT%...` 对应下表中的 19 项。表中“实验文件”为 `run/sediment_parameter.nml`，“标准文件”为 `run/standard_sediment_parameter.nml`；“源码回退值”是**参数文件未覆盖时**模块内部使用的默认值，不代表本次运行一定使用该值。[参数读取源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L407-L751)

| `DEF_SEDIMENT%` 字段 | 实验文件 | 标准文件 | 未覆盖时源码值 | 物理作用、单位及约束 |
|---|---:|---:|---:|---|
| `nsed` | 3 | 3 | 以 NetCDF `sed_n` 为准 | 粒级数；如显式设置必须与 routing NetCDF `sed_n` 相同 |
| `grain_diameter` | 2e-6, 2e-5, 1e-4 | 2e-5, 2e-4, 2e-3 | 2e-4, 2e-3, 2e-2 | 各粒级粒径 [m]；必须为正，最多支持 100 类 |
| `grain_density` | 2650 | 2650 | 2650 | 颗粒密度 [kg/m³]；程序内部换成 g/cm³；必须大于水密度 |
| `water_density` | 1000 | 1000 | 1000 | 水密度 [kg/m³]；内部换成 g/cm³ |
| `porosity` | 0.4 | 0.4 | 0.4 | 床层孔隙率，`0 ≤ porosity < 1` |
| `ndeposit_layers` | 5 | 5 | 5 | 沉积层数，正整数；床层分层结构 |
| `ignore_depth_m` | **0.05** | **0.01** | 0.05 | 浅水床水交换/坡面泥沙直入床阈值 [m]，`≥0` |
| `active_layer_depth` | **0.005** | **0.05** | 0.05 | 活动层厚度 [m]，`>0`；小于最大粒径时仅警告，不自动修改 |
| `viscosity` | 1e-6 | 1e-6 | 1e-6 | 运动黏度 [m²/s]，决定沉速 |
| `von_karman` | 0.4 | 0.4 | 0.4 | von Kármán 常数，影响起悬/近床剖面 |
| `settling_multiplier` | 1 | 1 | 1 | 无量纲沉速倍率，`>0` |
| `yield_coefficient` | 0.01 | 0.01 | 0.01 | 坡面产沙经验系数，`≥0`；单位依模型雨强/坡度约定 |
| `slope_exponent` | 2 | 2 | 2 | 产沙对 `sed_slope` 的非线性敏感性，`≥0` |
| `precipitation_exponent` | 2 | 2 | 2 | 产沙对逐步有效雨强的非线性敏感性，`≥0` |
| `unit_conversion` | 1e-6 | 1e-6 | 1e-6 | 坡面产沙公式的换算因子，`≥0`，须与原始公式的单位配套 |
| `cfl_adv` | 0.5 | 0.5 | 0.5 | **目前不参与平流子步控制**；仍需满足 `0<cfl_adv≤1` 并写入 restart 元数据 |
| `max_timestep_s` | **1800** | **1800** | 3600 | 泥沙形态更新最大时段 [s]，`>0` |
| `bed_depth` | 未指定 | 未指定 | **10** | 初始命名沉积床层总深度 [m]，须覆盖活动层+各命名层 |
| `max_concentration` | 未指定 | 未指定 | **0.01** | 悬浮固体总体积分数限幅 [m³/m³]，`0<值≤1` |

**[S] 高优先级提醒**：

- 当前仓库中的“实验/标准”粒径差异达到 **10–20 倍**，`ignore_depth_m` 为 **0.05 vs 0.01 m**，`active_layer_depth` 为 **0.005 vs 0.05 m**。对比实验之前必须明确**真正载入的文件**。
- 若要显式配置上限和床层深度，在 `&nl_colm_sediment_parameter` 中添加，例如：

```fortran
DEF_SEDIMENT%max_concentration = 0.01
DEF_SEDIMENT%bed_depth         = 10.0
```

- `MAX_SED_CONC`/`SED_DEFAULT_MAX_CONC` 是编译期**默认常量**，实际各条限幅路径使用 `sed_max_conc`，其值可以被 `max_concentration` 覆盖；**不要误认为所有路径永远固定 0.01**。
- `sed_frc`（粒级比例）和 `sed_slope`（坡面层坡度）来自 `DEF_UnitCatchment_file` 指向的 NetCDF，**不在上述 19 个标量 namelist 参数中**；`sed_frc` 每格正规化为和 1，零和时回退到各粒级均等比例。[静态数据源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L2500-L2552)
- 生效值以启动日志的 `SEDIMENT EFFECTIVE PARAMETERS`、粒级直径和 `settling_velocity(k)` 为准，而不是仅靠编辑器里打开的 nml。[启动输出源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L5209-L5272)

### 5.2 其他 tracer 元数据（不是 19 个泥沙物理参数）

`&nl_colm_tracer_parameter` 还包含 `DEF_TRACER%unit_kind='volume_fraction'`、`mol_weight`、`ref_ratio`、`init_delta`、`reactive_decay_rate`。其中 `init_delta` 对粒子模型是预留通用字段；`mol_weight/ref_ratio` 不应作为泥沙物理调参手柄。`&nl_colm_tracer_forcing` 当前 `forcing_num=0`，外部 `erosion_yield` forcing 角色**尚未实现**，坡面源项来自模型内部降水，而不是该预留块。[参数文件](../../run/sediment_parameter.nml)

### 5.3 不能直接在 sediment namelist 中设置的编译常数 [S]

| 源码常数 | 当前值 | 影响及提醒 |
|---|---:|---|
| `SED_BEDLOAD_COEFF` | 17 | 推移质关系的固定系数 |
| `SED_PRECIP_THRESHOLD_MM_DAY` | 10 mm/day | 强降水产沙开启阈值；调雨强指数时须考虑这一门槛 |
| `EXCH_SHEARVEL_MIN` | 1e-4 m/s | 静水沉降/剖面修正转变阈值 |
| `EXCH_SHEARVEL_BLEND` | 2e-4 m/s | 沉降剖面平滑过渡上界 |
| `EXCH_ZD_MAX` | 100 | 剖面因子的数值保护 |
| `SED_BALANCE_ABS_TOL` | 1e-10 | 局部固体体积守恒绝对误差容差 |
| `SED_BALANCE_REL_TOL` | 1e-10 | 局部固体体积守恒相对容差 |
| `SED_NEAR_DRY_DEPTH` | 1e-4 m | **当前泥沙源文件只定义并打印**；未发现它被用作此模块实际平流门槛，不应当成生效的可调参数 |
| `SED_ADV_DIAG_PERIODIC` | 5 | 初期平流统计的打印次数门槛，不是物理参数 |

[S] 数值代码明确使用的是 `ignore_depth_m` 控制交换/坡面入床浅水门槛；HYDRO 可能有独立近干判断，不可与仅打印的 `SED_NEAR_DRY_DEPTH` 混为一谈。

### 5.4 基础水动力与实验配置（不同于泥沙参数）

`run/Amazon_0p25.nml` 示例：陆面时间步 `DEF_simulation_time%timestep=1800 s`，`DEF_GRIDBASED_ROUTING_MAX_DT=3600 s`，泥沙参数文件中 `max_timestep_s=1800 s`，三者控制**不同的时间尺度**。该例设置日尺度 history、每月存放一个文件。`DEF_HIST_vars_namelist` 和 `DEF_TRACER_PARAM_FILES` 在当前仓库示例中使用服务器**绝对路径**，甚至可能指向不同的旧 checkout；迁移运行目录时必须重设，不能假定仓库所见文件就是运行时实际读取的文件。

[T] 项目既有报告显示“全球 unitcatchment 文件 + 区域 namelist”可能导致 unitcatchment 层仍有区域外活动格网。**history 仅为区域**不意味着 routing 计算也仅限区域；需用实际 `numucat`、单元经纬度、运行时间和空间域校验，不把这一历史观察写成所有版本必然行为。

## 6. 诊断变量与 history 完整字典

### 6.1 启用方式 [S]

在主 namelist 中选择 history 变量配置文件；若采用关闭所有默认 history 字段的方式，显式选择泥沙诊断：

```fortran
! Main namelist
DEF_HIST_vars_out_default = .false.
DEF_HIST_vars_namelist    = '/absolute/path/to/history_sediment.nml'
DEF_HIST_FREQ             = 'DAILY'
DEF_HIST_groupby          = 'MONTH'

! In &nl_colm_history, e.g. run/history_sediment.nml:
DEF_hist_vars%sedcon   = .true.
DEF_hist_vars%sedout   = .true.
DEF_hist_vars%bedout   = .true.
DEF_hist_vars%sedinp   = .true.
DEF_hist_vars%netflw   = .true.
DEF_hist_vars%sedlayer = .true.
DEF_hist_vars%shearvel = .true.
DEF_hist_vars%seddiag  = .true.
```

`seddiag` 是**history 输出开关**，不等于 `CoLMDEBUG` 编译选项。前者负责 NetCDF 过程变量，后者负责控制台高级日志、耗时、极值及内置站点信息。[history 实现](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L4768-L4938)

### 6.2 基础输出（`k=1..nsed`）

| 输出变量 | 单位 | 来源与含义 | 调参用途 |
|---|---|---|---|
| `f_sedcon_k` | m³/m³ | 按粒级、时间加权平均的固体体积分数 | 转 SSC；检查大小与季节性 |
| `f_sedout_k` | m³ solid/s | 悬浮泥沙有向输运通量 | 转悬浮输沙量 SSL，关注正负向 |
| `f_bedout_k` | m³ solid/s | 推移质有向输运通量 | 诊断床载负荷；不并入 SSC/悬沙 SSL |
| `f_sedinp_k` | m³ solid/s | 坡面泥沙输入 | 供沙是否过量 |
| `f_netflw_k` | m³ solid/s | 河床—水体净交换，**正值=净起悬**，负值=净入床 | 判断交换方向及直接入床贡献 |
| `f_layer_k` | m³ **bulk bed** | 活动层床沙体积（含孔隙） | 床层耗竭/补充；对比固体库存须乘 `1-λ` |
| `f_shearvel` | m/s | 摩阻速度 | 解释临界剪切、起悬与洪峰响应 |

### 6.3 `seddiag` 完整输出（`k=1..nsed`，共 10 组 × 粒级 + 1 总量）

| NetCDF 变量 | 单位 | 准确定义/限制 |
|---|---|---|
| `f_es_raw_k` | m³ solid/s | 床层起悬公式给出的原始 Es；尚未扣除床层可用量约束 |
| `f_d_raw_k` | m³ solid/s | 物理沉降 D；**已限制为不超过当前悬沙库存/dt** |
| `f_es_eff_k` | m³ solid/s | 库存限制后实际执行的总有效起悬 |
| `f_d_eff_k` | m³ solid/s | 所有有效入床：物理沉降 + 限幅/变干/保护区贡献 + 浅水外源直接入床 |
| `f_capdep_adv_k` | m³ solid/s | 正水体载体平流前/后的浓度上限转床 |
| `f_capdep_exch_k` | m³ solid/s | 床水交换之后的浓度上限转床 |
| `f_capdep_hill_k` | m³ solid/s | 湿河段坡面输入之后的浓度上限转床 |
| `f_capdep_other_k` | m³ solid/s | 其他分汊及可见/保护堤防的上限转床；相关路径未启用时为零 |
| `f_sedsto_k` | m³ **solid** | 本窗口按时间加权的**可见悬沙固体库存**；不含受保护悬沙 |
| `f_seddep_solid_k` | m³ **solid** | `(1-porosity)×sum(seddep_k, 所有命名层)`，不含活动层及保护床 |
| `f_ssc_precap_max` | mg/L | 本 history 窗口**限幅前**的最大**可见水体总 SSC**，跨粒级合计；非平均、非单粒级 |

**[S] 诊断时间约定**：

- 交换、平流及坡面输入的 **flux** 为 `sum(实际转移固体体积)/sum(dt)`。CapDep 在发生时已记录实际固体体积，**不能再乘一次 dt**。
- `f_es_raw/f_d_raw`、`f_es_eff/f_d_eff` 是平均 flux；`f_sedsto/f_seddep_solid/f_layer` 是每个形态步**床层重分配完成后的**时间加权状态采样。
- `f_ssc_precap_max` 是窗口**最大值**，不是平均值；没有正的可见水体载体时为缺测；缺少诊断累计时长时相关均值为缺测。
- **`CapDep` 属于 `D_eff` 的一部分，不可把 `f_capdep_*` 直接加到 `f_d_eff` 上再次计算总沉积。** `D_eff` 还包含 CapDep 没有计入的干载体/浅水直接入床。
- 上述 10 组每组每粒级一个变量；`nsed=3` 时合计 **31 个 `seddiag` 过程字段**。

依据：[history 变量与统计定义](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L194-L228)、[累积](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L1770-L1837)、[写出](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L4906-L4938)、[README](../../main/TRACER/README.md#sediment-process-history-diagnostics)。

### 6.4 单位换算：对接 OpenBench 的 Q / SSC / SSL [S]

对 `k=1..nsed` 求和（默认 `ρ_s=2650 kg/m³`），注意历史日均通量与日均浓度的区别：

```text
SSC [mg/L]   = sum_k(f_sedcon_k [m3/m3]) × rho_s [kg/m3] × 1000
             = sum_k(f_sedcon_k) × 2,650,000    (rho_s=2650)

SSL [t/day]  = sum_k(f_sedout_k [m3 solid/s]) × rho_s [kg/m3] × 86400 / 1000
             = sum_k(f_sedout_k) × 228,960     (rho_s=2650)
```

上述 SSL 为**有向悬浮输沙**。若评估只允许正的下游通量，必须先明确站点与河段流向定义、逆流处理和观测的 SSL 定义；不可随意对所有时段取绝对值。`f_bedout` 是推移质，不属于悬浮输沙量。每日 `Q×SSC` 通常**不等于**已按内部时间步聚合的 `f_sedout`，因为时间协方差、存储变化及不同聚合方法可能存在差异。

**[O]** 项目曾发现 OpenBench 结果的 SSC/SSL 大幅高估，但此前源码审计未发现上述体积单位在模型与 OpenBench 中被直接重复换算。仍应逐个实验核对实际 model registry/`compute` 公式与输出元数据，切勿仅凭历史讨论认定“单位永远正确”。

### 6.5 控制台日志与固定站点诊断 [S]

`CoLMDEBUG` 会增加 `Sediment timing`、`SED_PERF`（不同 worker 的 min/mean/max）、`Sediment precip diag`、极值 `MAX_SEDCON/MAX_SEDOUT`、守恒诊断及内置 4 个亚马孙参考站点的值：

- Serrinha
- Porto Velho
- Manacapuru
- Obidos

这些站点在 `initialize_sediment_diag_stations` 中采用**经度余弦修正的经纬距离最近格网中心**寻找 `ucat`，并经 MPI 选定编号。这个功能主要帮助日志定位，不带完整的河网连通性、上游面积/COMID/Pfaf 等水文匹配约束，不能直接视作正式验证的匹配证据。[源码](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L303-L403)

如需减少日志，仅调整对应 `CoLMDEBUG` 控制台 `WRITE` 分支；**不要**为了关闭站点日志误删 `seddiag` history、`sediment_diag_accumulate` 或进程间通信。`SED_NEAR_DRY_DEPTH` 打印并不等于该常数实际控制泥沙过程。

## 7. 亚马孙验证背景：事实、异常与调参解释

下面来自项目先前的运行/对话，**不是本手册所固定 commit 的复现实验**：

| 状态 | 观察/问题 | 下一步应索取的证据 |
|---|---|---|
| **[O]** | 2001–2002 日尺度、4 站 OpenBench 评估中，部分站点 Q 模拟尚可，但 SSC、SSL 有几十倍到更高的高估 | 保留每站 paired n、Q/SSC/SSL 配对序列、数据 QC、观测量纲 |
| **[O]** | 部分站点 SSC 模拟接近约 26,500 mg/L；曾报告 Porto Velho/Manacapuru 等站点存在显著高估 | 历史期 `f_ssc_precap_max`、`f_capdep_*`、`f_es_*/f_d_*` |
| **[H]** | 高估可能与产沙源过强、粒级沉速、河床反复起悬、有效水量过小/站点错配或时间配对相关 | 与干/湿季 `sedinp`、`Es`、`D`、`sedsto`、Q、深度联查 |
| **[T]** | 不应在查清因果前简单调低 `max_concentration` 作为“提高技能分”的手段 | 做参数扰动和物理过程预算的对照实验 |

**面向未来的最低限度评估**：Q、SSC、SSL 分开报告 KGE、NSE、PBIAS、相关性、有效样本数、均值比以及丰枯季表现；先检查 Q 水动力和站点—河网匹配，再解释 SSC/SSL。指标（尤其负 NSE/KGE）不能独立揭示物理原因，必须与 `Es/D`、CapDep 和库存联合看。

建议判读次序：

1. **站点与单位**：确认 OpenBench 采用对应河段、同一有效观测日、参考 QC 过滤，Q、SSC、SSL 的单位、求和粒级与有向流通量一致。
2. **水动力**：比对 Q、水体储量/水深、`shearvel`；如果 Q 本身偏差巨大，先不要将所有 SSL 偏差归咎于泥沙物理参数。
3. **坡面供沙**：`f_sedinp_k` 与降水、坡度、`yield_coefficient` 的响应；核查 `sed_slope` 输入分布。
4. **起悬/沉积**：比较 `Es_raw/Es_eff/D_raw/D_eff` 按粒级、按丰枯季的时间序列，检查床层库存受限是否主导。
5. **限幅机制**：`f_ssc_precap_max` 相比允许浓度、CapDep 的阶段组成及 `f_sedsto/f_layer/f_seddep_solid` 的演化。
6. **实验修改**：每轮原则上只改一个物理参数或一组明确耦合参数，保存 Git SHA、nml SHA、输入、restart、OpenBench 配置和成对评价结果。

**[H] 参数优先级不是按数值大小排序**：如果 `f_sedinp` 主导，先核实降水阈值/坡度/产沙系数；如果 `Es` 主导，先核实水动力、床沙粒级与起悬闭合；如果 `D_raw` 很高但 `D_eff` 与 CapDep 仍表现异常，先查床层限量、限幅时序及河床库存；`cfl_adv` 在本版不应被当作有效的物理调参目标。

## 8. Restart、版本与可复现性要求

**[S]** 当前模块在 restart 中保存 `sedsto`、`sedcon`、活动层、沉积层、必要的分汊/堤防库存与累积量，以及粒径、沉速、孔隙率、源项系数、`max_concentration`、床层深度、有效 `sed_frc`/`sed_slope` 等**参数与网格元数据**。读取时会校验配置，某些旧格式的非零悬沙浓度不能无损转换成当前主状态而被拒绝。开启 BIF/levee 时还使用特定 restart schema；不能因为两个实验使用同一个 restart 文件名就假定物理可兼容。[校验](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L3768-L3886)、[写出](../../main/TRACER/MOD_Tracer_Particle_Sediment.F90#L4400-L4467)

建议每轮实验至少记录：

```yaml
experiment_id: exp001
model_repo: jellywzx/CoLM202X_SYSU
model_commit: <exact SHA>
model_dirty_diff: <path or clean>
main_namelist: <path + sha256>
sediment_parameter_namelist: <actual resolved path + sha256>
history_namelist: <actual resolved path + sha256>
unitcatchment_file: <path + version/hash>
forcing: <product + time range + path/hash>
restart: <cold start or provenance + compatibility check>
history: <time resolution + period + mode>
station_matching: <method + river-network/area checks + mapping version>
reference: <dataset + QC rule + temporal pairing>
openbench: <commit/config + result directory>
metrics: <paired sample count + Q/SSC/SSL stationwise statistics>
process_diagnostics: <summary for Es/D, CapDep, storage, mass balance>
decision: <retain/reject hypothesis and next experiment>
```

不要把模型强制停止或质量守恒失败的产物当作成功完成的实验。运行脚本应保留真实退出状态，并在开始时打印真正运行的 Git commit、差异及**实际读取的**各个 namelist 和参数文件。

## 9. 代码已知边界与下一轮待核问题

| 标记 | 问题 | 下一轮工作 |
|---|---|---|
| **[S]** | CoLM 水动力先于泥沙，`sedsto`/泥沙床层不作为独立动力反馈去实时重塑河道地形 | 做形态反馈研究时必须单独开发/验证 |
| **[S]** | 有 BIF/levee 的保护区泥沙库存转移，但不等同于已经验证的完整洪泛平原侵蚀—沉积动力学 | 亚马孙应检查洪泛平原滞留/再悬浮的代表性 |
| **[S]** | 降水产沙为 unit-catchment 均一简化方案，`erosion_yield` 外部 forcing 当前是预留而未启用 | 对比空间显式侵蚀、土地利用与供沙连通性假设 |
| **[T]** | `sed_slope` 的原始单位及极值分布未由本文验证 | 审计 `grid_routing_data_15min_sed.nc` 及其构建脚本 |
| **[S]** | `SED_NEAR_DRY_DEPTH` 目前仅定义与打印；`cfl_adv` 当前不参与平流子步控制 | 不要据名称误判实际生效路径；修改代码后重新审计 |
| **[T]** | 历史 2001–2002 OpenBench 结果未与本版本 SHA/输入文件做完整锁定 | 将首轮基线独立写到 `CALIBRATION_VALIDATION.md`/`experiments/` |
| **[T]** | 内置亚马孙最近中心站点不等于河网一致匹配 | 引入同河段、河网连接和上游面积等匹配检验 |

## 10. 维护时建议执行的检查

已有测试入口可按代码修改范围选择运行；本文件**未声称这些测试在本次文档编辑中执行或通过**：

```bash
python -m pytest -q \
  tests/test_sediment_static.py \
  tests/test_sediment_max_concentration_paths.py \
  tests/test_sediment_critical_shear_runtime.py \
  tests/test_sediment_process_history.py \
  tests/test_sediment_substep_diagnostics.py \
  tests/test_sediment_hydro_handoff.py \
  tests/test_sediment_combined_ledger.py
```

MPI/BIF/levee 的行为另需对应测试、真实 MPI 作业与典型亚马孙站点时序验证，不能由静态测试代替。若 `model_commit`、参数、输出字段单位、restart schema 或关键算子顺序变化，优先更新 ``5–8` 与此文档的代码基线，再进行下一轮模型校准。

### 主要源码锚点（固定于本文件审计 SHA）

- [模块顶部与常数 1–257 行](https://github.com/jellywzx/CoLM202X_SYSU/blob/32ea6ce72b0a8ced75f5859d47c3d6cb3632cd07/main/TRACER/MOD_Tracer_Particle_Sediment.F90#L1-L257)
- [配置读取/校验 407–751 行](https://github.com/jellywzx/CoLM202X_SYSU/blob/32ea6ce72b0a8ced75f5859d47c3d6cb3632cd07/main/TRACER/MOD_Tracer_Particle_Sediment.F90#L407-L751)
- [形态学步进与调用顺序 1050–1363 行](https://github.com/jellywzx/CoLM202X_SYSU/blob/32ea6ce72b0a8ced75f5859d47c3d6cb3632cd07/main/TRACER/MOD_Tracer_Particle_Sediment.F90#L1050-L1363)
- [降水、静态参数与沉速 2375–2552 行](https://github.com/jellywzx/CoLM202X_SYSU/blob/32ea6ce72b0a8ced75f5859d47c3d6cb3632cd07/main/TRACER/MOD_Tracer_Particle_Sediment.F90#L2375-L2552)
- [临界剪切/平流/推移质 2676–3203 行](https://github.com/jellywzx/CoLM202X_SYSU/blob/32ea6ce72b0a8ced75f5859d47c3d6cb3632cd07/main/TRACER/MOD_Tracer_Particle_Sediment.F90#L2676-L3203)
- [床水交换、限幅与床层 3206–3501 行](https://github.com/jellywzx/CoLM202X_SYSU/blob/32ea6ce72b0a8ced75f5859d47c3d6cb3632cd07/main/TRACER/MOD_Tracer_Particle_Sediment.F90#L3206-L3501)
- [history 变量/单位/`seddiag` 4768–4938 行](https://github.com/jellywzx/CoLM202X_SYSU/blob/32ea6ce72b0a8ced75f5859d47c3d6cb3632cd07/main/TRACER/MOD_Tracer_Particle_Sediment.F90#L4768-L4938)

*文档 v0.1（2026-10-15）：以当前源码定义为主，整合项目 2026 年 9–10 月亚马孙调参与 OpenBench 诊断背景；后续实验结果请存放独立的调参/验证日志，不要覆盖这里的代码事实。*
