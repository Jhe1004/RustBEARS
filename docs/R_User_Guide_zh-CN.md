# RustBEARS R 软件包完整用户手册

<p align="left">
  <a href="R_User_Guide.md">English</a> | <strong>简体中文</strong> | <a href="../README_zh-CN.md">返回中文主页</a>
</p>

`RustBEARS` 完整保留了 `BioGeoBEARS` 的全部 **545** 个公开 R 函数、S4 类与绘图工具，底层由内存安全的编译型 Rust 引擎驱动，无需安装外部 Fortran（`rexpokit`）或 C++（`cladoRcpp`）扩展包。

---

## 1. 环境准备与安装

### 前置依赖
- **R** `>= 4.0.0`
- **Rust 工具链**（`rustc >= 1.82` 与 `cargo`）：可从 <https://rustup.rs/> 一键安装。
  - *Windows 特别提示*：Rtools 使用 GNU 工具链，首次在 Windows 上编译前请运行一次 `rustup target add x86_64-pc-windows-gnu`。

### 从 GitHub 直接安装

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_github("Jhe1004/RustBEARS", subdir = "r_package/RustBEARS")
```

*(说明：`MultinomialCI` 包已在 CRAN 归档，`RustBEARS` 已在 `DESCRIPTION` 的 `Remotes` 字段中内置其 CRAN 官方归档源地址，使用 `remotes::install_github()` 或 `pak::pkg_install("Jhe1004/RustBEARS/r_package/RustBEARS")` 时会自动下载安装。)*

### 全局精度与线程选项

```r
library(RustBEARS)

# 数值精度模式："StrictF64"（默认，64 位 IEEE-754 双精度）或 "LegacyF32"（模拟原版 32 位混合精度舍入）
options(RustBEARS.precision = "StrictF64")

# 可选：指定底层 Rust 并行线程数
options(RustBEARS.threads = 4L)
```

---

## 2. 原版 BioGeoBEARS 脚本零改动迁移

您现有的任何 `BioGeoBEARS` 分析脚本都无需重写——只需将开头的 `library(BioGeoBEARS)` 替换为 `library(RustBEARS)` 即可直接运行：

```r
library(RustBEARS)

trfn   <- system.file("extdata/Psychotria_5.2_ultrametric.newick", package = "RustBEARS")
geogfn <- system.file("extdata/Psychotria_geog.data", package = "RustBEARS")

BioGeoBEARS_run_object <- define_BioGeoBEARS_run(
  trfn = trfn,
  geogfn = geogfn,
  max_range_size = 4,
  num_cores_to_use = 4L
)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

res_DEC <- bears_optim_run(BioGeoBEARS_run_object)
plot_BioGeoBEARS_results(res_DEC, plotwhat = "pie")
```

默认情况下，`BioGeoBEARS_run_object` 采用 `likelihood_conditioning = "legacy_unconditioned"`，确保原有脚本产出与原版一致的未条件化对数似然与祖先状态概率。

---

## 3. 现生末端存活条件似然与六大模型配置

针对所有采样物种均为现生物种的超度量系统发育树，只需在优化前设置 `likelihood_conditioning <- "tip_survival"` 与 `define_root_prior()`，即可开启精确的末端存活条件似然，消除空态吸收造成的零灭绝偏差：

```r
run <- define_BioGeoBEARS_run(
  trfn = trfn,
  geogfn = geogfn,
  max_range_size = 4,
  calc_ancprobs = TRUE,
  print_optim = FALSE
)
run$likelihood_conditioning <- "tip_survival"
run$root_prior <- define_root_prior("uniform_nonnull")  # 或使用 custom_weights 自定义根状态权重
run <- readfiles_BioGeoBEARS_run(run)

# 1. DEC 模型
fit_dec <- bears_optim_run(run)

# 2. DEC+J 模型
run_decj <- run
run_decj$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
run_decj$BioGeoBEARS_model_object@params_table["j", c("init", "est")] <- 0.01
fit_decj <- bears_optim_run(run_decj)

# 3. DIVALIKE 模型（若再将 j 设为 free 即为 DIVALIKE+J）
run_diva <- run
run_diva$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
run_diva$BioGeoBEARS_model_object@params_table["s", c("init", "est")] <- 0.0
run_diva$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "2-j"
run_diva$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/2"
run_diva$BioGeoBEARS_model_object@params_table["v", "type"] <- "ysv*1/2"
run_diva$BioGeoBEARS_model_object@params_table["mx01v", "type"] <- "fixed"
run_diva$BioGeoBEARS_model_object@params_table["mx01v", c("init", "est")] <- 0.5
fit_diva <- bears_optim_run(run_diva)

# 4. BAYAREALIKE 模型（若再将 j 设为 free 即为 BAYAREALIKE+J）
run_bay <- run
run_bay$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
run_bay$BioGeoBEARS_model_object@params_table["s", c("init", "est")] <- 0.0
run_bay$BioGeoBEARS_model_object@params_table["v", "type"] <- "fixed"
run_bay$BioGeoBEARS_model_object@params_table["v", c("init", "est")] <- 0.0
run_bay$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "1-j"
run_bay$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/1"
run_bay$BioGeoBEARS_model_object@params_table["y", "type"] <- "1-j"
run_bay$BioGeoBEARS_model_object@params_table["mx01y", "type"] <- "fixed"
run_bay$BioGeoBEARS_model_object@params_table["mx01y", c("init", "est")] <- 0.9999
fit_bay <- bears_optim_run(run_bay)
```

### 按年代涌现岛屿的多时间层建模

当岛屿或区域由古老地质时期向现代逐步涌现（即较老时间层允许的区域在较年轻时间层中继续保留）时，`tip_survival` 直接支持时间分层（`timeperiods`）、允许区域掩膜（`areas_allowed`）、邻接矩阵及距离修饰因子：

```r
run_strat <- run
run_strat$timesfn <- system.file("extdata/examples/Psychotria_M3strat/BGB/timeperiods.txt", package = "RustBEARS")
run_strat$areas_allowed_fn <- system.file("extdata/examples/Psychotria_M3strat/BGB/areas_allowed.txt", package = "RustBEARS")
run_strat <- readfiles_BioGeoBEARS_run(run_strat)
run_strat <- section_the_tree(inputs = run_strat, make_master_table = TRUE, plot_pieces = FALSE)

fit_strat_dec <- bears_optim_run(run_strat)
```

---

## 4. 全树均匀化（Uniformization）随机建图

调用 `runBSM()` 并指定 `bsm_engine = "native_whole_tree"`，即可使用免除拒绝试次上限的全树端点条件均匀化引擎：

```r
# 抽取固定张数（例如 500 张随机演化历史）
maps <- runBSM(
  fit_dec, NULL,
  bsm_engine = "native_whole_tree",
  maxnum_maps_to_try = 500L,
  seedval = 42L
)

# 或采用自适应精度模式（达到目标相对蒙特卡洛标准误 target_rmcse 后自动停止）
maps_adaptive <- runBSM(
  fit_dec, NULL,
  bsm_engine = "native_whole_tree",
  target_rmcse = 0.05,
  min_maps = 200L,
  batch_size = 200L,
  max_maps = 5000L,
  seedval = 42L
)
```

---

## 5. 五大不确定性量化函数速查与示例

| 函数 | 核心用途 | 常用关键参数 |
| :--- | :--- | :--- |
| `profile_ci(fit, ...)` | 剖面似然置信区间（自动区分零边界最优、结构不可行零点与搜索上限平坦脊右删失） | `params = c("d", "e")`, `level = 0.95`, `n_grid = 16L` |
| `profile_surface(fit, ...)` | 二维联合剖面似然曲面扫描 | `param_x = "d"`, `param_y = "e"`, `grid_x = 12L`, `grid_y = 12L` |
| `bootstrap_lrt(null_fit, alt_fit, ...)` | 参数自助法似然比检验（含零模型热启动、保守 p 值、蒙特卡洛包络与断点续跑） | `n_boot = 199L`, `seedval = 42L`, `checkpoint_file = "lrt.rds"`, `resume = FALSE` |
| `model_average(fits, ...)` | 按 AIC/AICc 权重整合候选模型的祖先状态概率与随机建图配额 | `criterion = "AICc"`, `num_maps = 500L`, `seedval = 42L` |
| `fit_over_trees(run, trees, ...)` | 跨外部后验演化树集重拟合模型与随机建图，输出演化支单系支持率与缺失拟合上下界 | `models = c("DEC", "DEC+J")`, `likelihood_conditioning = "tip_survival"`, `bsm_maps_per_tree = 200L`, `checkpoint_file = "trees.rds"` |
| `bsm_diagnostics(bsm_obj, ...)` | 计算事件计数与路线流量的 MCSE、独立树 $t$ 区间、MCMC 有效样本量（`ESS >= 50`）检查及路线等效性检验 | `tree_target = "fixed_sample"`（或 `"posterior_population"`）, `tree_sampling = "independent"` / `"mcmc"`, `equivalence_delta = 0.25` |

### 完整工作流示例

```r
# 1. 剖面似然 95% 置信区间及可视化
prof <- profile_ci(fit_dec, params = c("d", "e"), n_grid = 16L)
print(as.data.frame(prof))
plot(prof)

# 2. 参数自助法似然比检验（中断后设 resume = TRUE 即可从 checkpoint_file 续跑）
lrt <- bootstrap_lrt(
  fit_dec, fit_decj,
  n_boot = 99L,
  seedval = 2026L,
  checkpoint_file = tempfile(fileext = ".rds")
)
print(lrt)

# 3. Akaike 多模型加权平均与随机建图精度诊断
avg <- model_average(
  list(DEC = fit_dec, DECJ = fit_decj),
  criterion = "AICc",
  num_maps = 400L,
  seedval = 42L
)
diag <- bsm_diagnostics(avg$joint_bsm, tree_target = "fixed_sample", equivalence_delta = 0.25)
print(diag$event_summary)
```

在 R 控制台中可随时输入 `?profile_ci`、`?bootstrap_lrt`、`?model_average`、`?fit_over_trees`、`?bsm_diagnostics` 查看内置帮助文档，或直接运行仓库自带脚本 `source("examples/quickstart_r.R")` 与 `source("examples/uncertainty_workflow.R")`。
