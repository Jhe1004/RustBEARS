# RustBEARS

<p align="left">
  <a href="README.md">English</a> | <strong>简体中文</strong>
</p>

[![License: GPL v3+](https://img.shields.io/badge/License-GPL%20v3%2B-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/Release-v1.2.0-success.svg)](CHANGELOG.md)
[![Rust Core](https://img.shields.io/badge/Rust-1.82%2B-orange.svg)](rust/Cargo.toml)
[![R Package](https://img.shields.io/badge/R-4.0%2B-276DC3.svg)](r_package/RustBEARS/DESCRIPTION)
[![Python Package](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](pyproject.toml)

**RustBEARS** 是一套面向概率历史生物地理学的内存安全、多核并行 **Rust 高性能计算引擎**，同时提供 **与 BioGeoBEARS 完全兼容的 R 软件包**（`RustBEARS`）以及 **原生 Python 软件包与命令行工具**（`rustbears`），无需编译安装外部 Fortran 77（`rexpokit`）或 C++（`cladoRcpp`）扩展包。

---

## 核心特性

1. **100% 兼容 BioGeoBEARS 且大幅提速**
   - 完整保留 [`BioGeoBEARS`](http://phylo.wikidot.com/biogeobears)（Matzke 2013, 2014）在 R 端的全部 **545 个公开函数、S4 对象结构与绘图工具**——您现有的 R 分析脚本只需把 `library(BioGeoBEARS)` 替换为 `library(RustBEARS)` 即可直接运行，无需重写代码。
   - 在来自 77 项已发表研究的 1,339 组配置中精确复现原版对数似然（中位偏差 `9.95e-14`）与祖先状态概率（`<= 1.50e-14`），同时将 64 至 256 状态的最大似然优化在单核下提速 **2.8–17.3 倍**、在八核下相对原版最快配置提速 **4.1–36.1 倍**（严格八核对八核比较达 **10.9–205.5 倍**）。
2. **现生末端存活条件似然（Tip-Survival Conditioning）**
   - 修正在单层及岛屿按年代涌现的多时间层模型中因未扣除谱系中途灭绝概率而导致的“零灭绝偏差”与创始者事件物种形成（`+J`）假阳性偏高问题，同时完整保留单区域物种的沿枝局地绝迹过程。
3. **无试次上限的均匀化（Uniformization）随机建图**
   - 用端点条件均匀化马尔可夫桥取代传统的逐枝正向试错拒绝抽样（30 项实证研究中位数提速 **170.9 倍**），彻底解决新涌现岛屿或罕见垫脚石路径上因接受率极低而卡死失败的问题，并保证每张随机历史按编号跨多线程严格确定性可复现。
4. **一体化频率学派不确定性量化套件**
   - 内置五大诊断模块：**剖面似然置信区间**（含零边界探测与高灭绝平坦似然脊预警）、**参数自助法似然比检验**（支持断点续跑）、**Akaike 多模型加权平均**、**外部贝叶斯后验演化树集传递**，以及 **随机建图蒙特卡洛标准误（MCSE）与 MCMC 有效样本量（ESS）精度检验**。

---

## 文档导航

| 文档 | 简体中文 | English |
| :--- | :--- | :--- |
| **仓库主页与快速入门** | **README_zh-CN.md**（本页） | [README.md](README.md) |
| **R 软件包完整用户手册**（`library(RustBEARS)`） | [docs/R_User_Guide_zh-CN.md](docs/R_User_Guide_zh-CN.md) | [docs/R_User_Guide.md](docs/R_User_Guide.md) |
| **Python API 与命令行完整手册**（`import rustbears` / `rustbears`） | [docs/Python_and_CLI_Guide_zh-CN.md](docs/Python_and_CLI_Guide_zh-CN.md) | [docs/Python_and_CLI_Guide.md](docs/Python_and_CLI_Guide.md) |
| **统计似然形式与模型适用边界** | [docs/Statistical_Methods_and_Contracts.md](docs/Statistical_Methods_and_Contracts.md) | [docs/Statistical_Methods_and_Contracts.md](docs/Statistical_Methods_and_Contracts.md) |
| **可直接运行的示例脚本** | [`examples/quickstart_r.R`](examples/quickstart_r.R) · [`examples/uncertainty_workflow.R`](examples/uncertainty_workflow.R) · [`examples/quickstart_python.py`](examples/quickstart_python.py) | 同左 |

---

## 安装指南

### 准备工作：安装 Rust 编译器工具链

首次安装前，请先通过 [rustup.rs](https://rustup.rs/) 安装标准 Rust 工具链（要求 `rustc >= 1.82` 与 `cargo`）：

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

*(本仓库已内置全部第三方 Rust 依赖源码 `rust/vendor`，编译时无需联网下载 Rust 依赖。若在 Windows 上配合 Rtools 编译 R 包，请先运行一次 `rustup target add x86_64-pc-windows-gnu`。)*

---

### 方式一：安装 R 软件包（`RustBEARS`，要求 R >= 4.0）

**在 R 中直接从 GitHub 一键安装**（会自动处理 CRAN 依赖及已归档的 `MultinomialCI` 包）：

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_github("Jhe1004/RustBEARS", subdir = "r_package/RustBEARS")
```

*(也可以使用 `pak` 安装：`pak::pkg_install("Jhe1004/RustBEARS/r_package/RustBEARS")`。)*

**或者克隆仓库到本地后安装：**

```bash
git clone https://github.com/Jhe1004/RustBEARS.git
cd RustBEARS
R CMD INSTALL r_package/RustBEARS
```

---

### 方式二：安装 Python 软件包与命令行工具（`rustbears`，要求 Python >= 3.10）

**直接从 GitHub 一键安装：**

```bash
python3 -m pip install "git+https://github.com/Jhe1004/RustBEARS.git"
```

**或者克隆仓库到本地后安装：**

```bash
git clone https://github.com/Jhe1004/RustBEARS.git
cd RustBEARS
python3 -m pip install .
```

验证命令行工具是否就绪：

```bash
rustbears fit --help
```

---

## 5 分钟快速上手

### 1. R 语言快速示例（使用内置的夏威夷茜草属 *Psychotria* 数据集）

```r
library(RustBEARS)

# 1. 读取包内自带的示例演化树与地理分布文件
trfn   <- system.file("extdata/Psychotria_5.2_ultrametric.newick", package = "RustBEARS")
geogfn <- system.file("extdata/Psychotria_geog.data", package = "RustBEARS")

run <- define_BioGeoBEARS_run(
  trfn = trfn,
  geogfn = geogfn,
  max_range_size = 4,
  num_cores_to_use = 1L,
  print_optim = FALSE
)

# 2. 可选：启用现生末端存活条件似然（若不设置此项，默认保持原版 BioGeoBEARS 历史兼容模式）
run$likelihood_conditioning <- "tip_survival"
run$root_prior <- define_root_prior("uniform_nonnull")
run <- readfiles_BioGeoBEARS_run(run)

# 3. 拟合 DEC 与 DEC+J 模型
fit_dec <- bears_optim_run(run)

run_j <- run
run_j$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
run_j$BioGeoBEARS_model_object@params_table["j", c("init", "est")] <- 0.01
fit_decj <- bears_optim_run(run_j)

# 4. 计算剖面似然 95% 置信区间与 Akaike 多模型加权平均
ci_dec <- profile_ci(fit_dec, params = c("d", "e"))
print(ci_dec)

avg <- model_average(list(DEC = fit_dec, DECJ = fit_decj), criterion = "AICc", num_maps = 200L, seedval = 42L)
print(avg$weights)
```

👉 原版脚本无缝迁移、岛屿涌现时间分层建模、全树均匀化随机建图 `runBSM()`、参数自助法 `bootstrap_lrt()` 及后验树集传递 `fit_over_trees()` 的完整说明请见 **[docs/R_User_Guide_zh-CN.md](docs/R_User_Guide_zh-CN.md)**。

---

### 2. Python API 与命令行快速示例

```python
import rustbears as rb

tree = "((K_sp1:2.0,O_sp1:2.0):2.0,((O_sp2:1.0,M_sp1:1.0):1.5,(M_sp2:1.0,H_sp1:1.0):1.5):1.5);"
geog = {
    "K_sp1": "1100", "O_sp1": "0100", "O_sp2": "0110",
    "M_sp1": "0010", "M_sp2": "0011", "H_sp1": "0001",
}

# 在现生末端存活条件似然下并行拟合六大经典模型
models = ["DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"]
fits = [
    rb.BioGeoBEARS(
        tree=tree,
        geog=geog,
        model=m,
        max_range_size=4,
        likelihood_conditioning="tip_survival",
        root_prior=rb.RootPriorSpec("uniform_nonnull"),
    ).fit(num_threads=4)
    for m in models
]

print(rb.compare_models(fits))

# Akaike 多模型加权均匀化随机建图与精度诊断
avg = rb.model_average(fits, criterion="AIC")
maps = avg.run_bsm(num_maps=200, num_threads=4, seed=42)
diag = rb.bsm_diagnostics(maps)
print(maps.summary_counts["event_means"])
```

**命令行一键批量拟合与随机建图：**

```bash
rustbears fit \
  --tree tree.newick \
  --geog geography.data \
  --models DEC,DEC+J,DIVALIKE,DIVALIKE+J,BAYAREALIKE,BAYAREALIKE+J \
  --max-range-size 4 \
  --likelihood-conditioning tip_survival \
  --maps 200 \
  --threads 4 \
  --out-dir results_biogeography
```

👉 Python 全部参数详解、时间分层与距离矩阵、断点续跑及 CLI 用法请见 **[docs/Python_and_CLI_Guide_zh-CN.md](docs/Python_and_CLI_Guide_zh-CN.md)**。

---

## 浮点精度与三种似然形式

- **浮点运算精度**：默认采用 64 位 IEEE-754 双精度（R 端 `StrictF64`，Python 端 `precision="strict_f64"`）。如需逐位复现历史文献中 `cladoRcpp` / `EXPOKIT` 的 32 位混合精度舍入行为，可在 R 中设置 `options(RustBEARS.precision = "LegacyF32")` 或在 Python 中传入 `precision="legacy_f32"`。
- **三种似然形式**：
  1. `"legacy_unconditioned"` *（向后兼容默认项）*：与原版 `BioGeoBEARS` 完全一致的未条件化剪枝似然（支持含化石末端树、非单调区域沉没掩膜及性状耦合模型）。
  2. `"tip_survival"`：对“所有采样现生末端均存活至今且占据非空分布区”作精确条件化校正（支持六大模型族、各类空间/距离/面积修饰因子，以及区域按年代单调涌现的多时间层模型）。
  3. `"unconditional_complete"`：在全部 $S^N$ 种树端结局（含空态末端）上严格归一化的生成概率分布，用于全状态空间正向模拟与概率加和核验。

---

## 仓库目录结构

```text
├── rust/                       # 底层共享 Rust 计算引擎 (bgb-core, bgb-numerics, bgb-clado, bgb-Step, bgb-bsm, bgb-ffi)
├── r_package/RustBEARS/        # BioGeoBEARS 兼容 R 软件包 (`library(RustBEARS)`)
├── rustbears/                  # 原生 Python 软件包与命令行工具 (`import rustbears` / `rustbears`)
├── docs/                       # R、Python、命令行及统计方法的中英双语完整文档
├── examples/                   # 开箱即用的 R 与 Python 快速示例及不确定性脚本
├── tests/                      # Rust、R 与 Python 回归测试及统计校准测试集
├── scripts/                    # 源码同步、发布打包与基准测试脚本
├── CITATION.cff                # 学术共同引用元数据
└── LICENSE                     # GNU 通用公共许可证 v3.0 或更高版本 (GPL-3.0-or-later)
```

---

## 学术引用与开源协议

`RustBEARS` 遵循 **GNU 通用公共许可证 v3.0 或更高版本（`GPL-3.0-or-later`）** 开源发布（详见 [`LICENSE`](LICENSE)、[`NOTICE-Rust-port.md`](NOTICE-Rust-port.md) 与 [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md)）。

使用 `RustBEARS` 发表研究成果时，请**同时引用** Nicholas J. Matzke 的 `BioGeoBEARS` 奠基论文与 `RustBEARS`（可在 R 中运行 `citation("RustBEARS")` 或查阅 [`CITATION.cff`](CITATION.cff)）：

- Matzke, N. J. (2013). Probabilistic historical biogeography: new models for founder-event speciation, imperfect detection, and fossils allow improved accuracy and model-testing. *Frontiers of Biogeography*, 5(4), 242–248. https://doi.org/10.21425/F5FBG19694
- Matzke, N. J. (2014). Model selection in historical biogeography reveals that founder-event speciation is a crucial process in island clades. *Systematic Biology*, 63(6), 951–970. https://doi.org/10.1093/sysbio/syu056
- He, J. (2026). RustBEARS: Survival-conditioned likelihoods, Uniformization stochastic mapping, and uncertainty quantification for historical biogeography (v1.2.0). https://github.com/Jhe1004/RustBEARS

如遇安装、编译或运行问题，欢迎在 [GitHub Issues](https://github.com/Jhe1004/RustBEARS/issues) 提交反馈。
