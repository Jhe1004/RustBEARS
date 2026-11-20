# RustBEARS Python API 与命令行完整用户手册

<p align="left">
  <a href="Python_and_CLI_Guide.md">English</a> | <strong>简体中文</strong> | <a href="../README_zh-CN.md">返回中文主页</a>
</p>

`rustbears` Python 软件包（要求 `Python >= 3.10`）提供简洁直观的面向对象分析接口、多进程/多线程并行调度、支持中断续跑的 JSONL 检查点，以及开箱即用的 `rustbears` 命令行工具。

---

## 1. 安装方法

请先确保系统已安装 Rust 工具链（`rustc >= 1.82` 与 `cargo`，见 <https://rustup.rs/>），然后执行：

```bash
python3 -m pip install "git+https://github.com/Jhe1004/RustBEARS.git"
```

或者克隆仓库到本地后安装：

```bash
git clone https://github.com/Jhe1004/RustBEARS.git
cd RustBEARS
python3 -m pip install .
```

---

## 2. 模型拟合与现生末端存活条件似然

`rb.BioGeoBEARS`（同义别名 `rb.RustBEARS`）的 `tree` 参数既可传入 Newick 文件路径也可直接传入 Newick 字符串；`geog` 参数既可传入 PHYLIP 格式地理分布文件路径，也可直接传入 Python 字典。

```python
import rustbears as rb

tree = "((K_sp1:2.0,O_sp1:2.0):2.0,((O_sp2:1.0,M_sp1:1.0):1.5,(M_sp2:1.0,H_sp1:1.0):1.5):1.5);"
geog = {
    "K_sp1": "1000", "O_sp1": "0100", "O_sp2": "0100",
    "M_sp1": "0010", "M_sp2": "0010", "H_sp1": "0001",
}

# 内置模型名："DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"
common_opts = dict(
    tree=tree,
    geog=geog,
    max_range_size=4,
    precision="strict_f64",                 # 或设为 "legacy_f32" 复现 32 位历史混合精度
    likelihood_conditioning="tip_survival", # 或 "legacy_unconditioned"（默认） / "unconditional_complete"
    root_prior=rb.RootPriorSpec("uniform_nonnull"),
)

fit_dec  = rb.BioGeoBEARS(model="DEC", **common_opts).fit(num_threads=4)
fit_decj = rb.BioGeoBEARS(model="DEC+J", **common_opts).fit(num_threads=4)

print(rb.compare_models([fit_dec, fit_decj]))
print(fit_dec.summary())
fit_dec.plot("ancestral_dec.pdf", plot_what="pie")
```

### 空间距离修饰因子与多时间层岛屿涌现建模

可直接向 `rb.BioGeoBEARS` 传入时间切片与空间矩阵文件路径：

```python
strat_fit = rb.BioGeoBEARS(
    tree="Psychotria_5.2_ultrametric.newick",
    geog="Psychotria_geog.data",
    model="DEC+J",
    max_range_size=4,
    timeperiods="timeperiods.txt",
    areas_allowed="areas_allowed.txt",          # tip_survival 原生支持区域按年代单调涌现！
    distances="distances.txt",                  # 可选：+x 地理距离依赖扩散
    dispersal_multipliers="multipliers.txt",    # 可选：手动扩散乘子矩阵
    likelihood_conditioning="tip_survival",
    root_prior=rb.RootPriorSpec("uniform_nonnull"),
).fit(num_threads=4)
```

---

## 3. 多线程确定性均匀化随机建图

每个拟合结果对象 `FitResult`（以及多模型平均结果 `ModelAverageResult`）均可直接执行全树均匀化随机建图：

```python
# 抽取固定张数（无论 num_threads 设为几核，每张地图按编号均严格确定性一致）
maps = fit_dec.run_bsm(num_maps=500, num_threads=4, seed=42)
print(maps.summary_counts["event_means"])
rb.plot_bsm(maps, "bsm_events.pdf")

# 按目标相对蒙特卡洛标准误（target_rmcse）自动追加抽样的自适应模式
maps_adap = fit_dec.run_bsm_adaptive(
    target_rmcse=0.05,
    min_maps=200,
    batch_size=200,
    max_maps=2000,
    num_threads=4,
    seed=42,
)
```

---

## 4. 五大频率学派不确定性量化工作流

### 4.1 剖面似然置信区间与二维似然曲面

```python
prof = fit_dec.profile_ci(params=["d", "e"], n_grid=16, num_threads=4)
print(prof.to_dataframe())
prof.plot("profile_dec.png")

surf = fit_dec.profile_surface("d", "e", grid_x=12, grid_y=12, num_threads=4)
```

- 在 `prof.to_dataframe()` 的 `status_lower` 与 `status_upper` 列中可查看边界诊断状态（如 `at_physical_zero` 物理零边界最优、`infeasible_zero_approached` 不可行零点、`censored_at_box_upper` 搜索上限平坦似然脊右删失预警）。

### 4.2 参数自助法似然比检验（支持断点续跑）

```python
lrt = rb.bootstrap_lrt(
    fit_dec,
    fit_decj,
    n_boot=199,
    seed=2026,
    checkpoint_file="bootstrap_dec_vs_decj.jsonl",
    parallel=rb.ParallelConfig(workers=4, native_threads=1),
)
print(lrt.summary())
# 如中途中断，保持参数一致并加上 `resume=True` 即可无缝续跑！
```

### 4.3 Akaike 多模型加权平均与建图精度诊断

```python
avg = rb.model_average([fit_dec, fit_decj], criterion="AIC")  # 样本量充足时推荐 "AICc"
print(avg.weights_df)

avg_maps = avg.run_bsm(num_maps=400, num_threads=4, seed=42)
diag = rb.bsm_diagnostics(avg_maps, tie_margin=0.25, tree_target="fixed_sample")
print(diag.event_summary_df)
```

### 4.4 外部后验演化树集不确定性传递（`fit_over_trees`）

将来自 BEAST 或 RevBayes 等软件的后验系统发育树集传入 `rb.fit_over_trees`，一次性传递拓扑、模型、参数与随机建图误差：

```python
multi = rb.fit_over_trees(
    trees="posterior_sample.trees",   # Newick/Nexus 文件路径或 Newick 字符串列表
    geog="geography.data",
    models=["DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"],
    summary_tree="mcc_reference.newick",
    likelihood_conditioning="tip_survival",
    root_prior=rb.RootPriorSpec("uniform_nonnull"),
    criterion="AICc",
    bsm_maps_per_tree=200,
    checkpoint_file="posterior_propagation.jsonl",
    parallel=rb.ParallelConfig(workers=4, native_threads=1),
)
print(multi.clade_summary_df)

# 计算面向后验树总体的建图精度区间（MCMC 树链自动启用 Geyer IPS 有效样本量 ESS >= 50 门禁）
pop_diag = rb.bsm_diagnostics(
    multi.bsm_result,
    tree_target="posterior_population",
    tree_dependence="mcmc",  # 若为独立抽样树集则设为 "independent"（采用 Student-t 区间）
    tie_margin=0.25,
)
```

---

## 5. 命令行工具（`rustbears` CLI）

无需编写任何代码，直接在终端中完成多模型拟合、模型比较表导出、祖先状态绘图与随机建图：

```bash
rustbears fit \
  --tree tree.newick \
  --geog geography.data \
  --models DEC,DEC+J,DIVALIKE,DIVALIKE+J,BAYAREALIKE,BAYAREALIKE+J \
  --max-range-size 4 \
  --likelihood-conditioning tip_survival \
  --precision strict_f64 \
  --maps 200 \
  --threads 4 \
  --seed 12345 \
  --out-dir results_out
```

### `rustbears fit` 参数一览

| 参数 | 说明 |
| :--- | :--- |
| `--tree PATH` | 输入 Newick 系统发育树文件（必填） |
| `--geog PATH` | 输入 PHYLIP 格式离散地理分布文件（必填） |
| `--models LIST` | 逗号分隔的模型列表，如 `DEC,DEC+J,DIVALIKE,BAYAREALIKE`（默认：`DEC`） |
| `--max-range-size INT` | 可选的最大同时占据区域数上限 |
| `--likelihood-conditioning` | 似然形式：`legacy_unconditioned`（默认）、`tip_survival` 或 `unconditional_complete` |
| `--precision` | 精度模式：`strict_f64`（默认 64 位双精度）或 `legacy_f32`（32 位历史混合精度） |
| `--use-ambiguities` | 展开末端地理分布中的 `?` 通配模糊编码 |
| `--timeperiods PATH` | 可选的地质时间层边界文件 |
| `--areas-allowed PATH` | 可选的各时间层允许区域矩阵 |
| `--areas-adjacency PATH` | 可选的各时间层区域邻接矩阵 |
| `--distances PATH` | 可选的区域间地理距离矩阵（自动启用 `+x` 指数参数） |
| `--envdistances PATH` | 可选的环境距离矩阵（自动启用 `+n` 指数参数） |
| `--dispersal-multipliers PATH` | 可选的手动扩散乘子矩阵（自动启用 `+w` 指数参数） |
| `--maps INT` | 每个模型抽取的均匀化随机建图张数（默认 `0` 表示不跑建图） |
| `--threads INT` | 并行 CPU 线程数（默认：`1`） |
| `--plot {pie,text}` | 祖先状态重建 PDF 图形样式（默认：`pie`） |
| `--out-dir PATH` | 输出目录（自动生成 `model-comparison.csv` 及各模型子目录图表） |
