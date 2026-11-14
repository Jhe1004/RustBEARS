# Python 不确定性分析工作流

这里的短例子用于检查接口与数据流。小规模自举、单棵小树以及成功运行都不构成覆盖率、第一类错误率或功效验证。

## 拟合与概率契约

```python
import rustbears as rb

tree = "((T1:1,T2:1):1,((T3:.5,T4:.5):.5,(T5:.5,T6:.5):.5):1);"
geog = {"T1": "10", "T2": "11", "T3": "01",
        "T4": "10", "T5": "01", "T6": "11"}
options = dict(tree=tree, geog=geog,
               likelihood_conditioning="tip_survival",
               root_prior=rb.RootPriorSpec("uniform_nonnull"))
null_fit = rb.BioGeoBEARS(model="DEC", **options).fit()
alt_fit = rb.BioGeoBEARS(model="DEC+J", **options).fit()
```

`tip_survival`（A）是在全部叶端非空条件下的概率似然；`unconditional_complete`（B）包括内部空态向空态子节点确定性传递，需要显式设置 `allow_null_tips=True`。`legacy_unconditioned`（C）保留历史评分及其原始根权重尺度。C 的默认正向模拟会拒绝执行，必须指定生成条件；历史评分的自举可以在指定的生成分布下校准。

`RootPriorSpec("custom_weights", {"A": 2, "B": 1, "AB": .5})` 使用精确的状态标签；也可使用整数区域元组，例如 `{(0,): 2, (1,): 1, (0, 1): .5}`。A/B 对权重归一化，C 保留自定义权重尺度；采样使用单独的归一化概率。输入字典会被复制，后续改动不会影响已构建模型。

当前新契约覆盖单层标准六模型，以及距离指数 `+x`、人工扩散乘子指数 `+w`。分层、性状、固定内部观测、根边、化石 hook、模糊地理编码等组合会明确报错。

## 参数剖面

```python
profile = alt_fit.profile_ci(params=["d"], n_grid=8, num_threads=2)
print(profile.to_dataframe())
figure = profile.plot("profile.png")
surface = alt_fit.profile_surface("d", "j", grid_x=4, grid_y=4)
```

A/B 默认返回 `asymptotic_confidence_interval`；C 返回 `profile_support_interval`。默认阈值始终是双侧 χ²₁ 参考值，不因点估计触及零点自动换成混合分布阈值。`include_onesided_zero_ref=True` 另外给出单侧零点参考上限，不替代双侧区间。

下限 `at_physical_zero` 表示精确零点在支持集合内；`infeasible_zero_approached` 表示零点结构不可行，实际下限为正。`censored_at_box_upper` 和 `upper=None` 表示当前搜索盒没有找到交点，不能据此判定数学不可识别。返回的是包含最大值的连通区间；网格发现远端重新进入支持集合时，`disconnected_support_detected=True`。

`method="bootstrap_plugin"` 从全局点估计生成数据，在每次模拟中重拟合完整模型和固定原参数值的受限模型。它是局部插入式校准，不是逐候选值检验反演。自举校准失败会明确报错，不会静默删除失败样本后计算分位数。

## 模型自举、失败与续跑

```python
test = rb.bootstrap_lrt(null_fit, alt_fit, n_boot=19, num_threads=2,
                        seed=20261115, checkpoint_file="demo-lrt.jsonl")
print(test.summary())
# 同一任务中断后，保留原计划次数、种子、输入与优化设置：
# test = rb.bootstrap_lrt(null_fit, alt_fit, n_boot=19, num_threads=2,
#                        seed=20261115, checkpoint_file="demo-lrt.jsonl", resume=True)
```

示例的 19 次重复仅展示流程。每个重复实际重拟合零模型和备择模型，并把零模型解作为备择初值之一。历史模板的数值下限可在新自举流程中扩到物理零；用户自定义的严格正下限会拒绝不满足嵌套的比较。最终实际边界见 `test.effective_parameter_bounds`。

`pval_failure_bounds` 是缺失重复落在任意尾部位置的代数界；`pval_mc_ci95` 是 `[CP_L(r,B), CP_U(r+failed,B)]` 蒙特卡洛包络。即使没有失败，后者也不会退化成零宽度。所有失败留在 `replicates_df` 和计划分母里。

改变目标契约必须设置 `target_conditioning=...` 与 `refit_observed=True`，程序会先对相同观测数据重新拟合两个模型。C 保持评分时还必须给出 `simulation_conditioning="tip_survival"` 或 `"unconditional_complete"`。

检查点逐项核对输入、根权重、参数约束、精度、种子、软件及原生库等信息。变更任务必须使用新文件。`calc_surrogate_power=True` 给出的 `surrogate_power_fixed_threshold` 使用固定原零模型阈值；它不是完整嵌套自举检验的功效。

## 模型平均、树集合与随机映射

```python
average = rb.model_average([null_fit, alt_fit], criterion="AICc")
print(average.weights_df)
maps = average.run_bsm(num_maps=200, seed=7)
diagnostics = rb.bsm_diagnostics(maps)
adaptive = alt_fit.run_bsm_adaptive(min_maps=200, batch_size=200,
                                  max_maps=600, seed=8)
```

权重是频率学 AIC/AICc 权重，不是模型后验概率。各模型仍使用条件点估计。`fit_over_trees` 面向外部给定的树样本，地理数据不会反向更新这些树的后验权重；有缺失候选模型的树整棵隔离，并保留缺失质量界限。详细参数可通过 `help(rb.fit_over_trees)` 查看。

外层并行统一使用 `ParallelConfig(workers=..., native_threads=1)`；存在多个外层 worker 时禁止再开多线程内核。BLAS 在线程池整个生命周期统一限制为 1。线程数变化不会改变任务身份和随机流。

## 方法验证命令

脚本默认只写预注册清单，不运行模拟。`smoke` 与 `formal` 清单和目录不可混用；正式运行尚需项目成本和计算资源门禁放行。

```bash
python scripts/run_uncertainty_validation.py --analysis v1 --mode smoke \
  --scenario interior --output validation-smoke-v1 --execute
python scripts/run_uncertainty_validation.py --analysis v2b --mode smoke \
  --scenario null_e02 --output validation-smoke-v2b --execute
python scripts/run_uncertainty_validation.py --analysis v2c --mode smoke \
  --scenario dominant --output validation-smoke-v2c --execute

# 仅生成正式预注册，不启动正式计算：
python scripts/run_uncertainty_validation.py --analysis v1 --mode formal \
  --scenario double_zero --n-tips 40 --areas 4 --output planned-v1-double-zero
python scripts/run_uncertainty_validation.py --analysis v2b --mode formal \
  --scenario alt_j003 --finite-b-extrapolation --output planned-v2b-alt-j003
```

正式 V1 为每情景 200 个外层重复；V2b 为 200 个外层重复、主要 B=199，有限 B 对照使用最大 B=399 的共同内层样本；V2c 为每情景 300 个重复，保留真实内部节点状态用于 Brier、对数评分及可靠性分箱。固定阈值近似功效和完整嵌套拒绝率分别保存。使用同一批外层数据的有限 B 外推及其配对误差应单独报告，不能当作精确的无限 B 结论。[方法来源：Boos & Zhang (2000)](https://doi.org/10.1080/01621459.2000.10474226)。
