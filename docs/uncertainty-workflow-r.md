# RustBEARS v1.2 R 不确定性工作流

这是开发版教程。运行示例检查接口，不代表完成了参数覆盖率或经验数据的正式推断。
完整可执行示例是 `examples/uncertainty_workflow.R`。安装本地R包后，在仓库根目录执行：

```r
source('examples/uncertainty_workflow.R')
```

## 先声明似然口径

在 `define_BioGeoBEARS_run()` 返回的对象上设置
`run$likelihood_conditioning <- 'tip_survival'` 和
`run$root_prior <- define_root_prior('uniform_nonnull')`，再读取文件和拟合。
Contract A表示所有叶端非空的条件概率；Contract B使用
`'unconditional_complete'`，还需设置 `run$allow_null_tips <- TRUE`；
Contract C的 `'legacy_unconditioned'` 保留历史评分尺度。

v1.2的R默认精度是StrictF64；复查旧结果时显式使用
`options(RustBEARS.precision='LegacyF32')`。模型平均、自举和剖面分析不会
把不同契约或不同根权重的结果暗中混用。已拟合后修改输入文件，需要重新拟合。

## 剖面与自举

```r
intervals <- profile_ci(dec, params='e', include_onesided_zero_ref=TRUE)
test <- bootstrap_lrt(dec, decj, n_boot=399L, seedval=42L,
                      checkpoint_file='bootstrap.rds')
test$pval_failure_bounds
test$pval_mc_ci95
```

每个剖面点重优化伴随参数。端点状态区分物理零、结构不可行零和搜索盒删失。
C契约得到支持区间；A/B的渐近或局部插入式区间也不能保证边界小样本覆盖率。
自举每次重拟合零模型和备择模型；失败重复保留在计划分母中。
`pval_failure_bounds`是优化失败造成的缺失界，`pval_mc_ci95`是有限模拟次数的
Clopper–Pearson包络，两者不能互相替代。固定阈值surrogate不是嵌套检验功效。
从C切换到A/B，必须设置 `refit_observed=TRUE`，程序会真正重拟合观测数据。
C生成模拟必须显式设置 `simulation_conditioning`。

## 模型与树集合

`model_average(list(dec,decj), criterion='AICc', num_maps=500L)`在同一
数据、树和契约下平均。AICc使用末端数作为明确的样本量约定；若模型参数过多
使其未定义，使用明确选择的AIC。权重是条件Akaike权重，不是模型贝叶斯后验。

`fit_over_trees(run, trees, models=c('DEC','DEC+J'),
likelihood_conditioning='tip_survival')`接受Newick集合、Nexus文件或phylo列表。
每棵树使用完整模型集。单个模型重试后仍失败，整棵树隔离；超过2%默认抛出错误。
保守返回模式是 `on_excess_failure='conservative_bounds'`。返回完整案例概率均值，
同时给出原目标树集上的概率缺失上下界。分支根部没有入枝，相关概率保持NA。
`monophyletic`仅匹配相同后代集合；`mrca`还报告外溢末端数。

检查点存储完整清单、节点数组和事件记录；改动输入或软件实现会拒绝混用缓存。
R目前外层串行、内层使用 `num_threads` 预算。该实现避免嵌套超订，不承诺外层加速。

## BSM精度与停止

```r
maps <- runBSM(dec, NULL, bsm_engine='native_whole_tree', seedval=42L,
               target_rmcse=.05, min_maps=200L, batch_size=200L, max_maps=5000L)
diagnostics <- bsm_diagnostics(maps, tree_target='fixed_sample')
diagnostics$stop_reason
diagnostics$route_rank_comparisons_df
```

固定树样本是默认目标，MCSE只衡量该给定树样本内的映射误差。
树总体目标需显式给 `tree_target='posterior_population'` 和
`tree_dependence='independent'` 或 `'mcmc_chain'`。MCMC按每个统计量计算ESS。
路线差区间整体落在等效带内才为near_tie，整体位于等效带一侧才为resolved，
其余为unresolved；罕见事件标为not_assessed。区间是逐比较探索性诊断，
不能声称同时覆盖、排名选择后覆盖或自适应停止后的95%覆盖保证。

当前新增概率契约与不确定性工具支持标准六模型、空间/环境/面积修饰符（`+x`、`+n`、`+w`、`+u`、
`areas_adjacency`），以及多层地质时间切片（`timeperiods`）与区域从古至今单调涌现（`areas_allowed`）；
非单调硬清零沉没掩码、化石钩、固定内部观测、性状耦合等组合仍明确拒绝并保留在 `legacy_unconditioned` 下运行。
