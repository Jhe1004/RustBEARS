# Maintenance and supported scope

RustBEARS v1.2 is currently a development release. The project plan commits to
at least five years of active maintenance after the first stable v1.2 release.
That period begins on the recorded release date, not the date of this draft.
Jian He is the RustBEARS maintainer. The project preserves attribution and
co-citation of the original BioGeoBEARS framework by Nicholas J. Matzke.

Maintenance covers installation, compatibility with supported R and Python
versions, reproducible numerical defects, and documentation. Issue reports
should include the release or commit, operating system, precision mode,
probability contract, complete model settings, and a minimal example. Each
release will publish its actual operating-system test results. The current
three-platform CI configuration still requires a successful remote run.

The first stable v1.2 scope is single-layer rooted binary trees with a fixed
range space: DEC, DIVALIKE, BAYAREALIKE, their +J alternatives, and validated
distance/manual-multiplier modifiers. The new uncertainty APIs reject
unvalidated time-stratified, direct-ancestor fossil, and trait-coupled models.
These restrictions concern the new statistical APIs; the historical interface
retains its separate compatibility contract.

Likelihood profiles, bootstrap tests, model averaging, external-tree summaries,
and mapping diagnostics have different inferential targets. The tutorials and
returned metadata record those targets. Fixed-tree Monte Carlo precision,
variation across external trees, and parameter uncertainty remain separate.
An empirical stopping rule does not provide an anytime-valid confidence band.

## 维护范围

项目承诺在首个稳定 v1.2 版本发布后至少主动维护五年，起算日以实际发布记录为准。
维护内容包括安装与环境适配、可复现的数值问题和文档修订。问题报告应附版本或提交号、
操作系统、精度模式、概率契约、完整模型设置和最小复现输入。每次发行保留实际的跨平台
检查结果；当前 CI 配置尚不能证明三个操作系统均已通过。

新统计模块首版支持单层、有根二叉树和固定状态空间。时间分层、直接祖先化石和性状
耦合扩展须经各自验证后才解锁。维护承诺不改变这些科学适用边界，也不将条件推断解释
为分子与地理数据的联合贝叶斯后验。
