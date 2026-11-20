# Statistical Formulations, Precision Modes, and Scope

<p align="left">
  <a href="../README.md">Back to English README</a> | <a href="../README_zh-CN.md">返回中文主页</a>
</p>

## 1. Three Likelihood Formulations on an $N$-Tip Phylogeny

Over $r$ discrete areas with state space $\mathcal{S} = \{\emptyset\} \cup \mathcal{S}_+$ of size $S = |\mathcal{S}|$, single-area states contract into the absorbing null range $\emptyset$ at rate $e \ge 0$. `RustBEARS` explicitly distinguishes three likelihood formulations:

1. **`"legacy_unconditioned"` (Compatibility Default)**
   - Sets cladogenetic transitions from $\emptyset$ to zero and combines root conditional likelihoods $\mathbf{L}_{\text{root}}$ with historical root weights without conditioning on tip survival.
   - Retains exact parity with upstream `BioGeoBEARS` (disabling $\emptyset$ recovers `DEC*` of Massana et al. 2015).
   - Supports the full historical feature set, including fossil tips, non-monotone area submergence masks, detection models, and trait-coupled models.

2. **`"tip_survival"` (Recommended for Contemporary Extant Clades)**
   - Retains single-area extirpation into $\emptyset$ along branches while conditioning observed tip ranges $Y = (Y_1, \dots, Y_N)$ on every sampled extant tip occupying a non-empty present-day range:
     $$\ln L_{\text{surv}}(\theta) = \ln\!\left(\sum_{i \in \mathcal{S}_+^{(M)}} \pi_i^{(M)} L_{\text{root}, i}\right) - \ln P_{\text{surv}}(\theta),$$
     where $P_{\text{surv}}(\theta) = \sum_{i \in \mathcal{S}_+^{(M)}} \pi_i^{(M)} \tilde{L}_{\text{root}, i}$ is evaluated by a second pruning pass initialized with ones across allowed present-day non-empty states and zero at $\emptyset$.
   - Eliminates the unconditioned survivorship penalty that drives DEC extinction estimates toward zero and inflates false selection of `+J`.

3. **`"unconditional_complete"` (Generative Verification Mode)**
   - Sets $C(\emptyset \to \emptyset, \emptyset) = 1$ under a normalized root prior $\pi$ on $\mathcal{S}_+$ (requires `allow_null_tips = True` / `TRUE`), summing strictly to 1 across all $S^N$ tip outcomes (including extinct tips). Under uniform weights over the $K$ allowed non-empty root states, its log-likelihood on extant non-empty tips differs from unit-weight legacy pruning by $-\ln K$.

---

## 2. Chronologically Monotone Area Emergence in Stratified Models

`tip_survival` and `unconditional_complete` extend across multi-stratum models ($0 = \tau_0 < \tau_1 < \dots < \tau_M$) when areas emerge progressively toward the present ($\mathcal{A}_M \subseteq \mathcal{A}_{M-1} \subseteq \dots \subseteq \mathcal{A}_1$). Under monotone emergence:
- Every range valid in older stratum $m$ remains valid entering younger stratum $m-1$ ($\mathcal{S}_+^{(m)} \subseteq \mathcal{S}_+^{(m-1)}$), conserving forward probability mass across stratum boundaries without ad hoc renormalization;
- Root priors $\pi^{(M)}$ are normalized strictly over non-empty states allowed in the oldest stratum;
- Non-monotone masks in which previously allowed areas vanish abruptly ($\mathcal{A}_m \not\subseteq \mathcal{A}_{m-1}$) without a physical extirpation transition are rejected outside `"legacy_unconditioned"`.

---

## 3. Numerical Precision Modes

- **`StrictF64` / `"strict_f64"` (Default)**: Pure 64-bit IEEE-754 double precision across rate-matrix construction, cladogenetic weight normalization, dense Padé and sparse Krylov matrix exponentials, and log-space pruning.
- **`LegacyF32` / `"legacy_f32"` (Optional Historical Audit Mode)**: Emulates the 32-bit single-precision rounding points of historical `cladoRcpp` and `EXPOKIT` helpers to reproduce published `BioGeoBEARS` outputs to machine precision (`median |ΔlnL| = 9.95e-14` across 1,339 configurations from 77 studies).
