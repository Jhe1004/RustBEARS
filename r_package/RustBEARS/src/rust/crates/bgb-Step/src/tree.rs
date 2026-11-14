mod stratified;
use bgb_core::{unpack_states, PrecisionMode};
use bgb_numerics::{ActionPlan, CompiledQTemplate, CsrMatrix, Matrix};
use rayon::prelude::*;
use std::collections::VecDeque;
use std::sync::{Arc, Mutex};
pub use stratified::CompiledStratifiedContext;

fn sum(v: &[f64]) -> f64 {
    #[cfg(feature = "r-blas")]
    unsafe {
        unsafe extern "C" {
            fn bgb_sum(x: *const f64, n: usize) -> f64;
        }
        bgb_sum(v.as_ptr(), v.len())
    }
    #[cfg(not(feature = "r-blas"))]
    {
        v.iter().sum()
    }
}
fn prod4(a: f64, b: f64, c: f64, d: f64) -> f64 {
    #[cfg(feature = "r-blas")]
    unsafe {
        unsafe extern "C" {
            fn bgb_prod4(a: f64, b: f64, c: f64, d: f64) -> f64;
        }
        bgb_prod4(a, b, c, d)
    }
    #[cfg(not(feature = "r-blas"))]
    {
        a * b * c * d
    }
}
fn missing() -> f64 {
    f64::from_bits(0x7ff00000000007a2)
}
fn normalized(v: &[f64]) -> Vec<f64> {
    let s = sum(v);
    v.iter().map(|x| x / s).collect()
}
fn transpose(q: &Matrix) -> Matrix {
    q.transpose()
}
fn mv(q: &Matrix, v: &[f64], forward: bool) -> Vec<f64> {
    assert_eq!(v.len(), q.dim(), "propagator vector dimension mismatch");
    #[cfg(feature = "r-blas")]
    {
        let mut out = vec![0.0; q.dim()];
        unsafe {
            unsafe extern "C" {
                fn bgb_dgemv(n: i32, a: *const f64, v: *const f64, out: *mut f64, trans: i32);
            }
            bgb_dgemv(
                q.dim() as i32,
                q.as_slice().as_ptr(),
                v.as_ptr(),
                out.as_mut_ptr(),
                forward as i32,
            );
        }
        out
    }
    #[cfg(not(feature = "r-blas"))]
    {
        if forward {
            transpose(q).mv(v)
        } else {
            q.mv(v)
        }
    }
}
struct Propagator {
    qt: Option<Matrix>,
    action: Option<ActionPlan>,
    forward: Option<ActionPlan>,
    legacy_dense: bool,
    fused: bool,
    cache: Mutex<VecDeque<(u64, Arc<Matrix>)>>,
}
impl Propagator {
    fn new(q: Matrix, mode: PrecisionMode, sparse: bool, fused: bool) -> Self {
        let qt = transpose(&q);
        let legacy_dense = mode == PrecisionMode::LegacyF32 && !sparse;
        let (action, forward) = if legacy_dense {
            (None, None)
        } else {
            (Some(ActionPlan::new(&q)), Some(ActionPlan::new(&qt)))
        };
        Self {
            qt: if legacy_dense { Some(qt) } else { None },
            action,
            forward,
            legacy_dense,
            fused,
            cache: Mutex::new(VecDeque::new()),
        }
    }
    fn from_csr(q: CsrMatrix, fused: bool) -> Result<Self, String> {
        let forward = ActionPlan::from_csr(q.transpose())?;
        let action = ActionPlan::from_csr(q)?;
        Ok(Self {
            qt: None,
            action: Some(action),
            forward: Some(forward),
            legacy_dense: false,
            fused,
            cache: Mutex::new(VecDeque::new()),
        })
    }
    fn apply(&self, t: f64, v: &[f64], forward: bool) -> Result<Vec<f64>, String> {
        if !self.legacy_dense {
            return if forward {
                self.forward.as_ref().unwrap().apply(t, v)
            } else {
                self.action.as_ref().unwrap().apply(t, v)
            };
        }
        let key = t.to_bits();
        let cached = self
            .cache
            .lock()
            .map_err(|_| "propagator cache poisoned")?
            .iter()
            .find(|x| x.0 == key)
            .map(|x| x.1.clone());
        let p = if let Some(p) = cached {
            p
        } else {
            // Match the original R transpose-before-Pade / row-major return.
            let p = Arc::new(transpose(
                &self.qt.as_ref().unwrap().legacy_exp(t, self.fused)?,
            ));
            let mut cache = self.cache.lock().map_err(|_| "propagator cache poisoned")?;
            if cache.len() >= 8 {
                cache.pop_front();
            }
            cache.push_back((key, p.clone()));
            p
        };
        Ok(mv(&p, v, forward))
    }
}
#[derive(Clone)]
struct Step {
    anc: usize,
    left: usize,
    right: usize,
    le: usize,
    re: usize,
}
pub struct CompiledModelContext {
    n: usize,
    tips: usize,
    nodes: usize,
    root: usize,
    steps: Vec<Step>,
    levels: Vec<Vec<usize>>,
    lengths: Vec<f64>,
    edge_of: Vec<Option<usize>>,
    parent: Vec<Option<usize>>,
    sister: Vec<Option<usize>>,
    min: f64,
    null: usize,
    root_edge: Option<f64>,
    kind: u32,
    method: u32,
    trait_len: usize,
    fixed: Vec<Vec<f64>>,
    mask: Vec<bool>,
    prior: Vec<f64>,
    tip_likes: Vec<f64>,
    down: Vec<f64>,
    cond: Vec<f64>,
    bottom: Vec<f64>,
    up: Vec<f64>,
    upbot: Vec<f64>,
    totals: Vec<f64>,
    pool: Arc<rayon::ThreadPool>,
    q_template: Option<CompiledQTemplate>,
    clado_template: Option<bgb_clado::CompiledCladoTemplate>,
}
fn index(x: f64, limit: usize) -> Result<usize, String> {
    if !x.is_finite() || x < 0.0 || x.fract() != 0.0 || x >= limit as f64 {
        Err("invalid tree/state index".into())
    } else {
        Ok(x as usize)
    }
}
impl CompiledModelContext {
    /// Static buffers: settings, edges (two columns, zero based), lengths,
    /// column-major tips, fixed nodes/likes, allowed states, root prior.
    pub fn new(a: &[&[f64]], threads: usize) -> Result<Self, String> {
        Self::new_with_pool(a, threads, None)
    }
    fn new_with_pool(
        a: &[&[f64]],
        threads: usize,
        shared_pool: Option<Arc<rayon::ThreadPool>>,
    ) -> Result<Self, String> {
        if a.len() != 13 || a[0].len() != 10 {
            return Err("tree context input count".into());
        }
        let c = a[0];
        let tips = c[0] as usize;
        let n = c[1] as usize;
        let nodes = c[2] as usize;
        if tips < 2
            || n == 0
            || nodes != 2 * tips - 1
            || a[1].len() != 2 * (nodes - 1)
            || a[2].len() != nodes - 1
            || (a[3].len() != tips * n && !(shared_pool.is_some() && a[3].is_empty()))
        {
            return Err("invalid binary tree/context dimensions".into());
        }
        if a[2].iter().any(|x| !x.is_finite()) || !c[3].is_finite() {
            return Err("invalid branch length".into());
        }
        let ne = nodes - 1;
        let mut steps = Vec::new();
        let mut depths = vec![0; nodes];
        let mut seen = vec![false; nodes];
        seen[..tips].fill(true);
        let mut levels: Vec<Vec<usize>> = Vec::new();
        let mut edge_of = vec![None; nodes];
        let mut parent = vec![None; nodes];
        let mut sister = vec![None; nodes];
        for e in (0..ne).step_by(2) {
            let anc = index(a[1][e], nodes)?;
            if index(a[1][e + 1], nodes)? != anc {
                return Err("tree is not binary pruningwise".into());
            }
            let l = index(a[1][e + ne], nodes)?;
            let r = index(a[1][e + 1 + ne], nodes)?;
            if anc < tips
                || seen[anc]
                || !seen[l]
                || !seen[r]
                || l == r
                || edge_of[l].is_some()
                || edge_of[r].is_some()
            {
                return Err("invalid pruning order".into());
            }
            let depth = depths[l].max(depths[r]) + 1;
            depths[anc] = depth;
            seen[anc] = true;
            while levels.len() < depth {
                levels.push(Vec::new());
            }
            levels[depth - 1].push(steps.len());
            edge_of[l] = Some(e);
            edge_of[r] = Some(e + 1);
            parent[l] = Some(anc);
            parent[r] = Some(anc);
            sister[l] = Some(r);
            sister[r] = Some(l);
            steps.push(Step {
                anc,
                left: l,
                right: r,
                le: e,
                re: e + 1,
            });
        }
        let root = steps.last().ok_or("empty tree")?.anc;
        if seen.iter().any(|x| !*x) || edge_of[root].is_some() {
            return Err("disconnected tree".into());
        }
        let mut fixed = vec![Vec::new(); nodes];
        let nf = a[4].len();
        if a[5].len() != nf * n {
            return Err("fixed likelihood dimensions".into());
        }
        for k in 0..nf {
            let node = index(a[4][k], nodes)?;
            fixed[node] = (0..n).map(|j| a[5][k + nf * j]).collect();
        }
        if !a[6].is_empty() && a[6].len() != n || !a[7].is_empty() && a[7].len() != n {
            return Err("root prior or allowed state dimensions".into());
        }
        let mask = if a[6].is_empty() {
            vec![true; n]
        } else {
            a[6].iter().map(|&x| x != 0.0).collect()
        };
        let pool = match shared_pool {
            Some(pool) => pool,
            None => Arc::new(
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads.max(1))
                    .build()
                    .map_err(|e| e.to_string())?,
            ),
        };
        if a[12].len() != 2 {
            return Err("model template dimensions".into());
        }
        let q_template = if a[9].is_empty() {
            None
        } else {
            Some(CompiledQTemplate::new(
                &unpack_states(a[8], a[9])?,
                a[12][0] as usize,
            )?)
        };
        let clado_template = if a[11].is_empty() {
            None
        } else {
            Some(compile_clado(a[10], a[11], a[12][1] as usize)?)
        };
        Ok(Self {
            n,
            tips,
            nodes,
            root,
            steps,
            levels,
            lengths: a[2].to_vec(),
            edge_of,
            parent,
            sister,
            min: c[3],
            null: (c[4] != 0.0) as usize,
            root_edge: if c[5] != 0.0 && c[6] > 0.0 {
                Some(c[6])
            } else {
                None
            },
            kind: c[7] as u32,
            method: c[8] as u32,
            trait_len: c[9] as usize,
            fixed,
            mask,
            prior: a[7].to_vec(),
            tip_likes: if a[3].is_empty() {
                vec![0.0; tips * n]
            } else {
                a[3].to_vec()
            },
            down: vec![0.0; nodes * n],
            cond: vec![0.0; nodes * n],
            bottom: Vec::new(),
            up: Vec::new(),
            upbot: Vec::new(),
            totals: vec![1.0; nodes],
            pool,
            q_template,
            clado_template,
        })
    }
    fn row<'a>(&self, x: &'a [f64], node: usize) -> &'a [f64] {
        &x[node * self.n..(node + 1) * self.n]
    }
    fn mask(&self, x: &mut [f64]) {
        for (a, &yes) in x.iter_mut().zip(&self.mask) {
            if !yes {
                *a = 0.0;
            }
        }
    }
    fn hook(&self, s: &Step) -> bool {
        self.lengths[s.le] < self.min || self.lengths[s.re] < self.min
    }
    /// Dynamic buffers: settings[COO,sparse,ancprobs,full], Q, rows, cols,
    /// event ancestor/left/right/weight, row denominators, explicit speciation matrix.
    pub fn evaluate(
        &mut self,
        a: &[&[f64]],
        mode: PrecisionMode,
        fused: bool,
    ) -> Result<Vec<Vec<f64>>, String> {
        if a.len() == 23 {
            let [power] = a[22] else {
                return Err("invalid branch length exponent".into());
            };
            if !power.is_finite() || *power <= 0.0 {
                return Err("invalid branch length exponent".into());
            }
            let scaled = self.lengths.iter().map(|x| x.powf(*power)).collect();
            let original = std::mem::replace(&mut self.lengths, scaled);
            let result = self.evaluate(&a[..22], mode, fused);
            self.lengths = original;
            return result;
        }
        if a.len() != 22 || a[0].len() != 4 {
            return Err("tree evaluation input count".into());
        }
        let c = a[0];
        let ancprobs = c[2] != 0.0;
        let full = c[3] != 0.0;
        let prop = build_propagator(self.n, self.q_template.as_ref(), a, mode, fused)?;
        let events = build_clado(
            self.n,
            self.null,
            self.kind,
            self.method,
            self.clado_template.as_ref(),
            a,
            mode,
        )?;
        let (rootcond, rootbottom) = self.downpass(&prop, &events, c[1] != 0.0, ancprobs)?;
        if self.totals.iter().any(|x| *x < 0.0) {
            return Err("negative likelihood requires original R warning semantics".into());
        }
        let logs: Vec<_> = self.totals.iter().map(|x| x.ln()).collect();
        let ln_l = sum(&logs);
        if ancprobs {
            self.uppass(&prop, &events, None)?;
        }
        if !full {
            return Ok(vec![vec![ln_l]]);
        }
        let nr = self.nodes + self.root_edge.is_some() as usize;
        let mut down = column_major(&self.down, self.nodes, self.n);
        let mut cond = column_major(&self.cond, self.nodes, self.n);
        if self.root_edge.is_some() {
            down = append_row(&down, self.nodes, self.n, &normalized(&rootcond));
            cond = append_row(&cond, self.nodes, self.n, &rootcond);
        }
        let mut out = vec![
            vec![ln_l, nr as f64],
            self.totals.clone(),
            down,
            cond,
            rootbottom,
        ];
        if ancprobs {
            out.push(column_major(&self.bottom, self.nodes, self.n));
            out.push(column_major(&self.upbot, self.nodes, self.n));
            out.push(column_major(&self.up, self.nodes, self.n));
        }
        Ok(out)
    }
    fn downpass(
        &mut self,
        prop: &Propagator,
        events: &Clado,
        sparse: bool,
        ancprobs: bool,
    ) -> Result<(Vec<f64>, Vec<f64>), String> {
        self.down.fill(0.0);
        self.cond.fill(0.0);
        self.totals.resize(self.nodes, 1.0);
        self.totals.fill(1.0);
        for i in 0..self.tips {
            for j in 0..self.n {
                let x = self.tip_likes[i + self.tips * j];
                self.down[i * self.n + j] = x;
                self.cond[i * self.n + j] = x;
            }
            let s = sum(&self.down[i * self.n..(i + 1) * self.n]);
            self.totals[i] = s / s;
        }
        if ancprobs {
            self.bottom.resize(self.nodes * self.n, missing());
            self.bottom.fill(missing());
        }
        for level in &self.levels {
            let records = self.pool.install(|| {
                level
                    .par_iter()
                    .map(|&k| {
                        let st = &self.steps[k];
                        let mut l =
                            prop.apply(self.lengths[st.le], self.row(&self.down, st.left), false)?;
                        let mut r =
                            prop.apply(self.lengths[st.re], self.row(&self.down, st.right), false)?;
                        if sparse {
                            for x in l.iter_mut().chain(r.iter_mut()) {
                                if *x < 0.0 {
                                    *x = 0.0;
                                }
                            }
                        }
                        self.mask(&mut l);
                        self.mask(&mut r);
                        let hook = self.hook(st);
                        let mut v = if hook || self.kind == 0 {
                            l.iter().zip(&r).map(|(a, b)| a * b).collect::<Vec<_>>()
                        } else {
                            events.split(&l, &r)?
                        };
                        if !hook && self.method == 3 && self.trait_len > 0 && self.kind == 1 {
                            let count = v[self.null..].iter().filter(|x| x.is_nan()).count();
                            if count == self.trait_len - 1 {
                                for x in &mut v[self.null..] {
                                    if x.is_nan() {
                                        *x = 0.0;
                                    }
                                }
                            }
                        }
                        if !hook && self.kind != 0 && !self.fixed[st.anc].is_empty() {
                            for (x, f) in v.iter_mut().zip(&self.fixed[st.anc]) {
                                *x *= f;
                            }
                        }
                        if !hook && self.kind == 1 {
                            self.mask(&mut v);
                        }
                        if st.anc == self.root && !self.prior.is_empty() {
                            for (x, f) in v.iter_mut().zip(&self.prior) {
                                *x *= f;
                            }
                        }
                        if v.iter().any(|x| x.is_nan()) {
                            return Err(format!("node_nan:{}", st.anc + 1));
                        }
                        Ok((
                            k,
                            v,
                            if ancprobs { normalized(&l) } else { Vec::new() },
                            if ancprobs { normalized(&r) } else { Vec::new() },
                        ))
                    })
                    .collect::<Result<Vec<_>, String>>()
            })?;
            // Commit in original step order; no parallel floating-point reduction.
            for (k, v, l, r) in records {
                let st = &self.steps[k];
                let total = sum(&v);
                self.totals[st.anc] = total;
                self.cond[st.anc * self.n..(st.anc + 1) * self.n].copy_from_slice(&v);
                for j in 0..self.n {
                    self.down[st.anc * self.n + j] = v[j] / total;
                }
                if ancprobs {
                    self.bottom[st.left * self.n..(st.left + 1) * self.n].copy_from_slice(&l);
                    self.bottom[st.right * self.n..(st.right + 1) * self.n].copy_from_slice(&r);
                }
            }
        }
        let mut rootbottom = self.row(&self.down, self.root).to_vec();
        let mut rootcond = Vec::new();
        if let Some(t) = self.root_edge {
            rootcond = prop.apply(t, &rootbottom, false)?;
            let total = sum(&rootcond);
            self.totals.push(total);
            rootbottom = rootcond.iter().map(|x| x / total).collect();
        }
        if self.mask.iter().any(|x| !*x) {
            self.mask(&mut rootbottom);
            rootbottom = normalized(&rootbottom);
        }
        Ok((rootcond, rootbottom))
    }
    fn uppass(
        &mut self,
        p: &Propagator,
        c: &Clado,
        root_seed: Option<(&[f64], &[f64])>,
    ) -> Result<(), String> {
        let fill = if root_seed.is_some() { 0.0 } else { missing() };
        self.up.resize(self.nodes * self.n, fill);
        self.up.fill(fill);
        self.upbot.resize(self.nodes * self.n, fill);
        self.upbot.fill(fill);
        let root = self.root * self.n..(self.root + 1) * self.n;
        if let Some((top, bottom)) = root_seed {
            self.up[root.clone()].copy_from_slice(top);
            self.upbot[root].copy_from_slice(bottom);
        } else {
            self.up[root].fill(1.0 / self.n as f64);
        }
        // Reverse wavefronts keep all mothers complete before their children.
        // Each record owns its arithmetic; commits never reduce floats across nodes.
        for level in self.levels.iter().rev() {
            let records = self.pool.install(|| {
                level
                    .par_iter()
                    .rev()
                    .map(|&k| {
                        let st = &self.steps[k];
                        let hook = self.hook(st);
                        let (mut top, mut bot) = if st.anc == self.root {
                            (
                                self.row(&self.up, st.anc).to_vec(),
                                self.row(&self.upbot, st.anc).to_vec(),
                            )
                        } else {
                            let mom = self.parent[st.anc].ok_or("missing mother")?;
                            let sister = self.sister[st.anc].ok_or("missing sister")?;
                            let sister_left = sister < st.anc;
                            let bot = c.uppass(
                                self.row(&self.up, mom),
                                if sister_left {
                                    Some(self.row(&self.bottom, sister))
                                } else {
                                    None
                                },
                                if sister_left {
                                    None
                                } else {
                                    Some(self.row(&self.bottom, sister))
                                },
                                hook,
                                !sister_left,
                            )?;
                            let top = p.apply(
                                self.lengths
                                    [self.edge_of[st.anc].ok_or("missing ancestor edge")?],
                                &bot,
                                true,
                            )?;
                            (top, bot)
                        };
                        let mut tips = Vec::with_capacity(2);
                        for (node, other, is_left, edge) in [
                            (st.left, st.right, true, st.le),
                            (st.right, st.left, false, st.re),
                        ] {
                            if node < self.tips {
                                let mut tipbot = c.uppass(
                                    &top,
                                    if is_left {
                                        None
                                    } else {
                                        Some(self.row(&self.bottom, other))
                                    },
                                    if is_left {
                                        Some(self.row(&self.bottom, other))
                                    } else {
                                        None
                                    },
                                    hook,
                                    is_left,
                                )?;
                                let mut tiptop = p.apply(self.lengths[edge], &tipbot, true)?;
                                self.mask(&mut tipbot);
                                self.mask(&mut tiptop);
                                tips.push((node, tiptop, tipbot));
                            }
                        }
                        // Preserve R's timing: tips use the ancestor before its mask/fixlikes.
                        self.mask(&mut top);
                        self.mask(&mut bot);
                        for (j, f) in self.fixed[st.anc].iter().enumerate() {
                            top[j] *= f;
                        }
                        Ok((k, top, bot, tips))
                    })
                    .collect::<Result<Vec<_>, String>>()
            })?;
            for (k, top, bot, tips) in records {
                let st = &self.steps[k];
                self.up[st.anc * self.n..(st.anc + 1) * self.n].copy_from_slice(&top);
                self.upbot[st.anc * self.n..(st.anc + 1) * self.n].copy_from_slice(&bot);
                for (node, top, bot) in tips {
                    self.up[node * self.n..(node + 1) * self.n].copy_from_slice(&top);
                    self.upbot[node * self.n..(node + 1) * self.n].copy_from_slice(&bot);
                }
                for node in [st.left, st.right] {
                    for j in 0..self.n {
                        if !self.mask[j] {
                            self.up[node * self.n + j] = 0.0;
                            self.upbot[node * self.n + j] = 0.0;
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

fn compile_clado(
    flat: &[f64],
    lens: &[f64],
    na: usize,
) -> Result<bgb_clado::CompiledCladoTemplate, String> {
    let states = unpack_states(flat, lens)?;
    let ones = vec![1.0; na * na];
    let params = [1.0; 4];
    let max = vec![na as f64; na];
    let model = bgb_clado::Model {
        states: &states,
        parameters: &params,
        d: &ones,
        priors: [&ones, &ones, &ones, &ones],
        maxmin: &max,
        na,
        mode: PrecisionMode::StrictF64,
    };
    bgb_clado::CompiledCladoTemplate::new(&model)
}
fn build_propagator(
    n: usize,
    template: Option<&CompiledQTemplate>,
    a: &[&[f64]],
    mode: PrecisionMode,
    fused: bool,
) -> Result<Propagator, String> {
    let c = a[0];
    if a[10].len() != 4 {
        return Err("model parameter settings".into());
    }
    let prop = if a[10][0] != 0.0 {
        let tpl = template.ok_or("missing compiled Q template")?;
        if a[10][2] != 0.0 {
            let q = tpl.refresh_csr(a[11], a[12], a[13], a[10][1] != 0.0, a[10][3], mode)?;
            Propagator::from_csr(q, fused)?
        } else {
            let q = tpl.refresh(a[11], a[12], a[13], a[10][1] != 0.0, false, a[10][3], mode)?;
            Propagator::new(Matrix::new(n, q[0].clone())?, mode, c[1] != 0.0, fused)
        }
    } else if c[0] != 0.0 {
        let q = CsrMatrix::from_triplets(n, a[2], a[3], a[1], true)?;
        if c[1] != 0.0 {
            Propagator::from_csr(q, fused)?
        } else {
            Propagator::new(q.dense(), mode, false, fused)
        }
    } else {
        Propagator::new(Matrix::new(n, a[1].to_vec())?, mode, c[1] != 0.0, fused)
    };
    Ok(prop)
}
fn build_clado(
    n: usize,
    null: usize,
    kind: u32,
    method: u32,
    template: Option<&bgb_clado::CompiledCladoTemplate>,
    a: &[&[f64]],
    mode: PrecisionMode,
) -> Result<Clado, String> {
    let compiled_columns;
    let compiled_rows;
    let (ea, el, er, ew, den) = if !a[14].is_empty() {
        let tpl = template.ok_or("missing compiled cladogenesis template")?;
        if a[21].len() != 1 {
            return Err("cladogenesis area count".into());
        }
        let model = bgb_clado::Model {
            states: &tpl.states,
            parameters: a[14],
            d: a[15],
            priors: [a[16], a[17], a[18], a[19]],
            maxmin: a[20],
            na: a[21][0] as usize,
            mode,
        };
        let refreshed = tpl.refresh(&model)?;
        let actual = refreshed
            .iter()
            .flat_map(|x| [x.0, x.1, x.2])
            .max()
            .map_or(0, |x| x + 1);
        if actual != n - null {
            return Err("nonstandard cladogenesis row count".into());
        }
        let mut columns = vec![Vec::with_capacity(refreshed.len()); 4];
        for (anc, l, r, w) in refreshed {
            columns[0].push(anc as f64);
            columns[1].push(l as f64);
            columns[2].push(r as f64);
            columns[3].push(w);
        }
        compiled_rows = bgb_clado::rowsums(&columns[0], &columns[3], actual)?;
        compiled_columns = columns;
        (
            &compiled_columns[0][..],
            &compiled_columns[1][..],
            &compiled_columns[2][..],
            &compiled_columns[3][..],
            &compiled_rows[..],
        )
    } else {
        (a[4], a[5], a[6], a[7], a[8])
    };
    Clado::new(n, null, ea, el, er, ew, den, a[9], kind, method, mode)
}
fn column_major(x: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    let mut v = vec![0.0; rows * cols];
    for i in 0..rows {
        for j in 0..cols {
            v[i + rows * j] = x[i * cols + j];
        }
    }
    v
}
fn append_row(x: &[f64], rows: usize, cols: usize, row: &[f64]) -> Vec<f64> {
    let mut v = vec![0.0; (rows + 1) * cols];
    for j in 0..cols {
        v[j * (rows + 1)..j * (rows + 1) + rows].copy_from_slice(&x[j * rows..(j + 1) * rows]);
        v[j * (rows + 1) + rows] = row[j];
    }
    v
}

struct Clado {
    n: usize,
    offset: usize,
    events: Vec<(usize, usize, usize, f64)>,
    up_rows: Vec<f64>,
    hook_states: usize,
    by_left: Vec<Vec<usize>>,
    by_right: Vec<Vec<usize>>,
    denom: Vec<f64>,
    matrix: Vec<f64>,
    kind: u32,
    method: u32,
    mode: PrecisionMode,
}
impl Clado {
    fn new(
        n: usize,
        offset: usize,
        a: &[f64],
        l: &[f64],
        r: &[f64],
        w: &[f64],
        denom: &[f64],
        matrix: &[f64],
        kind: u32,
        method: u32,
        mode: PrecisionMode,
    ) -> Result<Self, String> {
        if a.len() != l.len() || a.len() != r.len() || a.len() != w.len() {
            return Err("cladogenesis COO dimensions".into());
        }
        let m = n - offset;
        let mut events = Vec::with_capacity(a.len());
        for k in 0..a.len() {
            events.push((index(a[k], m)?, index(l[k], m)?, index(r[k], m)?, w[k]));
        }
        if kind == 1 && denom.len() != m {
            return Err("nonstandard cladogenesis row count".into());
        }
        if kind == 3 && matrix.len() != m * m * m {
            return Err("explicit speciation dimensions".into());
        }
        let mut up_rows = vec![0.0; n];
        let mut by_left = vec![Vec::new(); m];
        let mut by_right = vec![Vec::new(); m];
        for (k, &(a, i, j, w)) in events.iter().enumerate() {
            up_rows[a] += w;
            by_left[i].push(k);
            by_right[j].push(k);
        }
        Ok(Self {
            n,
            offset,
            events,
            up_rows,
            hook_states: m,
            by_left,
            by_right,
            denom: denom.to_vec(),
            matrix: matrix.to_vec(),
            kind,
            method,
            mode,
        })
    }
    fn split(&self, l: &[f64], r: &[f64]) -> Result<Vec<f64>, String> {
        let m = self.n - self.offset;
        let mut out = vec![0.0; self.n];
        if self.kind == 3 {
            for a in 0..m {
                let terms: Vec<_> = (0..m * m)
                    .map(|k| {
                        l[self.offset + k / m] * r[self.offset + k % m] * self.matrix[a + m * k]
                    })
                    .collect();
                out[a + self.offset] = sum(&terms);
            }
            return Ok(out);
        }
        for &(a, i, j, w) in &self.events {
            if self.method == 2 {
                let cell = self.mode.round(
                    l[i + self.offset] * r[j + self.offset] * self.mode.round(w) / self.denom[a],
                );
                out[a + self.offset] = self.mode.round(out[a + self.offset] + cell);
            } else {
                out[a + self.offset] +=
                    self.mode.round(l[i + self.offset] * r[j + self.offset] * w);
            }
        }
        if self.method != 2 {
            for i in 0..m {
                out[i + self.offset] = self.mode.round(out[i + self.offset] / self.denom[i]);
            }
        }
        Ok(out)
    }
    fn uppass(
        &self,
        anc: &[f64],
        left: Option<&[f64]>,
        right: Option<&[f64]>,
        hook: bool,
        want_left: bool,
    ) -> Result<Vec<f64>, String> {
        let count = if hook {
            self.hook_states
        } else {
            self.events.len()
        };
        let mut probs = Vec::with_capacity(count);
        for k in 0..count {
            let (a, i, j, w) = if hook { (k, k, k, 1.0) } else { self.events[k] };
            let v = prod4(
                anc[a + self.offset],
                left.map_or(1.0, |l| l[i + self.offset]),
                right.map_or(1.0, |r| r[j + self.offset]),
                w / if hook { 1.0 } else { self.up_rows[a] },
            );
            probs.push(v);
        }
        let total = sum(&probs);
        if total.is_nan() {
            return Err("STOP ERROR in calc_uppass_scenario_probs_new2(): NAs found in sum_relprob_each_split_scenario".into());
        }
        for x in &mut probs {
            *x = if total > 0.0 { *x / total } else { 0.0 };
        }
        let total2 = sum(&probs);
        for x in &mut probs {
            *x /= total2;
        }
        let mut out = vec![0.0; self.n];
        // R sums each matching scenario subsequence in extended precision.
        for state in 0..self.n - self.offset {
            out[state + self.offset] = if hook {
                probs.get(state).copied().unwrap_or(0.0)
            } else {
                let indices = if want_left {
                    &self.by_left[state]
                } else {
                    &self.by_right[state]
                };
                #[cfg(feature = "r-blas")]
                unsafe {
                    unsafe extern "C" {
                        fn bgb_indexed_sum(x: *const f64, indices: *const usize, n: usize) -> f64;
                    }
                    bgb_indexed_sum(probs.as_ptr(), indices.as_ptr(), indices.len())
                }
                #[cfg(not(feature = "r-blas"))]
                {
                    indices.iter().map(|&k| probs[k]).sum()
                }
            };
        }
        Ok(out)
    }
}
