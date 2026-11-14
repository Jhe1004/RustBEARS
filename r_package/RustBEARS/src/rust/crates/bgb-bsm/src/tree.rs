//! Whole-tree conditional BSM. GPL-2.0-or-later.
use crate::{reject, Event, Generator, Xoshiro256PlusPlus};
use bgb_numerics::{ActionPlan, Matrix};
use rayon::prelude::*;
use std::{cell::RefCell, sync::OnceLock};

pub struct LayerInput<'a> {
    pub q: &'a [f64],
    pub allowed: &'a [f64],
    pub events: &'a [f64],
    pub dispersal: &'a [f64],
}
#[derive(Clone)]
pub struct SegmentInput {
    pub node: usize,
    pub layer: usize,
    pub length: f64,
    pub start_age: f64,
}
pub struct TreeInput<'a> {
    pub states: &'a [Vec<i32>],
    pub tips: &'a [f64],
    pub tip_count: usize,
    pub root: usize,
    pub steps: &'a [[usize; 3]],
    pub ages: &'a [f64],
    pub layers_at_node: &'a [usize],
    pub hooks: &'a [f64],
    pub fixes: &'a [f64],
    pub segments: &'a [SegmentInput],
    pub layers: &'a [LayerInput<'a>],
}

#[derive(Clone, Copy)]
struct Split {
    left: usize,
    right: usize,
    prior: f64,
    kind: usize,
}
struct Layer {
    generator: Generator,
    q: Vec<f64>,
    backward: ActionPlan,
    forward: ActionPlan,
    splits: Vec<Vec<Split>>,
    allowed: Vec<f64>,
    dispersal: Vec<f64>,
}
struct Segment {
    spec: SegmentInput,
    target: Vec<f64>,
    rows: Vec<OnceLock<Result<Vec<f64>, String>>>,
}

fn positive(weights: &[f64], rng: &mut Xoshiro256PlusPlus) -> Result<usize, String> {
    let mut sample = rng.uniform() * weights.iter().sum::<f64>();
    let mut last = None;
    for (i, &w) in weights.iter().enumerate() {
        if w > 0.0 {
            last = Some(i);
            sample -= w;
            if sample <= 0.0 {
                return Ok(i);
            }
        }
    }
    last.ok_or_else(|| "zero conditional state mass".into())
}

fn normalize(x: &mut [f64]) -> Result<(), String> {
    let total: f64 = x.iter().sum();
    if total > 0.0 && total.is_finite() {
        x.iter_mut().for_each(|v| *v /= total);
        Ok(())
    } else {
        Err("zero or invalid branch conditional likelihood".into())
    }
}

fn areas(s: &[i32]) -> impl Iterator<Item = usize> + '_ {
    s.iter().filter(|&&x| x >= 0).map(|&x| x as usize)
}

fn split_type(parent: &[i32], left: &[i32], right: &[i32]) -> usize {
    if left == parent && right == parent {
        4 // y
    } else if left.iter().chain(right).any(|x| !parent.contains(x)) {
        3 // j
    } else if left == parent || right == parent {
        5 // s
    } else {
        6 // v
    }
}

fn select_source(
    state: &[i32],
    new_area: usize,
    dmat: &[f64],
    na: usize,
    rng: &mut Xoshiro256PlusPlus,
) -> f64 {
    let ids: Vec<_> = areas(state).collect();
    if ids.is_empty() {
        return -1.0;
    }
    let weights: Vec<_> = ids.iter().map(|&x| dmat[x + new_area * na]).collect();
    positive(&weights, rng)
        .map(|j| ids[j] as f64)
        .unwrap_or_else(|_| ids[(rng.uniform() * ids.len() as f64) as usize] as f64)
}

fn bridge(
    layer: &Layer,
    from: usize,
    to: usize,
    t: f64,
    width: usize,
    rejection: bool,
    rng: &mut Xoshiro256PlusPlus,
) -> Result<Vec<Event>, String> {
    if !rejection || t == 0.0 {
        return layer.generator.sample(from, to, t, rng);
    }
    let source = RefCell::new(rng);
    let trial = reject(
        &layer.q,
        width,
        from,
        to,
        t,
        40_000,
        |rate| -source.borrow_mut().uniform().ln() / rate,
        |weights, _, _| positive(weights, &mut source.borrow_mut()).unwrap_or(0),
    )?;
    if trial.accepted {
        Ok(trial.events)
    } else {
        Err("rejection bridge exhausted 40000 attempts; use uniformization".into())
    }
}

// Owned geometry: no borrowed caller arrays survive preparation.
struct OwnedTree {
    states: Vec<Vec<i32>>,
    tips: Vec<f64>,
    tip_count: usize,
    root: usize,
    steps: Vec<[usize; 3]>,
    ages: Vec<f64>,
    layers_at_node: Vec<usize>,
    hooks: Vec<f64>,
    fixes: Vec<f64>,
    segments: Vec<SegmentInput>,
}

/// Immutable prepared probability tree with thread-safe lazy transition rows.
/// Mode 0 observes tips, 1 simulates nonempty tips, 2 simulates the complete
/// absorbing process. Internal mode 3 preserves the old whole-tree sampler.
pub struct PreparedTreeSampler {
    input: OwnedTree,
    layers: Vec<Layer>,
    branches: Vec<Vec<Segment>>,
    down: Vec<Vec<f64>>,
    branch_likes: Vec<Vec<f64>>,
    na: usize,
    conditioning: u32,
}

/// Legacy adapter: preserve old state order, arithmetic and RNG consumption.
pub fn sample_tree(
    input: TreeInput<'_>,
    maps: usize,
    seed: u64,
    threads: usize,
    rejection: bool,
) -> Result<Vec<Vec<f64>>, String> {
    PreparedTreeSampler::prepare_impl(input, &[], 3)?.draw_bsm(1, maps, seed, threads, rejection)
}

impl PreparedTreeSampler {
    pub fn prepare(
        input: TreeInput<'_>,
        root_weights: &[f64],
        conditioning: u32,
    ) -> Result<Self, String> {
        if conditioning > 2 {
            return Err("invalid sampler conditioning".into());
        }
        Self::prepare_impl(input, root_weights, conditioning)
    }

    fn prepare_impl(
        source: TreeInput<'_>,
        root_weights: &[f64],
        conditioning: u32,
    ) -> Result<Self, String> {
        let input = OwnedTree {
            states: source.states.to_vec(),
            tips: source.tips.to_vec(),
            tip_count: source.tip_count,
            root: source.root,
            steps: source.steps.to_vec(),
            ages: source.ages.to_vec(),
            layers_at_node: source.layers_at_node.to_vec(),
            hooks: source.hooks.to_vec(),
            fixes: source.fixes.to_vec(),
            segments: source.segments.to_vec(),
        };
        let n = input.states.len();
        let nn = input.ages.len();
        let nt = input.tip_count;
        let na = input
            .states
            .iter()
            .flat_map(|s| areas(s))
            .max()
            .ok_or("no areas")?
            + 1;
        if n == 0
            || nn != 2 * nt - 1
            || nt < 2
            || input.root >= nn
            || input.tips.len() != n * nt
            || input.steps.len() != nn - nt
            || input.layers_at_node.len() != nn
            || input.hooks.len() != nn
            || input.fixes.len() != nn * n
            || input
                .layers_at_node
                .iter()
                .any(|&l| l >= source.layers.len())
        {
            return Err("invalid BSM tree dimensions".into());
        }
        if conditioning != 3 {
            if source.layers.len() != 1
                || input.hooks.iter().any(|x| *x != 0.0)
                || input.fixes.iter().any(|x| *x != 1.0)
                || source.layers[0].allowed.iter().any(|x| *x != 1.0)
            {
                return Err("UnsupportedModelContractError: prepared Tier-1a sampler requires one layer, no hooks or fixed-node constraints".into());
            }
            if root_weights.len() != n
                || root_weights.iter().any(|x| !x.is_finite() || *x < 0.0)
                || root_weights.iter().sum::<f64>() <= 0.0
            {
                return Err("invalid sampler root weights".into());
            }
            if conditioning != 0
                && ((root_weights.iter().sum::<f64>() - 1.0).abs() > 1e-12
                    || input
                        .states
                        .first()
                        .map_or(true, |x| !x.iter().all(|i| *i < 0))
                    || (conditioning == 1 && root_weights[0] != 0.0))
            {
                return Err(
                    "simulation requires normalized root prior and an explicit empty state".into(),
                );
            }
        }
        let layers: Vec<_> = source
            .layers
            .iter()
            .map(|l| {
                if l.allowed.len() != n || l.dispersal.len() != na * na || l.events.len() % 4 != 0 {
                    return Err("invalid BSM layer dimensions".into());
                }
                let generator = Generator::new(l.q, n)?;
                let q = Matrix::new(n, l.q.to_vec())?;
                let qt = Matrix::new(n, (0..n * n).map(|k| l.q[k / n + n * (k % n)]).collect())?;
                let mut splits = vec![Vec::new(); n];
                let mut totals = vec![0.0; n];
                let parse_idx = |x: f64| -> Result<usize, String> {
                    if !x.is_finite() || x < 0.0 || x.fract() != 0.0 || x >= n as f64 {
                        Err("invalid BSM split weights".into())
                    } else {
                        Ok(x as usize)
                    }
                };
                for e in l.events.chunks_exact(4) {
                    let (p, a, b) = (parse_idx(e[0])?, parse_idx(e[1])?, parse_idx(e[2])?);
                    if !e[3].is_finite() || e[3] < 0.0 {
                        return Err("invalid BSM split weights".into());
                    }
                    if e[3] > 0.0 {
                        totals[p] += e[3];
                        splits[p].push(Split {
                            left: a,
                            right: b,
                            prior: e[3],
                            kind: split_type(&input.states[p], &input.states[a], &input.states[b]),
                        });
                    }
                }
                if conditioning != 3 && input.states[0].iter().all(|x| *x < 0) {
                    if !splits[0].is_empty() {
                        return Err("empty-state split is supplied by the complete contract".into());
                    }
                    splits[0].push(Split {
                        left: 0,
                        right: 0,
                        prior: 1.0,
                        kind: 4,
                    });
                    totals[0] = 1.0;
                }
                for (v, total) in splits.iter_mut().zip(totals) {
                    v.iter_mut().for_each(|s| s.prior /= total);
                }
                Ok(Layer {
                    generator,
                    q: l.q.to_vec(),
                    backward: ActionPlan::new(&q),
                    forward: ActionPlan::new(&qt),
                    splits,
                    allowed: l.allowed.to_vec(),
                    dispersal: l.dispersal.to_vec(),
                })
            })
            .collect::<Result<_, String>>()?;
        let mut branches = (0..nn).map(|_| Vec::new()).collect::<Vec<Vec<Segment>>>();
        for spec in &input.segments {
            if spec.node >= nn
                || spec.node == input.root
                || spec.layer >= layers.len()
                || spec.length < 0.0
                || !spec.length.is_finite()
                || !spec.start_age.is_finite()
            {
                return Err("invalid BSM branch segment".into());
            }
            branches[spec.node].push(Segment {
                spec: spec.clone(),
                target: Vec::new(),
                rows: (0..n).map(|_| OnceLock::new()).collect(),
            });
        }

        let mut ready = vec![false; nn];
        let mut down = vec![vec![0.0; n]; nn];
        let mut branch_likes = vec![vec![0.0; n]; nn];
        let propagate = |value: &[f64], branch: &mut Vec<Segment>| -> Result<Vec<f64>, String> {
            let mut result = value.to_vec();
            for segment in branch.iter_mut().rev() {
                let layer = &layers[segment.spec.layer];
                for (x, a) in result.iter_mut().zip(&layer.allowed) {
                    if *a == 0.0 {
                        *x = 0.0;
                    }
                }
                segment.target = result.clone();
                result = layer.backward.apply(segment.spec.length, &result)?;
                for (x, a) in result.iter_mut().zip(&layer.allowed) {
                    *x = if *a == 0.0 { 0.0 } else { x.max(0.0) };
                }
                normalize(&mut result)?;
            }
            Ok(result)
        };
        for i in 0..nt {
            for s in 0..n {
                down[i][s] = match conditioning {
                    1 => {
                        if s == 0 {
                            0.0
                        } else {
                            1.0
                        }
                    }
                    2 => 1.0,
                    _ => input.tips[i + nt * s],
                };
            }
            normalize(&mut down[i])?;
            branch_likes[i] = propagate(&down[i], &mut branches[i])?;
            ready[i] = true;
        }
        for &[node, left, right] in &input.steps {
            if node < nt
                || node >= nn
                || left >= nn
                || right >= nn
                || ready[node]
                || left == right
                || !ready[left]
                || !ready[right]
            {
                return Err("BSM steps must be a rooted binary pruning order".into());
            }
            let layer = &layers[input.layers_at_node[node]];
            for s in 0..n {
                let a = &branch_likes[left];
                let b = &branch_likes[right];
                let v = if input.hooks[node] != 0.0 {
                    a[s] * b[s]
                } else {
                    layer.splits[s]
                        .iter()
                        .map(|split| split.prior * a[split.left] * b[split.right])
                        .sum()
                };
                down[node][s] = v * input.fixes[node + nn * s];
            }
            normalize(&mut down[node])?;
            if node != input.root {
                branch_likes[node] = propagate(&down[node], &mut branches[node])?;
            }
            ready[node] = true;
        }
        if input.steps.last().map(|x| x[0]) != Some(input.root) || !ready.iter().all(|x| *x) {
            return Err("incomplete BSM tree".into());
        }

        if !root_weights.is_empty() {
            for (p, w) in down[input.root].iter_mut().zip(root_weights) {
                *p *= w;
            }
            if down[input.root].iter().sum::<f64>() <= 0.0 {
                return Err("zero root conditional mass".into());
            }
        }
        Ok(Self {
            input,
            layers,
            branches,
            down,
            branch_likes,
            na,
            conditioning,
        })
    }

    fn sample_one(
        &self,
        map: usize,
        seed: u64,
        rejection: bool,
        history: bool,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<usize>), String> {
        let input = &self.input;
        let n = input.states.len();
        let nn = input.ages.len();
        let layers = &self.layers;
        let branches = &self.branches;
        let branch_likes = &self.branch_likes;
        let down = &self.down;
        let na = self.na;
        let mut rng = Xoshiro256PlusPlus::from_key(seed, map as u64, 0, 0);
        let mut state = vec![0; nn];
        let mut clado = Vec::new();
        let mut ana = Vec::new();
        state[input.root] = positive(&down[input.root], &mut rng)?;
        for &[node, left, right] in input.steps.iter().rev() {
            let layer_id = input.layers_at_node[node];
            let layer = &layers[layer_id];
            let parent = state[node];
            let split = if input.hooks[node] != 0.0 {
                Split {
                    left: parent,
                    right: parent,
                    prior: 1.0,
                    kind: 4,
                }
            } else {
                let rows = &layer.splits[parent];
                let mass: Vec<_> = rows
                    .iter()
                    .map(|s| s.prior * branch_likes[left][s.left] * branch_likes[right][s.right])
                    .collect();
                rows[positive(&mass, &mut rng)?]
            };
            let mut source = -1.0;
            let mut destination = -1.0;
            if history && split.kind == 3 {
                if let Some(area) = areas(&input.states[split.left])
                    .chain(areas(&input.states[split.right]))
                    .find(|&a| !input.states[parent].contains(&(a as i32)))
                {
                    destination = area as f64;
                    source =
                        select_source(&input.states[parent], area, &layer.dispersal, na, &mut rng);
                }
            }
            if history && (self.conditioning == 3 || !input.states[parent].iter().all(|x| *x < 0)) {
                clado.extend([
                    map as f64,
                    node as f64,
                    left as f64,
                    right as f64,
                    parent as f64,
                    split.left as f64,
                    split.right as f64,
                    input.ages[node],
                    layer_id as f64,
                    split.kind as f64,
                    source,
                    destination,
                ]);
            }
            for (child, start) in [(left, split.left), (right, split.right)] {
                let mut current = start;
                for (j, segment) in branches[child].iter().enumerate() {
                    let lid = segment.spec.layer;
                    let g = &layers[lid];
                    let t = segment.spec.length;
                    let probs = segment.rows[current]
                        .get_or_init(|| {
                            let mut unit = vec![0.0; n];
                            unit[current] = 1.0;
                            let mut row = g.forward.apply(t, &unit)?;
                            for (x, &p) in row.iter_mut().zip(&segment.target) {
                                *x = x.max(0.0) * p;
                            }
                            Ok(row)
                        })
                        .as_ref()
                        .map_err(Clone::clone)?;
                    let end = positive(probs, &mut rng)?;
                    let mut branch_rng = Xoshiro256PlusPlus::from_key(
                        seed,
                        map as u64,
                        lid as u64 + 1,
                        ((child as u64) << 32) | (j as u64),
                    );
                    if history {
                        for ev in bridge(g, current, end, t, n, rejection, &mut branch_rng)? {
                            let before = &input.states[ev.from];
                            let after = &input.states[ev.to];
                            let gain = areas(after).find(|&a| !before.contains(&(a as i32)));
                            let loss = areas(before).find(|&a| !after.contains(&(a as i32)));
                            let (kind, src, dst) = match (gain, loss) {
                                (Some(a), Some(b)) => (2, b as f64, a as f64),
                                (Some(a), None) => (
                                    0,
                                    select_source(before, a, &g.dispersal, na, &mut rng),
                                    a as f64,
                                ),
                                (None, Some(a)) => (1, a as f64, -1.0),
                                _ => (7, -1.0, -1.0), // change in trait only
                            };
                            ana.extend([
                                map as f64,
                                child as f64,
                                lid as f64,
                                ev.from as f64,
                                ev.to as f64,
                                segment.spec.start_age - ev.time,
                                ev.time,
                                kind as f64,
                                src,
                                dst,
                            ]);
                        }
                    }
                    current = end;
                }
                state[child] = current;
            }
        }
        Ok((clado, ana, state))
    }

    fn validate_range(&self, start: u64, count: usize, threads: usize) -> Result<(), String> {
        if count == 0
            || count > 1_000_000
            || threads == 0
            || threads > 256
            || start == 0
            || start.checked_add(count as u64).is_none()
            || start + count as u64 - 1 > 9_007_199_254_740_991
        {
            return Err("invalid sampler range or threads".into());
        }
        Ok(())
    }

    pub fn draw_bsm(
        &self,
        start: u64,
        count: usize,
        seed: u64,
        threads: usize,
        rejection: bool,
    ) -> Result<Vec<Vec<f64>>, String> {
        self.validate_range(start, count, threads)?;
        if self.conditioning == 1 || self.conditioning == 2 {
            return Err("simulation handle cannot draw observed-data BSM".into());
        }
        let simulate = |map: usize| self.sample_one(map, seed, rejection, true);
        let range = start as usize..start as usize + count;
        let all: Vec<_> = if threads == 1 {
            range.map(simulate).collect()
        } else {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| e.to_string())?;
            pool.install(|| range.into_par_iter().map(simulate).collect())
        };
        let mut out = vec![Vec::new(), Vec::new()];
        for result in all {
            let (c, a, _) = result?;
            out[0].extend(c);
            out[1].extend(a);
        }
        Ok(out)
    }

    pub fn simulate_tips(
        &self,
        start: u64,
        count: usize,
        seed: u64,
        threads: usize,
    ) -> Result<Vec<Vec<f64>>, String> {
        self.validate_range(start, count, threads)?;
        if self.conditioning != 1 && self.conditioning != 2 {
            return Err("observed-data handle cannot simulate datasets".into());
        }
        let seed = seed ^ 0x5349_4D5F_4247_4233;
        let simulate = |map: usize| self.sample_one(map, seed, false, false).map(|x| x.2);
        let range = start as usize..start as usize + count;
        let all: Vec<_> = if threads == 1 {
            range.map(simulate).collect()
        } else {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| e.to_string())?;
            pool.install(|| range.into_par_iter().map(simulate).collect())
        };
        let mut tips = Vec::with_capacity(count * self.input.tip_count);
        let mut nodes = Vec::with_capacity(count * self.input.ages.len());
        for result in all {
            let states = result?;
            tips.extend(states[..self.input.tip_count].iter().map(|x| *x as f64));
            nodes.extend(states.iter().map(|x| *x as f64));
        }
        Ok(vec![tips, nodes])
    }
}
