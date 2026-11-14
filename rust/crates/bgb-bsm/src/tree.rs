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
struct Layer<'a> {
    generator: Generator,
    q: &'a [f64],
    backward: ActionPlan,
    forward: ActionPlan,
    splits: Vec<Vec<Split>>,
    allowed: &'a [f64],
    dispersal: &'a [f64],
}
struct Segment<'a> {
    spec: &'a SegmentInput,
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
    layer: &Layer<'_>,
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
        layer.q,
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

/// Input steps are in pruningwise order. Branch segments are ordered from
/// older to younger; a shared transition cache is populated per starting state.
pub fn sample_tree(
    input: TreeInput<'_>,
    maps: usize,
    seed: u64,
    threads: usize,
    rejection: bool,
) -> Result<Vec<Vec<f64>>, String> {
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
        || maps == 0
        || threads == 0
        || input
            .layers_at_node
            .iter()
            .any(|&l| l >= input.layers.len())
    {
        return Err("invalid BSM tree dimensions".into());
    }
    let layers: Vec<_> = input
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
            for (v, total) in splits.iter_mut().zip(totals) {
                v.iter_mut().for_each(|s| s.prior /= total);
            }
            Ok(Layer {
                generator,
                q: l.q,
                backward: ActionPlan::new(&q),
                forward: ActionPlan::new(&qt),
                splits,
                allowed: l.allowed,
                dispersal: l.dispersal,
            })
        })
        .collect::<Result<_, String>>()?;
    let mut branches = (0..nn).map(|_| Vec::new()).collect::<Vec<Vec<Segment>>>();
    for spec in input.segments {
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
            spec,
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
            for (x, a) in result.iter_mut().zip(layer.allowed) {
                if *a == 0.0 {
                    *x = 0.0;
                }
            }
            segment.target = result.clone();
            result = layer.backward.apply(segment.spec.length, &result)?;
            for (x, a) in result.iter_mut().zip(layer.allowed) {
                *x = if *a == 0.0 { 0.0 } else { x.max(0.0) };
            }
            normalize(&mut result)?;
        }
        Ok(result)
    };
    for i in 0..nt {
        for s in 0..n {
            down[i][s] = input.tips[i + nt * s];
        }
        normalize(&mut down[i])?;
        branch_likes[i] = propagate(&down[i], &mut branches[i])?;
        ready[i] = true;
    }
    for &[node, left, right] in input.steps {
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

    let simulate = |map: usize| -> Result<(Vec<f64>, Vec<f64>), String> {
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
            if split.kind == 3 {
                if let Some(area) = areas(&input.states[split.left])
                    .chain(areas(&input.states[split.right]))
                    .find(|&a| !input.states[parent].contains(&(a as i32)))
                {
                    destination = area as f64;
                    source =
                        select_source(&input.states[parent], area, layer.dispersal, na, &mut rng);
                }
            }
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
                    for ev in bridge(g, current, end, t, n, rejection, &mut branch_rng)? {
                        let before = &input.states[ev.from];
                        let after = &input.states[ev.to];
                        let gain = areas(after).find(|&a| !before.contains(&(a as i32)));
                        let loss = areas(before).find(|&a| !after.contains(&(a as i32)));
                        let (kind, src, dst) = match (gain, loss) {
                            (Some(a), Some(b)) => (2, b as f64, a as f64),
                            (Some(a), None) => (
                                0,
                                select_source(before, a, g.dispersal, na, &mut rng),
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
                    current = end;
                }
                state[child] = current;
            }
        }
        Ok((clado, ana))
    };
    let all: Vec<_> = if threads == 1 {
        (1..=maps).map(simulate).collect()
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| e.to_string())?;
        pool.install(|| (1..=maps).into_par_iter().map(simulate).collect())
    };
    let mut out = vec![Vec::new(), Vec::new()];
    for result in all {
        let (c, a) = result?;
        out[0].extend(c);
        out[1].extend(a);
    }
    Ok(out)
}
