use crate::Event;
use rayon::prelude::*;
use std::collections::VecDeque;

const NEG_INF: f64 = f64::NEG_INFINITY;
const MAX_STEPS: usize = 100_000;
const MAX_PLAN_CELLS: usize = 4_000_000;
const LOG_REL_TAIL: f64 = -54.0 * std::f64::consts::LN_2;

pub struct Xoshiro256PlusPlus([u64; 4]);

impl Xoshiro256PlusPlus {
    pub fn from_key(seed: u64, map: u64, stratum: u64, step: u64) -> Self {
        // Each lane is a permutation. A nonzero map gives a nonzero RNG state.
        fn mix(mut x: u64) -> u64 {
            x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
            x ^ (x >> 31)
        }
        let mut stream = Self([
            mix(seed ^ 0x9e3779b97f4a7c15),
            mix(map),
            mix(stratum),
            mix(step),
        ]);
        // Mix the four key lanes before drawing: the first raw xoshiro output
        // depends on only two state words.
        for _ in 0..16 {
            stream.uniform();
        }
        stream
    }

    pub fn uniform(&mut self) -> f64 {
        let s = &mut self.0;
        let output = s[0].wrapping_add(s[3]).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        ((output >> 12) as f64 + 0.5) * (1.0 / (1_u64 << 52) as f64)
    }
}

fn log_add(a: f64, b: f64) -> f64 {
    if a == NEG_INF {
        return b;
    }
    if b == NEG_INF {
        return a;
    }
    a.max(b) + (-(a - b).abs()).exp().ln_1p()
}

fn pick(weights: &[f64], rng: &mut Xoshiro256PlusPlus) -> usize {
    let max = weights.iter().copied().fold(NEG_INF, f64::max);
    let mass = |value: f64| (value - max).exp();
    let total: f64 = weights.iter().copied().map(mass).sum();
    let mut threshold = total * rng.uniform();
    let mut last = 0;
    for (index, &weight) in weights.iter().enumerate() {
        if weight.is_finite() {
            last = index;
            threshold -= mass(weight);
            if threshold <= 0.0 {
                break;
            }
        }
    }
    last
}

/// Sparse jump-chain rows. Preserve the diagonal waiting rate and normalize
/// the off-diagonal jump probabilities, as the upstream rejection sampler
/// does. The tolerance admits generator rows rounded in LegacyF32.
pub struct Generator {
    rows: Vec<Vec<(usize, f64)>>,
    rate: f64,
}

impl Generator {
    pub fn new(q: &[f64], n: usize) -> Result<Self, String> {
        if n == 0 || Some(q.len()) != n.checked_mul(n) {
            return Err("invalid generator dimensions".into());
        }
        let mut exits = Vec::with_capacity(n);
        let mut rates = Vec::with_capacity(n);
        let mut rows = Vec::with_capacity(n);
        for i in 0..n {
            let mut edges = Vec::new();
            let mut exit = 0.0_f64;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let value = q[i + j * n];
                if !value.is_finite() || value < 0.0 {
                    return Err("generator has a negative or non-finite off-diagonal rate".into());
                }
                if value > 0.0 {
                    edges.push((j, value));
                    exit += value;
                }
            }
            let diag = -q[i * (n + 1)];
            if !exit.is_finite()
                || !diag.is_finite()
                || diag < 0.0
                || (diag - exit).abs() > diag.max(exit) * (8.0 * f32::EPSILON as f64)
            {
                return Err(format!(
                    "generator row {} has waiting rate {diag} and outgoing total {exit}",
                    i + 1
                ));
            }
            exits.push(exit);
            rates.push(diag);
            rows.push(edges);
        }
        let rate = rates.iter().copied().fold(0.0_f64, f64::max);
        if rate > 0.0 {
            for (i, edges) in rows.iter_mut().enumerate() {
                for (_, value) in edges.iter_mut() {
                    *value = (value.ln() - exits[i].ln()) + (rates[i].ln() - rate.ln());
                }
                let stay = rate - rates[i];
                if stay > 0.0 {
                    edges.push((i, stay.ln() - rate.ln()));
                }
            }
        }
        Ok(Self { rows, rate })
    }

    fn reachable(&self, from: usize, to: usize) -> bool {
        let mut seen = vec![false; self.rows.len()];
        let mut todo = VecDeque::from([from]);
        seen[from] = true;
        while let Some(i) = todo.pop_front() {
            if i == to {
                return true;
            }
            for &(j, _) in &self.rows[i] {
                if !seen[j] {
                    seen[j] = true;
                    todo.push_back(j);
                }
            }
        }
        false
    }

    pub fn sample(
        &self,
        from: usize,
        to: usize,
        duration: f64,
        rng: &mut Xoshiro256PlusPlus,
    ) -> Result<Vec<Event>, String> {
        let n = self.rows.len();
        if from >= n || to >= n || !duration.is_finite() || duration < 0.0 {
            return Err("invalid bridge endpoint or branch length".into());
        }
        if duration == 0.0 || self.rate == 0.0 {
            return if from == to {
                Ok(Vec::new())
            } else {
                Err("unreachable endpoint".into())
            };
        }
        if !self.reachable(from, to) {
            return Err("unreachable endpoint".into());
        }
        let log_lambda = self.rate.ln() + duration.ln();
        let lambda = self.rate * duration;
        let limit = MAX_STEPS.min(MAX_PLAN_CELLS / n);
        if !lambda.is_finite() || lambda >= limit as f64 {
            return Err("uniformization work limit exceeded".into());
        }

        // Backward messages h[k][i] = log((R^k)[i,to]).
        // Weight for N=k is the Poisson term times exp(h[k][from]).
        // All terms omit the same exp(-lambda) factor. This cancels in
        // sampling and in the relative endpoint-tail bound.
        let mut current = vec![NEG_INF; n];
        current[to] = 0.0;
        let mut messages = Vec::new();
        let mut weights = Vec::new();
        let mut log_total = NEG_INF;
        let mut log_poisson = 0.0;
        for k in 0..=limit {
            let weight = log_poisson + current[from];
            weights.push(weight);
            log_total = log_add(log_total, weight);
            messages.push(current);
            let log_next = log_poisson + log_lambda - ((k + 1) as f64).ln();
            let ratio = lambda / (k + 2) as f64;
            if ratio < 1.0
                && log_total.is_finite()
                && log_next - (-ratio).ln_1p() <= log_total + LOG_REL_TAIL
            {
                let count = pick(&weights, rng);
                let mut times: Vec<_> = (0..count).map(|_| rng.uniform() * duration).collect();
                times.sort_unstable_by(f64::total_cmp);
                let mut state = from;
                let mut events = Vec::new();
                let mut next_weights = Vec::new();
                for (i, time) in times.into_iter().enumerate() {
                    let row = &self.rows[state];
                    let suffix = &messages[count - i - 1];
                    next_weights.clear();
                    next_weights.extend(row.iter().map(|&(next, log_p)| log_p + suffix[next]));
                    let next = row[pick(&next_weights, rng)].0;
                    if state != next {
                        events.push(Event {
                            from: state,
                            to: next,
                            time,
                        });
                    }
                    state = next;
                }
                debug_assert_eq!(state, to);
                return Ok(events);
            }
            log_poisson = log_next;
            let previous = &messages[k];
            current = self
                .rows
                .iter()
                .map(|row| {
                    row.iter().fold(NEG_INF, |acc, &(j, log_p)| {
                        log_add(acc, log_p + previous[j])
                    })
                })
                .collect();
        }
        Err("uniformization tail bound did not converge within the work limit".into())
    }
}

pub struct BridgeJob {
    pub generator: usize,
    pub from: usize,
    pub to: usize,
    pub duration: f64,
    pub map: u64,
    pub stratum: u64,
    pub step: u64,
    pub node: u64,
}

pub fn sample_batch(
    generators: &[Generator],
    jobs: &[BridgeJob],
    seed: u64,
    threads: usize,
) -> Result<Vec<Vec<Event>>, String> {
    let sample = |job: &BridgeJob| {
        if job.map == 0 {
            return Err("map index must be positive".into());
        }
        let generator = generators
            .get(job.generator)
            .ok_or("invalid generator index")?;
        let mut rng = Xoshiro256PlusPlus::from_key(seed, job.map, job.stratum, job.step);
        generator
            .sample(job.from, job.to, job.duration, &mut rng)
            .map_err(|error| {
                format!(
                    "{error} (map {}, stratum {}, node {}, start {}, end {}, length {})",
                    job.map,
                    job.stratum,
                    job.node,
                    job.from + 1,
                    job.to + 1,
                    job.duration
                )
            })
    };
    if threads == 1 || jobs.len() <= 1 {
        return jobs.iter().map(sample).collect();
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| e.to_string())?;
    // Collect in submission order so even diagnostics are scheduling invariant.
    let results: Vec<_> = pool.install(|| jobs.par_iter().map(sample).collect());
    results.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn very_rare_bridge_uses_reachability_without_probability_underflow() {
        // Total probability of this two-jump path is below binary64 range.
        let q = [-1e-200, 0.0, 0.0, 1e-200, -1e-200, 0.0, 0.0, 1e-200, 0.0];
        let g = Generator::new(&q, 3).unwrap();
        let mut rng = Xoshiro256PlusPlus::from_key(9, 1, 0, 1);
        let history = g.sample(0, 2, 1.0, &mut rng).unwrap();
        assert_eq!(history.len(), 2);
        assert_eq!((history[0].from, history[0].to, history[1].to), (0, 1, 2));
        assert!(history[0].time <= history[1].time);
        assert_eq!(
            g.sample(2, 0, 1.0, &mut rng).unwrap_err(),
            "unreachable endpoint"
        );
    }

    #[test]
    fn large_poisson_mean_has_a_closed_history() {
        let g = Generator::new(&[-800.0, 800.0, 800.0, -800.0], 2).unwrap();
        let mut rng = Xoshiro256PlusPlus::from_key(9, 1, 0, 1);
        let history = g.sample(0, 1, 1.0, &mut rng).unwrap();
        assert_eq!(history.len() % 2, 1);
        assert_eq!(history.last().unwrap().to, 1);
    }
}
