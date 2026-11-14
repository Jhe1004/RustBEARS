//! Nonnegative uniformization with reachability and a componentwise tail gate.
//! No test against EPS times the largest/global accumulated component is used.
#[derive(Clone)]
pub(super) struct PositivePlan {
    pub rate: f64,
    n: usize,
    entries: Vec<(usize, usize, f64)>,
    col_ptr: Vec<usize>,
    beta: f64,
    one_norm: bool,
}
impl PositivePlan {
    pub fn new(n: usize, entries: Vec<(usize, usize, f64)>, rate: f64) -> Option<Self> {
        let mut rows = vec![0.0; n];
        let mut cols = vec![0.0; n];
        let mut col_ptr = vec![0; n + 1];
        for &(i, j, x) in &entries {
            if !x.is_finite() || x < 0.0 {
                return None;
            }
            rows[i] += x;
            cols[j] += x;
            col_ptr[j + 1] += 1;
        }
        for j in 0..n {
            col_ptr[j + 1] += col_ptr[j];
        }
        let row_norm = rows.into_iter().fold(0.0, f64::max);
        let col_norm = cols.into_iter().fold(0.0, f64::max);
        let one_norm = col_norm < row_norm;
        let beta = row_norm.min(col_norm) * (1.0 + 64.0 * (n as f64 + 1.0) * f64::EPSILON);
        if !beta.is_finite() {
            return None;
        }
        Some(Self {
            rate,
            n,
            entries,
            col_ptr,
            beta,
            one_norm,
        })
    }
    fn reach(&self, v: &[f64]) -> Vec<bool> {
        let mut seen: Vec<_> = v.iter().map(|&x| x > 0.0).collect();
        let mut queue: Vec<_> = (0..self.n).filter(|&i| seen[i]).collect();
        let mut at = 0;
        while at < queue.len() {
            let j = queue[at];
            at += 1;
            for k in self.col_ptr[j]..self.col_ptr[j + 1] {
                let i = self.entries[k].0;
                if !seen[i] {
                    seen[i] = true;
                    queue.push(i);
                }
            }
        }
        seen
    }
    fn nonnegative(&self, t: f64, v: &[f64]) -> Result<(Vec<f64>, usize), String> {
        if self.rate == 0.0 || t == 0.0 {
            return Ok((v.to_vec(), 0));
        }
        let mean = self.rate * t;
        if mean == 0.0 {
            return Err("uniformization time product underflow".into());
        }
        let steps = (mean / 0.5).ceil().max(1.0) as usize;
        if steps > 1_000_000 {
            return Err("uniformization action scaling limit exceeded".into());
        }
        let z = self.rate * (t / steps as f64);
        let weight = (-z).exp();
        let effective = z * self.beta;
        let log_effective = effective.ln();
        let mut out = v.to_vec();
        let mut term = vec![0.0; self.n];
        let mut next = vec![0.0; self.n];
        let mut total = vec![0.0; self.n];
        let mut multiplies = 0;
        for _ in 0..steps {
            let max = out.iter().copied().fold(0.0, f64::max);
            if max == 0.0 {
                return Ok((out, multiplies));
            }
            let log_norm = max.ln()
                + if self.one_norm {
                    out.iter().map(|x| x / max).sum::<f64>().ln()
                } else {
                    0.0
                };
            let required = self.reach(&out);
            for i in 0..self.n {
                term[i] = out[i] * weight;
                total[i] = term[i];
            }
            if self.beta == 0.0 {
                out.copy_from_slice(&total);
                continue;
            }
            let mut log_coefficient = 0.0;
            let mut converged = false;
            for k in 1..=4096 {
                next.fill(0.0);
                for &(i, j, a) in &self.entries {
                    next[i] += a * term[j];
                }
                for i in 0..self.n {
                    next[i] *= z / k as f64;
                    total[i] += next[i];
                }
                multiplies += 1;
                std::mem::swap(&mut next, &mut term);
                if total.iter().any(|x| !x.is_finite()) {
                    return Err("exponential action overflow".into());
                }
                log_coefficient += log_effective - (k as f64).ln();
                let ratio = effective / (k + 2) as f64;
                if ratio >= 1.0 {
                    continue;
                }
                // e^-z ||v|| (z*beta)^(k+1)/(k+1)! / (1-ratio).
                let log_tail = log_norm - z + log_coefficient + log_effective
                    - ((k + 1) as f64).ln()
                    - (-ratio).ln_1p();
                let complete = required
                    .iter()
                    .zip(&total)
                    .all(|(&needed, &x)| !needed || x > 0.0);
                let smallest = total
                    .iter()
                    .copied()
                    .filter(|&x| x > 0.0)
                    .fold(f64::INFINITY, f64::min);
                let relative_gate =
                    complete && log_tail <= smallest.ln() + (f64::EPSILON / 64.0).ln();
                // An absolute underflow floor also permits truly unrepresentable
                // probabilities to remain zero; this is not a global EPS gate.
                let underflow_gate =
                    log_tail <= f64::from_bits(1).ln() - std::f64::consts::LN_2 - 8.0;
                if relative_gate || underflow_gate {
                    converged = true;
                    break;
                }
            }
            if !converged {
                return Err("positive series tail bound not attained".into());
            }
            out.copy_from_slice(&total);
        }
        Ok((out, multiplies))
    }
    pub fn apply(&self, t: f64, v: &[f64]) -> Result<(Vec<f64>, usize), String> {
        if v.iter().all(|&x| x >= 0.0) {
            return self.nonnegative(t, v);
        }
        let plus: Vec<_> = v.iter().map(|x| x.max(0.0)).collect();
        let minus: Vec<_> = v.iter().map(|x| (-x).max(0.0)).collect();
        let (a, ka) = self.nonnegative(t, &plus)?;
        let (b, kb) = self.nonnegative(t, &minus)?;
        Ok((a.iter().zip(b).map(|(a, b)| a - b).collect(), ka + kb))
    }
}
