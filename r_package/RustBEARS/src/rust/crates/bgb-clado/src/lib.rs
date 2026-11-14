//! Port of cladoRcpp's GPL algorithms; original author Nicholas J. Matzke.
//! Preserve event order and mixed-precision conversion sites, not only formulas.
// Indexed loops intentionally preserve the reference reduction order.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]
use bgb_core::PrecisionMode;
mod template;
pub use template::CompiledCladoTemplate;
pub type Events = Vec<(usize, usize, usize, f64)>;
pub struct Model<'a> {
    pub states: &'a [Vec<i32>],
    pub parameters: &'a [f64],
    pub d: &'a [f64],
    pub priors: [&'a [f64]; 4],
    pub maxmin: &'a [f64],
    pub na: usize,
    pub mode: PrecisionMode,
}
impl Model<'_> {
    pub fn validate(&self) -> Result<(), String> {
        let n = self.states.len();
        if n == 0
            || self.parameters.len() != 4
            || self.d.len() != self.na * self.na
            || self.priors.iter().any(|a| a.len() != self.na * self.na)
            || self.maxmin.len() < self.na
            || self
                .states
                .iter()
                .any(|s| s.is_empty() || s.iter().any(|&x| x < 0 || x >= self.na as i32))
        {
            return Err("invalid cladogenesis dimensions or state indices".into());
        }
        Ok(())
    }
    fn rate(&self, k: usize) -> f64 {
        self.mode.round(self.parameters[k])
    }
    fn prior(&self, k: usize, a: usize, b: usize) -> f64 {
        self.priors[k][a - 1 + self.na * (b - 1)]
    }
    fn distance(&self, a: &[i32], b: &[i32]) -> f64 {
        let mut sum = 0.0;
        for &x in a {
            for &y in b {
                sum = self
                    .mode
                    .round(sum + self.d[x as usize + self.na * y as usize]);
            }
        }
        let den = self.mode.round(a.len() as f64) * self.mode.round(b.len() as f64);
        self.mode.round(sum / den)
    }
    /// Tuple keeps the jump-distance factor separate to preserve multiplication order.
    fn dense_factor(&self, anc: usize, left: usize, right: usize) -> (f64, f64) {
        let a = &self.states[anc];
        let l = &self.states[left];
        let r = &self.states[right];
        let asz = a.len();
        let ls = l.len();
        let rs = r.len();
        if ls.min(rs) > self.maxmin[asz - 1] as usize {
            return (0.0, 1.0);
        }
        if ls == asz && rs == asz && l == r && l == a {
            return (self.rate(3) * self.prior(3, asz, asz), 1.0);
        }
        if ls == asz && l == a && (rs < asz || asz == 1) {
            if r.iter().all(|x| !a.contains(x)) {
                return (self.rate(2) * self.prior(2, asz, rs), self.distance(a, r));
            } else {
                return (self.rate(0) * self.prior(0, asz, rs), 1.0);
            }
        }
        if rs == asz && r == a && (ls < asz || asz == 1) {
            if l.iter().all(|x| !a.contains(x)) {
                return (self.rate(2) * self.prior(2, asz, ls), self.distance(a, l));
            } else {
                return (self.rate(0) * self.prior(0, asz, ls), 1.0);
            }
        }
        let mut merged = l.clone();
        merged.extend_from_slice(r);
        merged.sort_unstable(); // merge_int_vectors preserves duplicates.
        if merged == *a {
            return (self.rate(1) * self.prior(1, asz, ls.min(rs)), 1.0);
        }
        (0.0, 1.0)
    }
    pub fn dense(
        &self,
        left: &[f64],
        right: &[f64],
        denom: &[f64],
        kind: u32,
    ) -> Result<Vec<Vec<f64>>, String> {
        self.validate()?;
        let n = self.states.len();
        if left.len() != n || right.len() != n || (kind == 0 && denom.len() != n) {
            return Err("cladogenesis vector dimensions".into());
        }
        let mut rows = vec![0.0; n];
        let mut ilists = vec![Vec::new(); n];
        let mut jlists = vec![Vec::new(); n];
        let mut wlists = vec![Vec::new(); n];
        let cutoff = self.mode.round(1e-30);
        for a in 0..n {
            let mut sum = 0.0;
            for (j, &r) in right.iter().enumerate() {
                for (i, &l) in left.iter().enumerate() {
                    let (b, dist) = self.dense_factor(a, i, j);
                    let w = if kind == 0 {
                        self.mode.round(b * l * r * dist)
                    } else {
                        self.mode.round(b * dist)
                    };
                    sum = self.mode.round(sum + w);
                    if kind == 2 && w.abs() > cutoff {
                        ilists[a].push(i as f64);
                        jlists[a].push(j as f64);
                        wlists[a].push(w);
                    }
                }
            }
            rows[a] = if kind == 0 { sum / denom[a] } else { sum };
        }
        if kind != 2 {
            return Ok(vec![rows]);
        }
        let mut out = ilists;
        out.extend(jlists);
        out.extend(wlists);
        Ok(out)
    }
    pub fn faster(&self) -> Result<Events, String> {
        self.validate()?;
        let n = self.states.len();
        let mut out = Vec::new();
        let cut = self.mode.round(1e-30);
        let m = self.mode;
        let group_count = self.states.iter().map(|s| s.len()).max().unwrap_or(0);
        // Full sympatry, followed by subset events and their exact mirror ordering.
        if self.rate(3).abs() > cut {
            for (a, s) in self.states.iter().enumerate() {
                let k = s.len();
                if k <= self.maxmin[k - 1] as usize {
                    let w = m.round(self.rate(3) * self.prior(3, k, k));
                    if w.abs() > cut {
                        out.push((a, a, a, w));
                    }
                }
            }
        }
        let mirror = |out: &mut Events, events: &Events| {
            out.extend(events.iter().copied());
            out.extend(events.iter().map(|&(a, i, j, w)| (a, j, i, w)));
        };
        if self.rate(0).abs() > cut {
            let mut events = Vec::new();
            for (a, s) in self.states.iter().enumerate() {
                let k = s.len();
                if k < 2 {
                    continue;
                }
                for size in 1..=self.maxmin[k - 1] as usize {
                    if size >= group_count {
                        continue;
                    }
                    for (j, r) in self
                        .states
                        .iter()
                        .enumerate()
                        .filter(|(_, r)| r.len() == size)
                    {
                        if r != s && r.iter().all(|x| s.contains(x)) {
                            let w = m.round(self.rate(0) * self.prior(0, k, size));
                            if w.abs() > cut {
                                events.push((a, a, j, w));
                            }
                        }
                    }
                }
            }
            mirror(&mut out, &events);
        }
        if self.rate(1).abs() > cut {
            let mut equal = Vec::new();
            let mut unequal = Vec::new();
            for (a, s) in self.states.iter().enumerate() {
                let k = s.len();
                if k < 2 {
                    continue;
                }
                for size in 1..=self.maxmin[k - 1] as usize {
                    if size >= group_count {
                        continue;
                    }
                    for (j, r) in self
                        .states
                        .iter()
                        .enumerate()
                        .filter(|(_, r)| r.len() == size)
                    {
                        if r.iter().all(|x| s.contains(x)) {
                            let l: Vec<_> = s.iter().filter(|x| !r.contains(x)).copied().collect();
                            if l.is_empty() {
                                continue;
                            }
                            if let Some(i) = self.states.iter().position(|v| *v == l) {
                                let w = m.round(self.rate(1) * self.prior(1, k, size));
                                if w.abs() > cut {
                                    if l.len() == r.len() {
                                        equal.push((a, i, j, w));
                                    } else {
                                        unequal.push((a, i, j, w));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            out.extend(equal);
            mirror(&mut out, &unequal);
        }
        if self.rate(2).abs() > cut {
            let mut events = Vec::new();
            for (a, s) in self.states.iter().enumerate() {
                let k = s.len();
                if k > 2000 {
                    continue;
                }
                for size in 1..=self.maxmin[k - 1] as usize {
                    if size >= group_count {
                        continue;
                    }
                    for (j, r) in self
                        .states
                        .iter()
                        .enumerate()
                        .filter(|(_, r)| r.len() == size)
                    {
                        if r.iter().all(|x| !s.contains(x)) {
                            let dist = self.distance(s, r);
                            let w = m.round(self.rate(2) * self.prior(2, k, size) * dist);
                            if w.abs() > cut {
                                events.push((a, a, j, w));
                            }
                        }
                    }
                }
            }
            mirror(&mut out, &events);
        }
        let _ = n;
        Ok(out)
    }
}
pub fn checked_coo_index(value: f64, limit: usize, label: &str) -> Result<usize, String> {
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 || value >= limit as f64 {
        return Err(format!("COO {label} index out of range"));
    }
    Ok(value as usize)
}
pub fn rowsums(indices: &[f64], weights: &[f64], n: usize) -> Result<Vec<f64>, String> {
    if indices.len() != weights.len() {
        return Err("COO length mismatch".into());
    }
    let mut out = vec![0.0; n];
    for (&i, &w) in indices.iter().zip(weights) {
        let idx = checked_coo_index(i, n, "ancestor")?;
        out[idx] += w;
    }
    Ok(out)
}
pub fn split(
    left: &[f64],
    right: &[f64],
    anc: &[f64],
    li: &[f64],
    ri: &[f64],
    weights: &[f64],
    denom: &[f64],
    mode: PrecisionMode,
) -> Result<Vec<f64>, String> {
    let n = anc.len();
    if li.len() != n || ri.len() != n || weights.len() != n {
        return Err("COO buffer mismatch".into());
    }
    let mut out = vec![0.0; denom.len()];
    for k in 0..n {
        let a = checked_coo_index(anc[k], out.len(), "ancestor")?;
        let i = checked_coo_index(li[k], left.len(), "left")?;
        let j = checked_coo_index(ri[k], right.len(), "right")?;
        out[a] += mode.round(left[i] * right[j] * weights[k]);
    }
    for (i, x) in out.iter_mut().enumerate() {
        *x = mode.round(*x / denom[i]);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coo_helpers_reject_negative_nan_fractional_and_oob_indices() {
        for bad in [-1.0, -0.1, 0.5, f64::NAN, f64::INFINITY, 3.0] {
            assert!(rowsums(&[bad], &[1.0], 3).is_err());
            assert!(
                split(
                    &[1.0, 1.0],
                    &[1.0, 1.0],
                    &[bad],
                    &[0.0],
                    &[0.0],
                    &[1.0],
                    &[1.0, 1.0, 1.0],
                    PrecisionMode::StrictF64,
                )
                .is_err()
            );
            assert!(
                split(
                    &[1.0, 1.0],
                    &[1.0, 1.0],
                    &[0.0],
                    &[bad],
                    &[0.0],
                    &[1.0],
                    &[1.0, 1.0, 1.0],
                    PrecisionMode::StrictF64,
                )
                .is_err()
            );
        }
    }
}
