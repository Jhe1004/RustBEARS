//! State ordering and legacy arithmetic. GPL-2.0-or-later.
pub const M0_ABI_VERSION: u32 = 1;
/// Bit operations are an internal fast path; arbitrary public Vec states survive.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum StateBitmask {
    Small(u64),
    Medium(u128),
    General(Vec<i32>),
}
impl StateBitmask {
    pub fn new(state: &[i32]) -> Self {
        let mut bits = 0u128;
        for &x in state {
            if !(0..128).contains(&x) || bits & (1u128 << x) != 0 {
                return Self::General(state.to_vec());
            }
            bits |= 1u128 << x;
        }
        if bits <= u64::MAX as u128 {
            Self::Small(bits as u64)
        } else {
            Self::Medium(bits)
        }
    }
    fn bits(&self) -> Option<u128> {
        match self {
            Self::Small(x) => Some(*x as u128),
            Self::Medium(x) => Some(*x),
            _ => None,
        }
    }
    pub fn pop_count(&self) -> usize {
        match self {
            Self::General(x) => x.len(),
            _ => self.bits().unwrap().count_ones() as usize,
        }
    }
    pub fn subset_of(&self, other: &Self) -> bool {
        match (self.bits(), other.bits()) {
            (Some(a), Some(b)) => a & b == a,
            _ => self.values().iter().all(|x| other.values().contains(x)),
        }
    }
    pub fn disjoint(&self, other: &Self) -> bool {
        match (self.bits(), other.bits()) {
            (Some(a), Some(b)) => a & b == 0,
            _ => self.values().iter().all(|x| !other.values().contains(x)),
        }
    }
    fn values(&self) -> Vec<i32> {
        match self {
            Self::General(x) => x.clone(),
            _ => {
                let mut b = self.bits().unwrap();
                let mut out = Vec::new();
                while b != 0 {
                    let k = b.trailing_zeros();
                    out.push(k as i32);
                    b &= b - 1;
                }
                out
            }
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PrecisionMode {
    LegacyF32,
    StrictF64,
}
impl PrecisionMode {
    #[inline]
    pub fn round(self, x: f64) -> f64 {
        match self {
            Self::LegacyF32 => (x as f32) as f64,
            Self::StrictF64 => x,
        }
    }
}
pub fn choose(n: usize, k: usize) -> Result<usize, String> {
    if k > n {
        return Ok(0);
    }
    let k = k.min(n - k);
    let mut v = 1u128;
    for i in 0..k {
        v = v
            .checked_mul((n - i) as u128)
            .ok_or("combination overflow")?
            / (i + 1) as u128;
    }
    usize::try_from(v).map_err(|_| "combination overflow".into())
}
pub fn combinations(n: usize, k: usize, limit: usize) -> Result<Vec<Vec<i32>>, String> {
    if k == 0 || k > n {
        return Err("invalid combination size".into());
    }
    let count = choose(n, k)?;
    if count > limit {
        return Err("combination limit exceeded".into());
    }
    let mut out = Vec::with_capacity(count);
    let mut a: Vec<_> = (0..k).collect();
    loop {
        out.push(a.iter().map(|&x| x as i32).collect());
        let mut pos = k;
        while pos > 0 && a[pos - 1] == n - k + pos - 1 {
            pos -= 1
        }
        if pos == 0 {
            break;
        }
        a[pos - 1] += 1;
        for j in pos..k {
            a[j] = a[j - 1] + 1;
        }
    }
    Ok(out)
}
pub fn states(n: usize, max: usize, null: bool) -> Result<Vec<Vec<i32>>, String> {
    if n > 1024 {
        return Err("too many areas".into());
    }
    let mut out = Vec::new();
    if null {
        out.push(vec![-1])
    }
    for k in 1..=max.min(n) {
        out.append(&mut combinations(n, k, 10_000_000)?);
        if out.len() > 10_000_000 {
            return Err("state limit exceeded".into());
        }
    }
    Ok(out)
}
pub fn unpack_states(flat: &[f64], lengths: &[f64]) -> Result<Vec<Vec<i32>>, String> {
    let mut out = Vec::with_capacity(lengths.len());
    let mut off = 0usize;
    for &len_val in lengths {
        if !len_val.is_finite()
            || len_val < 1.0
            || len_val > (flat.len() as f64)
            || len_val.fract() != 0.0
        {
            return Err("empty or invalid state".into());
        }
        let n = len_val as usize;
        if off.checked_add(n).is_none_or(|end| end > flat.len()) {
            return Err("state buffer mismatch".into());
        }
        let slice = &flat[off..off + n];
        if n == 1 && (slice[0].is_nan() || slice[0] == -1.0) {
            out.push(vec![-1]);
        } else {
            let mut state = Vec::with_capacity(n);
            for &x in slice {
                if !x.is_finite() || x < 0.0 || x > (i32::MAX as f64) || x.fract() != 0.0 {
                    return Err("invalid state area index".into());
                }
                state.push(x as i32);
            }
            out.push(state);
        }
        off += n;
    }
    if off != flat.len() {
        return Err("state buffer mismatch".into());
    }
    Ok(out)
}
/// Clean-room Rust solver for discrete maximum-entropy distributions under linear
/// moment constraints via Improved Iterative Scaling (Della Pietra, Della Pietra
/// & Lafferty 1997, *IEEE TPAMI* 19:380-393). In `LegacyF32` compatibility mode,
/// per-state step differences are rounded to `f32` before testing `tol`.
pub fn maxent(
    states: &[f64],
    n: usize,
    k: usize,
    target: &[f64],
    prior: &[f64],
    tol: f64,
    mode: PrecisionMode,
) -> Result<Vec<Vec<f64>>, String> {
    if n == 0
        || k == 0
        || states.len() != n * k
        || target.len() != k
        || prior.len() != n
        || !tol.is_finite()
        || tol <= 0.0
    {
        return Err("invalid maxent dimensions or tolerance".into());
    }
    let mut weights = prior.to_vec();
    let mut feature_totals = vec![0.0; k];
    for feat in 0..k {
        let col = &states[feat * n..(feat + 1) * n];
        feature_totals[feat] = col.iter().copied().sum();
    }
    let mut expected_moments = vec![0.0; k];
    let mut step_scales = vec![0.0; k];
    let mut next_weights = vec![0.0; n];
    let mut iterations = 0usize;
    let mut max_step_delta = f64::INFINITY;
    while max_step_delta > tol {
        iterations += 1;
        if iterations > 1_000_000 {
            return Err("maxent iterative scaling did not converge within 1000000 iterations".into());
        }
        for feat in 0..k {
            let col = &states[feat * n..(feat + 1) * n];
            let moment: f64 = weights.iter().zip(col.iter()).map(|(&w, &s)| w * s).sum();
            expected_moments[feat] = moment;
            let rhs = target[feat];
            let fsum = feature_totals[feat];
            if !moment.is_finite()
                || moment <= 0.0
                || !rhs.is_finite()
                || rhs <= 0.0
                || !fsum.is_finite()
                || fsum <= 0.0
            {
                return Err("maxent infeasible non-positive moment or feature total".into());
            }
            step_scales[feat] = (rhs / moment).ln() / fsum;
        }
        let mut partition_sum = 0.0;
        for state_idx in 0..n {
            let mut linear_predictor = 0.0;
            for feat in 0..k {
                linear_predictor += step_scales[feat] * states[state_idx + n * feat];
            }
            let unnormalized = linear_predictor.exp() * weights[state_idx];
            next_weights[state_idx] = unnormalized;
            partition_sum += unnormalized;
        }
        if !partition_sum.is_finite() || partition_sum <= 0.0 {
            return Err("maxent partition function underflowed or overflowed".into());
        }
        max_step_delta = 0.0;
        for state_idx in 0..n {
            let normalized = next_weights[state_idx] / partition_sum;
            let delta = mode.round((normalized - weights[state_idx]).abs());
            if delta > max_step_delta {
                max_step_delta = delta;
            }
            weights[state_idx] = normalized;
        }
    }
    let neg_entropy: f64 = weights
        .iter()
        .copied()
        .filter(|&w| w > 0.0)
        .map(|w| w * w.ln())
        .sum();
    Ok(vec![
        weights,
        vec![-neg_entropy, iterations as f64],
        expected_moments,
    ])
}
#[cfg(test)]
mod tests {
    // Indexed loops intentionally preserve the reference reduction order.
    #![allow(clippy::needless_range_loop, clippy::too_many_arguments)]
    use super::*;
    #[test]
    fn ordering() {
        assert_eq!(
            combinations(4, 2, 10).unwrap(),
            vec![
                vec![0, 1],
                vec![0, 2],
                vec![0, 3],
                vec![1, 2],
                vec![1, 3],
                vec![2, 3]
            ]
        );
        assert_eq!(states(8, 8, true).unwrap().len(), 256);
    }
    #[test]
    fn mask_fast_and_general_are_equivalent() {
        for a in [
            vec![0, 63],
            vec![64, 127],
            vec![0, 128, 200],
            vec![-1],
            vec![1, 1],
        ] {
            for b in [vec![0, 63, 64, 127, 128, 200], vec![-1], vec![1, 2]] {
                assert_eq!(
                    StateBitmask::new(&a).subset_of(&StateBitmask::new(&b)),
                    a.iter().all(|x| b.contains(x))
                );
                assert_eq!(
                    StateBitmask::new(&a).disjoint(&StateBitmask::new(&b)),
                    a.iter().all(|x| !b.contains(x))
                );
            }
        }
    }
    #[test]
    fn unpack_states_rejects_fractional_and_invalid_area_indices() {
        assert_eq!(
            unpack_states(&[f64::NAN, 0.0, 1.0], &[1.0, 2.0]).unwrap(),
            vec![vec![-1], vec![0, 1]]
        );
        assert_eq!(
            unpack_states(&[-1.0, 0.0], &[1.0, 1.0]).unwrap(),
            vec![vec![-1], vec![0]]
        );
        for bad_flat in [
            &[0.5][..],
            &[-0.1][..],
            &[-2.0][..],
            &[f64::INFINITY][..],
            &[-1.0, 0.0][..],
            &[f64::NAN, 0.0][..],
        ] {
            let len = bad_flat.len() as f64;
            assert!(unpack_states(bad_flat, &[len]).is_err());
        }
        assert!(unpack_states(&[0.0], &[1.5]).is_err());
    }
}
