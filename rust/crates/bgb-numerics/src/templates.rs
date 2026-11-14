use bgb_core::{PrecisionMode, StateBitmask};
#[derive(Clone)]
enum Transition {
    Swap(usize, usize),
    Expand(Vec<usize>, usize),
    Extinct(usize),
}
/// Legal edges in exactly the historical source/destination enumeration order.
pub struct CompiledQTemplate {
    pub states: Vec<Vec<i32>>,
    pub areas: usize,
    edges: Vec<(usize, usize, Transition)>,
    row_ptr: Vec<usize>,
    col_idx: Vec<usize>,
    edge_slots: Vec<usize>,
    diag_slots: Vec<usize>,
}
impl CompiledQTemplate {
    pub fn new(states: &[Vec<i32>], areas: usize) -> Result<Self, String> {
        if areas == 0
            || states
                .iter()
                .any(|s| s.is_empty() || s.iter().any(|&x| x >= areas as i32))
        {
            return Err("invalid state index".into());
        }
        let masks: Vec<_> = states.iter().map(|s| StateBitmask::new(s)).collect();
        let mut edges = Vec::new();
        for (i, a) in states.iter().enumerate() {
            for (j, b) in states.iter().enumerate() {
                let kind = if a.len() == 1 && b.len() == 1 && a[0] >= 0 && b[0] >= 0 && a[0] != b[0]
                {
                    Some(Transition::Swap(a[0] as usize, b[0] as usize))
                } else if a.len() + 1 == b.len() && masks[i].subset_of(&masks[j]) {
                    if a.iter().any(|&x| x < 0) {
                        return Err("invalid empty-range expansion".into());
                    }
                    let extra = *b.iter().find(|x| !a.contains(x)).ok_or("missing area")?;
                    Some(Transition::Expand(
                        a.iter().map(|&x| x as usize).collect(),
                        extra as usize,
                    ))
                } else if b[0] < 0 && a[0] >= 0 && a.len() == 1 {
                    Some(Transition::Extinct(a[0] as usize))
                } else if a.len() == b.len() + 1 && masks[j].subset_of(&masks[i]) {
                    Some(Transition::Extinct(
                        *a.iter().find(|x| !b.contains(x)).ok_or("missing area")? as usize,
                    ))
                } else {
                    None
                };
                if let Some(kind) = kind {
                    edges.push((i, j, kind));
                }
            }
        }
        let n = states.len();
        let mut row_ptr = vec![0; n + 1];
        let mut col_idx = Vec::with_capacity(edges.len() + n);
        let mut edge_slots = Vec::with_capacity(edges.len());
        let mut diag_slots = vec![0; n];
        let mut k = 0;
        for i in 0..n {
            for j in 0..n {
                let edge = k < edges.len() && edges[k].0 == i && edges[k].1 == j;
                if edge || i == j {
                    let slot = col_idx.len();
                    col_idx.push(j);
                    if i == j {
                        diag_slots[i] = slot;
                    }
                    if edge {
                        edge_slots.push(slot);
                        k += 1;
                    }
                }
            }
            row_ptr[i + 1] = col_idx.len();
        }
        Ok(Self {
            states: states.to_vec(),
            areas,
            edges,
            row_ptr,
            col_idx,
            edge_slots,
            diag_slots,
        })
    }
    fn value(
        &self,
        kind: &Transition,
        d: &[f64],
        e: &[f64],
        amat: &[f64],
        mode: PrecisionMode,
    ) -> f64 {
        let na = self.areas;
        match kind {
            Transition::Swap(a, b) => mode.round(amat[a + na * b]),
            Transition::Extinct(a) => mode.round(e[*a]),
            Transition::Expand(from, to) => {
                let mut x = 0.0;
                for a in from {
                    x = mode.round(x + mode.round(d[a + na * to]));
                }
                x
            }
        }
    }
    /// Parameter-only sparse refresh. Topology is already sorted and includes
    /// every potential diagonal; filtered edges remain explicit zero slots.
    pub fn refresh_csr(
        &self,
        d: &[f64],
        e: &[f64],
        amat: &[f64],
        normalize: bool,
        min: f64,
        mode: PrecisionMode,
    ) -> Result<super::CsrMatrix, String> {
        let na = self.areas;
        let n = self.states.len();
        if d.len() != na * na || amat.len() != na * na || e.len() < na {
            return Err("invalid DEmat dimensions".into());
        }
        let threshold = mode.round(min);
        let mut values = vec![0.0; self.col_idx.len()];
        let mut full = vec![0.0; n];
        let mut last_source = None;
        let mut before_last = 0.0;
        for (k, (i, j, kind)) in self.edges.iter().enumerate() {
            let x = self.value(kind, d, e, amat, mode);
            if i != j && x >= threshold {
                values[self.edge_slots[k]] = x;
                last_source = Some(*i);
                before_last = full[*i];
                full[*i] = mode.round(full[*i] + x);
            }
        }
        let mut appended = false;
        if normalize {
            for i in 0..n {
                let sum = if mode == PrecisionMode::LegacyF32 && !appended && last_source == Some(i)
                {
                    before_last
                } else {
                    full[i]
                };
                if sum >= threshold {
                    values[self.diag_slots[i]] = -sum;
                    appended = true;
                }
            }
        }
        Ok(super::CsrMatrix {
            n,
            row_ptr: self.row_ptr.clone(),
            col_idx: self.col_idx.clone(),
            values,
        })
    }
    pub fn refresh(
        &self,
        d: &[f64],
        e: &[f64],
        amat: &[f64],
        normalize: bool,
        coo: bool,
        min: f64,
        mode: PrecisionMode,
    ) -> Result<Vec<Vec<f64>>, String> {
        let na = self.areas;
        let n = self.states.len();
        if d.len() != na * na || amat.len() != na * na || e.len() < na {
            return Err("invalid DEmat dimensions".into());
        }
        let threshold = mode.round(if coo { min } else { 1e-15 });
        let mut q = if coo { Vec::new() } else { vec![0.0; n * n] };
        let mut ia = Vec::new();
        let mut ja = Vec::new();
        let mut vals = Vec::new();
        for (i, j, kind) in &self.edges {
            let x = self.value(kind, d, e, amat, mode);
            if coo {
                if x >= threshold && i != j {
                    ia.push((j + 1) as f64);
                    ja.push((i + 1) as f64);
                    vals.push(x);
                }
            } else {
                q[i + n * j] = x;
            }
        }
        if coo {
            // Upstream initially excludes the final off-diagonal entry. After
            // the first diagonal append the excluded entry is that diagonal.
            let mut full = vec![0.0; n];
            let mut initial = vec![0.0; n];
            let count = vals.len();
            for k in 0..count {
                let i = ja[k] as usize - 1;
                full[i] = mode.round(full[i] + vals[k]);
                if k + 1 < count {
                    initial[i] = mode.round(initial[i] + vals[k]);
                }
            }
            let mut appended = false;
            if normalize {
                for i in 0..n {
                    let sum = if mode == PrecisionMode::LegacyF32 && !appended {
                        initial[i]
                    } else {
                        full[i]
                    };
                    if sum >= threshold {
                        ia.push((i + 1) as f64);
                        ja.push((i + 1) as f64);
                        vals.push(-sum);
                        appended = true;
                    }
                }
            }
            Ok(vec![ia, ja, vals])
        } else {
            if normalize {
                for i in 0..n {
                    let mut sum = 0.0;
                    for j in 0..n {
                        if i != j {
                            sum = mode.round(sum + q[i + n * j]);
                        }
                    }
                    q[i + n * i] = if mode == PrecisionMode::LegacyF32 && sum < threshold {
                        0.0
                    } else {
                        -sum
                    };
                }
            }
            Ok(vec![q])
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn template_matches_unoptimized_port_at_rounding_and_filter_boundaries() {
        for na in 1..=5 {
            let states = bgb_core::states(na, na, true).unwrap();
            let t = CompiledQTemplate::new(&states, na).unwrap();
            for seed in 0..4 {
                let d: Vec<_> = (0..na * na)
                    .map(|i| {
                        if seed == 0 {
                            0.0
                        } else {
                            (i as f64 + 0.331) * (seed as f64 + 0.71) * 1e-3
                        }
                    })
                    .collect();
                let e: Vec<_> = (0..na).map(|i| (i as f64 + 1.31) * 1e-4).collect();
                for mode in [PrecisionMode::LegacyF32, PrecisionMode::StrictF64] {
                    for coo in [false, true] {
                        for normalize in [false, true] {
                            let a = t.refresh(&d, &e, &d, normalize, coo, 1e-26, mode).unwrap();
                            let b =
                                crate::demat(&states, na, &d, &e, &d, normalize, coo, 1e-26, mode)
                                    .unwrap();
                            assert_eq!(
                                a.iter().flatten().map(|x| x.to_bits()).collect::<Vec<_>>(),
                                b.iter().flatten().map(|x| x.to_bits()).collect::<Vec<_>>()
                            );
                            if coo {
                                let csr =
                                    t.refresh_csr(&d, &e, &d, normalize, 1e-26, mode).unwrap();
                                let other = crate::CsrMatrix::from_triplets(
                                    states.len(),
                                    &a[0],
                                    &a[1],
                                    &a[2],
                                    true,
                                )
                                .unwrap();
                                assert_eq!(csr.dense().a, other.dense().a);
                            }
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn template_has_general_state_fallback_above_128_areas() {
        let na = 129;
        let states = vec![vec![-1], vec![0], vec![128], vec![0, 128]];
        let d = vec![0.031; na * na];
        let e = vec![0.023; na];
        let tpl = CompiledQTemplate::new(&states, na).unwrap();
        for mode in [PrecisionMode::LegacyF32, PrecisionMode::StrictF64] {
            let actual = tpl.refresh(&d, &e, &d, true, false, 1e-26, mode).unwrap();
            let old = crate::demat(&states, na, &d, &e, &d, true, false, 1e-26, mode).unwrap();
            assert_eq!(actual, old);
        }
    }
}
