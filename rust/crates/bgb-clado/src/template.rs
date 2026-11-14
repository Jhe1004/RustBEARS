use super::{Events, Model};
use bgb_core::StateBitmask;
/// All potential events, in the upstream y/s/v/j and mirror ordering.
/// Runtime weights are refreshed without re-enumerating state combinations.
pub struct CompiledCladoTemplate {
    pub states: Vec<Vec<i32>>,
    pub areas: usize,
    events: Vec<(usize, usize, usize, usize, usize)>,
}
impl CompiledCladoTemplate {
    pub fn new(m: &Model<'_>) -> Result<Self, String> {
        m.validate()?;
        let na = m.na;
        let states = m.states;
        let masks: Vec<_> = states.iter().map(|s| StateBitmask::new(s)).collect();
        let group_count = states.iter().map(|s| s.len()).max().unwrap_or(0);
        let mut events = Vec::new();
        for a in 0..states.len() {
            events.push((a, a, a, 3, a));
        }
        let mirror = |out: &mut Vec<(usize, usize, usize, usize, usize)>,
                      items: Vec<(usize, usize, usize, usize, usize)>| {
            out.extend(items.iter().copied());
            out.extend(items.into_iter().map(|(a, l, r, k, c)| (a, r, l, k, c)));
        };
        let mut subset = Vec::new();
        let mut equal = Vec::new();
        let mut unequal = Vec::new();
        let mut jumps = Vec::new();
        for (a, s) in states.iter().enumerate() {
            for size in 1..group_count {
                for (r, t) in states.iter().enumerate().filter(|(_, t)| t.len() == size) {
                    if s.len() >= 2 && t != s && masks[r].subset_of(&masks[a]) {
                        subset.push((a, a, r, 0, r));
                    }
                    if s.len() >= 2 && masks[r].subset_of(&masks[a]) {
                        let left: Vec<_> = s.iter().filter(|x| !t.contains(x)).copied().collect();
                        if !left.is_empty() {
                            if let Some(l) = states.iter().position(|x| *x == left) {
                                let event = (a, l, r, 1, r);
                                if left.len() == t.len() {
                                    equal.push(event);
                                } else {
                                    unequal.push(event);
                                }
                            }
                        }
                    }
                    if s.len() <= 2000 && masks[r].disjoint(&masks[a]) {
                        jumps.push((a, a, r, 2, r));
                    }
                }
            }
        }
        mirror(&mut events, subset);
        events.extend(equal);
        mirror(&mut events, unequal);
        mirror(&mut events, jumps);
        Ok(Self {
            states: m.states.to_vec(),
            areas: na,
            events,
        })
    }
    pub fn refresh(&self, m: &Model<'_>) -> Result<Events, String> {
        m.validate()?;
        if self.states != m.states || self.areas != m.na {
            return Err("cladogenesis template mismatch".into());
        }
        let cut = m.mode.round(1e-30);
        let mut out = Vec::with_capacity(self.events.len());
        for &(a, l, r, kind, child) in &self.events {
            let size = m.states[a].len();
            let k = m.states[child].len();
            if m.rate(kind).abs() <= cut || k > m.maxmin[size - 1] as usize {
                continue;
            }
            let w = if kind == 2 {
                m.mode.round(
                    m.rate(kind)
                        * m.prior(kind, size, k)
                        * m.distance(&m.states[a], &m.states[child]),
                )
            } else {
                m.mode.round(m.rate(kind) * m.prior(kind, size, k))
            };
            if w.abs() > cut {
                out.push((a, l, r, w));
            }
        }
        Ok(out)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use bgb_core::PrecisionMode;
    #[test]
    fn template_preserves_mirrors_prior_sizes_and_mixed_precision() {
        for na in 2..=5 {
            let s = bgb_core::states(na, na, false).unwrap();
            let d: Vec<_> = (0..na * na).map(|i| (i as f64 + 0.314) * 0.11).collect();
            let prior: Vec<_> = (0..na * na).map(|i| (i as f64 + 0.713) * 0.017).collect();
            for max in [1.0, 2.0, na as f64] {
                let maxmin = vec![max; na];
                for params in [
                    [0.0, 1.0, 0.0, 1.0],
                    [0.3, 0.7, 0.4, 0.1],
                    [1e-32, 0.0, 1.1, 0.0],
                ] {
                    for mode in [PrecisionMode::LegacyF32, PrecisionMode::StrictF64] {
                        let m = Model {
                            states: &s,
                            parameters: &params,
                            d: &d,
                            priors: [&prior, &prior, &prior, &prior],
                            maxmin: &maxmin,
                            na,
                            mode,
                        };
                        let a = m.faster().unwrap();
                        let b = CompiledCladoTemplate::new(&m).unwrap().refresh(&m).unwrap();
                        assert_eq!(a.len(), b.len());
                        for (x, y) in a.iter().zip(&b) {
                            assert_eq!(
                                (x.0, x.1, x.2, x.3.to_bits()),
                                (y.0, y.1, y.2, y.3.to_bits())
                            );
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn cladogenesis_template_has_general_state_fallback_above_128_areas() {
        let na = 129;
        let states = vec![vec![0], vec![128], vec![0, 128]];
        let d = vec![0.031; na * na];
        let p = vec![0.23; na * na];
        let limit = vec![1.0; na];
        let rates = [0.3, 0.7, 0.2, 0.9];
        for mode in [PrecisionMode::LegacyF32, PrecisionMode::StrictF64] {
            let m = Model {
                states: &states,
                parameters: &rates,
                d: &d,
                priors: [&p, &p, &p, &p],
                maxmin: &limit,
                na,
                mode,
            };
            let tpl = CompiledCladoTemplate::new(&m).unwrap();
            assert_eq!(m.faster().unwrap(), tpl.refresh(&m).unwrap());
        }
    }
}
