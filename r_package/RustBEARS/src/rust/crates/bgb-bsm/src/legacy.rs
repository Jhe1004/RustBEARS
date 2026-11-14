use crate::Event;

pub struct LegacyResult {
    pub attempts: usize,
    pub accepted: bool,
    pub events: Vec<Event>,
}

/// Run the native attempts before the final upstream R attempt. The caller
/// owns R's RNG and supplies its exact waiting-time and sample primitives.
pub fn reject(
    q: &[f64],
    n: usize,
    start: usize,
    end: usize,
    length: f64,
    max_attempts: usize,
    mut wait: impl FnMut(f64) -> f64,
    mut sample: impl FnMut(&[f64], &mut [f64], &mut [i32]) -> usize,
) -> Result<LegacyResult, String> {
    if n == 0
        || Some(q.len()) != n.checked_mul(n)
        || start >= n
        || end >= n
        || !length.is_finite()
        || length <= 0.0
    {
        return Err("invalid branch dimensions".into());
    }
    let rates: Vec<_> = (0..n).map(|i| -q[i * (n + 1)]).collect();
    for (i, rate) in rates.iter().copied().enumerate() {
        if !rate.is_finite() || rate < 0.0 || (rate > 0.0 && !rate.recip().is_finite()) {
            return Err("invalid departure rate".into());
        }
        let mut total = 0.0;
        for j in 0..n {
            let value = q[i + n * j];
            if !value.is_finite() {
                return Err("non-finite rate matrix".into());
            }
            total += value.max(0.0);
        }
        if rate > 0.0 && (!(total > 0.0) || !total.is_finite()) {
            return Err("invalid jump weights".into());
        }
    }

    const MAX_BRANCH_EVENTS: usize = 10_000;
    let mut probabilities = vec![0.0; n];
    let mut sorted = vec![0.0; n];
    let mut permutation = vec![0; n];
    let mut events = Vec::new();
    for attempts in 1..=max_attempts {
        events.clear();
        let mut state = start;
        let mut time = 0.0;
        let mut completed = false;
        loop {
            let rate = rates[state];
            if rate == 0.0 {
                completed = true;
                break;
            }
            time += wait(rate);
            if time > length {
                completed = true;
                break;
            }
            if events.len() >= MAX_BRANCH_EVENTS {
                break;
            }
            for (j, value) in probabilities.iter_mut().enumerate() {
                *value = q[state + n * j].max(0.0);
            }
            let next = sample(&probabilities, &mut sorted, &mut permutation);
            events.push(Event {
                from: state,
                to: next,
                time,
            });
            state = next;
        }
        if !completed {
            return Err(
                "rejection sampling exceeded 10000 transitions on a single branch; use uniformization"
                    .into(),
            );
        }
        if state == end {
            return Ok(LegacyResult {
                attempts,
                accepted: true,
                events,
            });
        }
    }
    Ok(LegacyResult {
        attempts: max_attempts,
        accepted: false,
        events: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reject_refuses_truncated_high_rate_branch_instead_of_accepting() {
        // Column-major 2x2 symmetric Q matrix with departure rate 20_000 on length 1.0:
        // deterministic waiting step 1 / 20_000 = 5e-5 reaches only t = 0.50005 after 10_001 draws.
        let r = 20_000.0;
        let q = [-r, r, r, -r];
        let res = reject(
            &q,
            2,
            0,
            0,
            1.0,
            10,
            |rate| 1.0 / rate,
            |weights, _, _| usize::from(weights[1] > weights[0]),
        );
        assert!(res.is_err(), "expected truncated high-rate rejection branch to fail");
    }
}
