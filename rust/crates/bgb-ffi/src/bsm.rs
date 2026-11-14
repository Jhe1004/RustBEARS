use super::Reply;
use bgb_bsm::{BridgeJob, Event, Generator, Xoshiro256PlusPlus};
use std::ffi::CString;
use std::panic::{catch_unwind, AssertUnwindSafe};

fn packed(events: Vec<Event>) -> Vec<f64> {
    events
        .into_iter()
        .flat_map(|e| [(e.from + 1) as f64, (e.to + 1) as f64, e.time])
        .collect()
}

fn key(x: f64) -> Result<u64, String> {
    if !(0.0..=9_007_199_254_740_991.0).contains(&x) || x.fract() != 0.0 {
        return Err("BSM key must be a nonnegative, exactly representable integer".into());
    }
    Ok(x as u64)
}

pub(super) fn bridges(a: &[&[f64]]) -> Result<Vec<Vec<f64>>, String> {
    if a.len() < 2 || a[0].len() != 2 || a[1].len() % 8 != 0 {
        return Err("invalid BSM batch buffers".into());
    }
    let seed = key(a[0][0])?;
    let threads = key(a[0][1])? as usize;
    if !(1..=256).contains(&threads) {
        return Err("BSM threads must be between 1 and 256".into());
    }
    let generators: Result<Vec<_>, _> = a[2..]
        .iter()
        .map(|q| {
            let n = (q.len() as f64).sqrt() as usize;
            Generator::new(q, n)
        })
        .collect();
    let jobs: Result<Vec<_>, String> = a[1]
        .chunks_exact(8)
        .map(|j| {
            Ok(BridgeJob {
                generator: key(j[0])?.checked_sub(1).ok_or("zero generator index")? as usize,
                from: key(j[1])?.checked_sub(1).ok_or("zero starting state")? as usize,
                to: key(j[2])?.checked_sub(1).ok_or("zero ending state")? as usize,
                duration: j[3],
                map: key(j[4])?,
                stratum: key(j[5])?,
                step: key(j[6])?,
                node: key(j[7])?,
            })
        })
        .collect();
    Ok(bgb_bsm::sample_batch(&generators?, &jobs?, seed, threads)?
        .into_iter()
        .map(packed)
        .collect())
}

pub(super) fn uniforms(a: &[&[f64]]) -> Result<Vec<Vec<f64>>, String> {
    if a.len() != 1 || a[0].len() % 4 != 0 {
        return Err("invalid BSM stream buffers".into());
    }
    let values = a[0]
        .chunks_exact(4)
        .map(|k| {
            if k[1] == 0.0 {
                return Err("map index must be positive".into());
            }
            let mut rng =
                Xoshiro256PlusPlus::from_key(key(k[0])?, key(k[1])?, key(k[2])?, key(k[3])?);
            Ok(rng.uniform())
        })
        .collect::<Result<_, String>>()?;
    Ok(vec![values])
}

/// # Safety
/// The matrix is column-major and lives through the call. The callbacks
/// perform arithmetic/R RNG only. They cannot allocate R objects or longjmp.
#[no_mangle]
pub unsafe extern "C" fn bgb_bsm_reject(
    q: *const f64,
    n: usize,
    start: usize,
    end: usize,
    duration: f64,
    attempts: usize,
    wait: extern "C" fn(f64) -> f64,
    pick: extern "C" fn(*const f64, *mut f64, *mut i32, usize) -> usize,
) -> *mut Reply {
    let call = catch_unwind(AssertUnwindSafe(|| {
        let count = n.checked_mul(n).ok_or("legacy matrix too large")?;
        let q = unsafe { std::slice::from_raw_parts(q, count) };
        let result = bgb_bsm::reject(
            q,
            n,
            start,
            end,
            duration,
            attempts,
            |rate| wait(rate),
            |weights, work, perm| {
                pick(
                    weights.as_ptr(),
                    work.as_mut_ptr(),
                    perm.as_mut_ptr(),
                    weights.len(),
                )
            },
        )?;
        Ok::<_, String>(vec![
            vec![result.attempts as f64, result.accepted as u8 as f64],
            packed(result.events),
        ])
    }));
    let reply = match call {
        Ok(Ok(buffers)) => Reply {
            buffers,
            error: None,
        },
        Ok(Err(error)) => Reply {
            buffers: Vec::new(),
            error: CString::new(error).ok(),
        },
        Err(_) => Reply {
            buffers: Vec::new(),
            error: CString::new("legacy BSM kernel panicked").ok(),
        },
    };
    Box::into_raw(Box::new(reply))
}

pub(super) fn whole_tree(a: &[&[f64]]) -> Result<Vec<Vec<f64>>, String> {
    if a.len() < 10
        || a[0].len() != 6
        || a[4].len() % 3 != 0
        || a[9].len() % 4 != 0
        || (a.len() - 10) % 4 != 0
    {
        return Err("invalid whole-tree BSM buffers".into());
    }
    let flags = a[0]
        .iter()
        .map(|&v| key(v).map(|x| x as usize))
        .collect::<Result<Vec<_>, _>>()?;
    let [tips, root, maps, _, threads, rejection] = flags[..] else {
        unreachable!()
    };
    if threads == 0 || threads > 256 || maps == 0 || maps > 1_000_000 || tips < 2 {
        return Err("invalid BSM size or threads".into());
    }
    let states = bgb_core::unpack_states(a[1], a[2])?;
    let indices = |v: f64| key(v).map(|x| x as usize);
    let steps = a[4]
        .chunks_exact(3)
        .map(|x| Ok([indices(x[0])?, indices(x[1])?, indices(x[2])?]))
        .collect::<Result<Vec<_>, String>>()?;
    let node_layers = a[6]
        .iter()
        .map(|&v| indices(v))
        .collect::<Result<Vec<_>, _>>()?;
    let segments = a[9]
        .chunks_exact(4)
        .map(|x| {
            Ok(bgb_bsm::SegmentInput {
                node: indices(x[0])?,
                layer: indices(x[1])?,
                length: x[2],
                start_age: x[3],
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let layers = a[10..]
        .chunks_exact(4)
        .map(|l| bgb_bsm::LayerInput {
            q: l[0],
            allowed: l[1],
            events: l[2],
            dispersal: l[3],
        })
        .collect::<Vec<_>>();
    bgb_bsm::sample_tree(
        bgb_bsm::TreeInput {
            states: &states,
            tips: a[3],
            tip_count: tips,
            root,
            steps: &steps,
            ages: a[5],
            layers_at_node: &node_layers,
            hooks: a[7],
            fixes: a[8],
            segments: &segments,
            layers: &layers,
        },
        maps,
        key(a[0][3])?,
        threads,
        rejection != 0,
    )
}
