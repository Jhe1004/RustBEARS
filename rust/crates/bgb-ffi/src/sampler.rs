//! Owned samplers; wrappers serialize native use and final release.
use super::{buffers, Reply};
use bgb_bsm::{LayerInput, PreparedTreeSampler, SegmentInput, TreeInput};
use std::ffi::{c_char, CString};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::Mutex;
pub struct SamplerHandle {
    sampler: Mutex<Option<PreparedTreeSampler>>,
    error: Option<CString>,
}
fn idx(x: f64) -> Result<usize, String> {
    if !x.is_finite() || x < 0.0 || x.fract() != 0.0 || x > 9_007_199_254_740_991.0 {
        Err("invalid sampler integer".into())
    } else {
        Ok(x as usize)
    }
}
fn prepare(a: &[&[f64]], root: &[f64], mode: u32) -> Result<PreparedTreeSampler, String> {
    if a.len() < 14
        || a[0].len() != 2
        || a[4].len() % 3 != 0
        || a[9].len() % 4 != 0
        || (a.len() - 10) % 4 != 0
    {
        return Err("invalid prepared sampler buffers".into());
    }
    let states = bgb_core::unpack_states(a[1], a[2])?;
    let steps = a[4]
        .chunks_exact(3)
        .map(|x| Ok([idx(x[0])?, idx(x[1])?, idx(x[2])?]))
        .collect::<Result<Vec<_>, String>>()?;
    let node_layers = a[6]
        .iter()
        .map(|x| idx(*x))
        .collect::<Result<Vec<_>, _>>()?;
    let segments = a[9]
        .chunks_exact(4)
        .map(|x| {
            Ok(SegmentInput {
                node: idx(x[0])?,
                layer: idx(x[1])?,
                length: x[2],
                start_age: x[3],
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let layers = a[10..]
        .chunks_exact(4)
        .map(|x| LayerInput {
            q: x[0],
            allowed: x[1],
            events: x[2],
            dispersal: x[3],
        })
        .collect::<Vec<_>>();
    PreparedTreeSampler::prepare(
        TreeInput {
            states: &states,
            tips: a[3],
            tip_count: idx(a[0][0])?,
            root: idx(a[0][1])?,
            steps: &steps,
            ages: a[5],
            layers_at_node: &node_layers,
            hooks: a[7],
            fixes: a[8],
            segments: &segments,
            layers: &layers,
        },
        root,
        mode,
    )
}
/// # Safety
/// Buffers remain readable for the declared lengths through this call. The
/// returned handle is uniquely owned and released exactly once.
/// `conditioning` is 0 (observed), 1 (tip survival), 2 (complete process),
/// or 4 (tip survival with an already conditioned root distribution).
/// Mode 3 is private to the legacy adapter and is rejected here.
#[no_mangle]
pub unsafe extern "C" fn bgb_sampler_prepare(
    ptrs: *const *const f64,
    lens: *const usize,
    n: usize,
    root_weights: *const f64,
    root_len: usize,
    conditioning: u32,
) -> *mut SamplerHandle {
    let r = catch_unwind(AssertUnwindSafe(|| {
        let a = unsafe { buffers(ptrs, lens, n) };
        let root = if root_len == 0 {
            &[][..]
        } else {
            if root_weights.is_null() {
                return Err("null root weights".to_string());
            }
            unsafe { std::slice::from_raw_parts(root_weights, root_len) }
        };
        prepare(&a, root, conditioning)
    }));
    let (sampler, error) = match r {
        Ok(Ok(s)) => (Some(s), None),
        Ok(Err(e)) => (None, Some(CString::new(e.replace('\0', " ")).unwrap())),
        Err(_) => (
            None,
            Some(CString::new("sampler preparation panicked").unwrap()),
        ),
    };
    Box::into_raw(Box::new(SamplerHandle {
        sampler: Mutex::new(sampler),
        error,
    }))
}
/// # Safety
/// The handle is live.
#[no_mangle]
pub unsafe extern "C" fn bgb_sampler_error(p: *const SamplerHandle) -> *const c_char {
    if p.is_null() {
        return std::ptr::null();
    }
    unsafe { (*p).error.as_ref().map_or(std::ptr::null(), |x| x.as_ptr()) }
}
fn reply(call: impl FnOnce() -> Result<Vec<Vec<f64>>, String>) -> *mut Reply {
    Box::into_raw(Box::new(match catch_unwind(AssertUnwindSafe(call)) {
        Ok(Ok(buffers)) => Reply {
            buffers,
            error: None,
        },
        Ok(Err(e)) => Reply {
            buffers: Vec::new(),
            error: Some(CString::new(e.replace('\0', " ")).unwrap()),
        },
        Err(_) => Reply {
            buffers: Vec::new(),
            error: Some(CString::new("sampler operation panicked").unwrap()),
        },
    }))
}
/// # Safety
/// Handle remains live through the call. The caller releases the returned Reply.
#[no_mangle]
pub unsafe extern "C" fn bgb_sampler_draw_bsm(
    p: *mut SamplerHandle,
    start: u64,
    count: usize,
    seed: u64,
    threads: usize,
    rejection: u32,
) -> *mut Reply {
    reply(|| {
        if p.is_null() {
            return Err("null sampler handle".into());
        }
        if rejection > 1 {
            return Err("invalid rejection flag".into());
        }
        let lock = unsafe { &*p }
            .sampler
            .lock()
            .map_err(|_| "sampler lock poisoned")?;
        lock.as_ref().ok_or("invalid sampler handle")?.draw_bsm(
            start,
            count,
            seed,
            threads,
            rejection != 0,
        )
    })
}
/// # Safety
/// Handle remains live through the call. The caller releases the returned Reply.
/// Output arrays are row-major replicate-by-tip and replicate-by-node indices.
#[no_mangle]
pub unsafe extern "C" fn bgb_sampler_simulate_tips(
    p: *mut SamplerHandle,
    start: u64,
    count: usize,
    seed: u64,
    threads: usize,
) -> *mut Reply {
    reply(|| {
        if p.is_null() {
            return Err("null sampler handle".into());
        }
        let lock = unsafe { &*p }
            .sampler
            .lock()
            .map_err(|_| "sampler lock poisoned")?;
        lock.as_ref()
            .ok_or("invalid sampler handle")?
            .simulate_tips(start, count, seed, threads)
    })
}
/// # Safety
/// Handle is uniquely owned and no operation is in progress. Wrappers clear
/// their pointer after release so explicit close and finalizers cannot double free.
#[no_mangle]
pub unsafe extern "C" fn bgb_sampler_release(p: *mut SamplerHandle) {
    if !p.is_null() {
        drop(unsafe { Box::from_raw(p) });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn buffers_fixture() -> Vec<Vec<f64>> {
        vec![
            vec![2.0, 2.0],
            vec![-1.0, 0.0, 1.0],
            vec![1.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
            vec![2.0, 0.0, 1.0],
            vec![0.0, 0.0, 1.0],
            vec![0.0; 3],
            vec![0.0; 3],
            vec![1.0; 9],
            vec![0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.1, 1.0, 0.0, -0.1, 0.0, 0.0, 0.0, -1.0],
            vec![1.0; 3],
            vec![1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 1.0],
            vec![1.0; 4],
        ]
    }

    fn owned_handle(mode: u32) -> *mut SamplerHandle {
        let a = buffers_fixture();
        let ptrs: Vec<_> = a.iter().map(|x| x.as_ptr()).collect();
        let lens: Vec<_> = a.iter().map(Vec::len).collect();
        let roots = [0.0, 0.4, 0.6];
        // All input arrays are dropped before the caller uses the returned handle.
        unsafe {
            bgb_sampler_prepare(
                ptrs.as_ptr(),
                lens.as_ptr(),
                a.len(),
                roots.as_ptr(),
                roots.len(),
                mode,
            )
        }
    }

    #[test]
    fn conditioned_root_is_not_weighted_by_survival_again() {
        unsafe {
            let fixed = owned_handle(4);
            assert!(bgb_sampler_error(fixed).is_null());
            let ordinary = owned_handle(1);
            assert!(bgb_sampler_error(ordinary).is_null());
            let fixed_reply = Box::from_raw(bgb_sampler_simulate_tips(fixed, 1, 20_000, 91, 2));
            let ordinary_reply =
                Box::from_raw(bgb_sampler_simulate_tips(ordinary, 1, 20_000, 91, 1));
            assert!(fixed_reply.error.is_none());
            assert!(ordinary_reply.error.is_none());
            let frequency = |r: &Reply| {
                r.buffers[1]
                    .chunks_exact(3)
                    .filter(|row| row[2] == 1.0)
                    .count() as f64
                    / 20_000.0
            };
            assert!((frequency(&fixed_reply) - 0.4).abs() < 0.015);
            let weighted =
                0.4 * (-0.2_f64).exp() / (0.4 * (-0.2_f64).exp() + 0.6 * (-2.0_f64).exp());
            assert!((frequency(&ordinary_reply) - weighted).abs() < 0.015);
            assert!(fixed_reply.buffers[0].iter().all(|x| *x != 0.0));
            assert!(fixed_reply.buffers[1]
                .chunks_exact(3)
                .all(|row| row[0] == row[2] && row[1] == row[2]));
            let misuse = Box::from_raw(bgb_sampler_draw_bsm(fixed, 1, 1, 91, 1, 0));
            assert!(misuse.error.is_some());
            bgb_sampler_release(fixed);
            bgb_sampler_release(ordinary);
        }
    }

    #[test]
    fn private_legacy_mode_is_not_exposed_by_preparation() {
        unsafe {
            let invalid = owned_handle(3);
            assert!(!bgb_sampler_error(invalid).is_null());
            bgb_sampler_release(invalid);
            bgb_sampler_release(std::ptr::null_mut());
        }
    }

    #[test]
    fn conditioned_root_batches_and_threads_keep_replicate_identity() {
        unsafe {
            let handle = owned_handle(4);
            assert!(bgb_sampler_error(handle).is_null());
            let all = Box::from_raw(bgb_sampler_simulate_tips(handle, 1, 101, 12345, 1));
            let first = Box::from_raw(bgb_sampler_simulate_tips(handle, 1, 37, 12345, 4));
            let second = Box::from_raw(bgb_sampler_simulate_tips(handle, 38, 64, 12345, 2));
            assert!(all.error.is_none() && first.error.is_none() && second.error.is_none());
            for i in 0..2 {
                let joined = [&first.buffers[i][..], &second.buffers[i][..]].concat();
                assert_eq!(all.buffers[i], joined);
            }
            bgb_sampler_release(handle);
        }
    }
}
