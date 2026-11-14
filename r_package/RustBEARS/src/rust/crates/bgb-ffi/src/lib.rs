//! R-free C ABI. No R longjmp or R allocation crosses a live Rust stack.
// Indexed loops intentionally preserve the reference reduction order.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]
mod bsm;
use bgb_core::{unpack_states, PrecisionMode};
use bgb_numerics::{ActionPlan, Matrix};
use std::cell::RefCell;
use std::ffi::{c_char, CString};
use std::panic::{catch_unwind, AssertUnwindSafe};
thread_local! {static ACTION_CACHE:RefCell<Vec<(Vec<f64>,usize,ActionPlan)>>=const{RefCell::new(Vec::new())};}
thread_local! {static SPARSE_ACTION_CACHE:RefCell<Vec<(Vec<u64>,ActionPlan)>>=const{RefCell::new(Vec::new())};}
thread_local! {static Q_TEMPLATES:RefCell<Vec<bgb_numerics::CompiledQTemplate>>=const{RefCell::new(Vec::new())};}
thread_local! {static C_TEMPLATES:RefCell<Vec<bgb_clado::CompiledCladoTemplate>>=const{RefCell::new(Vec::new())};}
pub struct Reply {
    buffers: Vec<Vec<f64>>,
    error: Option<CString>,
}
static CONTEXT_CREATED: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
static CONTEXT_DROPPED: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
static STRAT_CREATED: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
static STRAT_DROPPED: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
#[no_mangle]
pub extern "C" fn bgb_m0_abi_version() -> u32 {
    1
}
fn need(a: &[&[f64]], i: usize) -> Result<f64, String> {
    a.get(i)
        .and_then(|x| x.first())
        .copied()
        .ok_or_else(|| format!("missing scalar {i}"))
}
fn scalar_usize(x: f64, label: &str) -> Result<usize, String> {
    if !x.is_finite() || x < 0.0 || x > (i32::MAX as f64) || x.fract() != 0.0 {
        return Err(format!("invalid integer {label}"));
    }
    Ok(x as usize)
}
fn need_usize(a: &[&[f64]], i: usize, label: &str) -> Result<usize, String> {
    scalar_usize(need(a, i)?, label)
}
fn model<'a>(
    a: &'a [&'a [f64]],
    s: &'a [Vec<i32>],
    mode: PrecisionMode,
) -> Result<bgb_clado::Model<'a>, String> {
    if a.len() < 11 {
        return Err("missing clado buffers".into());
    }
    let m = bgb_clado::Model {
        states: s,
        parameters: a[2],
        d: a[3],
        priors: [a[4], a[5], a[6], a[7]],
        maxmin: a[8],
        na: need_usize(a, 9, "area count")?,
        mode,
    };
    m.validate()?;
    Ok(m)
}
fn run(op: u32, a: &[&[f64]], mode: PrecisionMode, fused: bool) -> Result<Vec<Vec<f64>>, String> {
    match op {
        1 => Ok(bgb_core::states(
            need_usize(a, 0, "area count")?,
            need_usize(a, 1, "max range size")?,
            need(a, 2)? != 0.0,
        )?
        .into_iter()
        .map(|s| s.into_iter().map(|x| x as f64).collect())
        .collect()),
        2 => Ok(vec![bgb_core::combinations(
            need_usize(a, 0, "n")?,
            need_usize(a, 1, "k")?,
            need_usize(a, 2, "limit")?,
        )?
        .into_iter()
        .flatten()
        .map(|x| x as f64)
        .collect()]),
        3 => {
            let s = unpack_states(a[0], a[1])?;
            if a[5].len() != 4 {
                return Err("missing DEmat options".into());
            }
            let na = scalar_usize(a[5][0], "area count")?;
            Q_TEMPLATES.with(|cell| {
                let mut cache = cell.borrow_mut();
                let index =
                    if let Some(i) = cache.iter().position(|x| x.states == s && x.areas == na) {
                        i
                    } else {
                        if cache.len() >= 4 {
                            cache.remove(0);
                        }
                        cache.push(bgb_numerics::CompiledQTemplate::new(&s, na)?);
                        cache.len() - 1
                    };
                cache[index].refresh(
                    a[2],
                    a[3],
                    a[4],
                    a[5][1] != 0.0,
                    a[5][2] != 0.0,
                    a[5][3],
                    mode,
                )
            })
        }
        4 => {
            let s = unpack_states(a[0], a[1])?;
            let m = model(a, &s, mode)?;
            let events = C_TEMPLATES.with(|cell| {
                let mut cache = cell.borrow_mut();
                let i = if let Some(i) = cache.iter().position(|x| x.states == s && x.areas == m.na)
                {
                    i
                } else {
                    if cache.len() >= 4 {
                        cache.remove(0);
                    }
                    cache.push(bgb_clado::CompiledCladoTemplate::new(&m)?);
                    cache.len() - 1
                };
                cache[i].refresh(&m)
            })?;
            let mut out = vec![Vec::new(); 4];
            for (anc, i, j, w) in events {
                out[0].push(anc as f64);
                out[1].push(i as f64);
                out[2].push(j as f64);
                out[3].push(w);
            }
            Ok(out)
        }
        5 => {
            let s = unpack_states(a[0], a[1])?;
            let m = model(a, &s, mode)?;
            m.dense(a[11], a[12], a[13], need_usize(a, 14, "clado kind")? as u32)
        }
        6 => Ok(vec![bgb_clado::rowsums(
            a[0],
            a[1],
            need_usize(a, 2, "state count")?,
        )?]),
        7 => Ok(vec![bgb_clado::split(
            a[0], a[1], a[2], a[3], a[4], a[5], a[6], mode,
        )?]),
        8 => {
            if a.len() < 7
                || a[2].len() != a[3].len()
                || a[2].len() != a[4].len()
                || a[2].len() != a[5].len()
            {
                return Err("COO buffer mismatch".into());
            }
            let mut out = vec![0.0; a[6].len()];
            for k in 0..a[2].len() {
                let anc = bgb_clado::checked_coo_index(a[2][k], out.len(), "ancestor")?;
                let i = bgb_clado::checked_coo_index(a[3][k], a[0].len(), "left")?;
                let j = bgb_clado::checked_coo_index(a[4][k], a[1].len(), "right")?;
                let w = mode.round(a[5][k]);
                let cell = mode.round(a[0][i] * a[1][j] * w / a[6][anc]);
                out[anc] = mode.round(out[anc] + cell);
            }
            Ok(vec![out])
        }
        9 => {
            if a[0].is_empty() || a[1].is_empty() {
                return Err("empty convolution input".into());
            }
            let mut out = vec![0.0; a[0].len() + a[1].len() - 1];
            for (i, &x) in a[0].iter().enumerate() {
                for (j, &y) in a[1].iter().enumerate() {
                    out[i + j] += x * y;
                }
            }
            Ok(vec![out])
        }
        10 => {
            let mut out = Vec::new();
            for &y in a[1] {
                for &x in a[0] {
                    out.push(x * y)
                }
            }
            Ok(vec![out])
        }
        11 => {
            let m = Matrix::new(need_usize(a, 1, "matrix order")?, a[0].to_vec())?;
            Ok(vec![if mode == PrecisionMode::LegacyF32 {
                m.legacy_exp(need(a, 2)?, fused)?.into_vec()
            } else {
                m.strict_exp(need(a, 2)?)?.into_vec()
            }])
        }
        12 => {
            let m = Matrix::new(need_usize(a, 1, "matrix order")?, a[0].to_vec())?;
            let dim = m.dim();
            ACTION_CACHE.with(|cell| {
                let mut cache = cell.borrow_mut();
                let index = if let Some(i) = cache.iter().position(|(data, n, _)| {
                    *n == dim
                        && data
                            .iter()
                            .zip(m.as_slice())
                            .all(|(x, y)| x.to_bits() == y.to_bits())
                }) {
                    i
                } else {
                    if cache.len() >= 8 {
                        cache.remove(0);
                    }
                    let plan = ActionPlan::new(&m);
                    cache.push((m.into_vec(), dim, plan));
                    cache.len() - 1
                };
                Ok(vec![cache[index].2.apply(need(a, 2)?, a[3])?])
            })
        }
        13 => bgb_core::maxent(
            a[0],
            need_usize(a, 1, "state count")?,
            need_usize(a, 2, "trait count")?,
            a[3],
            a[4],
            need(a, 5)?,
            mode,
        ),
        14 => {
            let n = need_usize(a, 3, "matrix order")?;
            let tr = need(a, 4)? != 0.0;
            if a[0].len() != a[1].len() || a[0].len() != a[2].len() {
                return Err("COO lengths differ".into());
            }
            let (row, col) = if tr { (a[1], a[0]) } else { (a[0], a[1]) };
            let mut order: Vec<_> = (0..row.len()).collect();
            order.sort_by(|&i, &j| row[i].total_cmp(&row[j]));
            let rows: Vec<_> = order.iter().map(|&i| row[i]).collect();
            let cols: Vec<_> = order.iter().map(|&i| col[i]).collect();
            let vals: Vec<_> = order.iter().map(|&i| a[2][i]).collect();
            let mut ptr = vec![1.0];
            // R indexes 1:length(ia) with a logical vector of length nnz-1.
            // Preserve its recycling, including the possible repeated end pointer.
            if rows.len() > 1 {
                for k in 0..rows.len() {
                    let logical_index = k % (rows.len() - 1);
                    if rows[logical_index] + 1.0 == rows[logical_index + 1] {
                        ptr.push((k + 2) as f64);
                    }
                }
            }
            ptr.push((rows.len() + 1) as f64);
            if ptr.len() < n + 1 {
                for r in 1..=n {
                    if !row.contains(&(r as f64)) {
                        ptr.push(r as f64)
                    }
                }
                ptr.sort_by(|x, y| x.total_cmp(y));
                if ptr.len() != n + 1 {
                    return Err(format!("STOP ERROR in coo2crs: the output crs_ia (indexes to which 'a' values where the rows change) should have nrows+1 entries. However, the function is producing length(crs_ia)={} and (nrows+1)={}. This indicates an error in the function dealing with missing rows. Check your inputs, n, and the function.",ptr.len(),n+1));
                }
            }
            Ok(vec![ptr, cols, vals, vec![n as f64, n as f64]])
        }
        15 => {
            let m = Matrix::new(need_usize(a, 1, "matrix order")?, a[0].to_vec())?;
            let dim = m.dim();
            let times = a[2];
            let v = a[3];
            if dim.checked_mul(times.len()) != Some(v.len())
                || times.iter().chain(v).any(|x| !x.is_finite())
            {
                return Err("batch dimensions mismatch".into());
            }
            if dim <= 128 && times.len() >= 32 {
                if let Some(plan) = bgb_numerics::EigenPlan::new(&m) {
                    return Ok(vec![
                        plan.apply_batch(times, v),
                        vec![2.0, plan.condition, plan.residual],
                    ]);
                }
            }
            let plan = ActionPlan::new(&m);
            let mut out = Vec::new();
            for (j, &t) in times.iter().enumerate() {
                out.extend(plan.apply(t, &v[j * dim..(j + 1) * dim])?);
            }
            Ok(vec![out, vec![1.0]])
        }
        16 => Ok(vec![bgb_numerics::legacy_arnoldi(
            need_usize(a, 1, "matrix order")?,
            need_usize(a, 4, "Krylov basis size")?,
            need(a, 2)?,
            a[3],
            need(a, 5)?,
            need(a, 6)?,
            a[8],
            a[9],
            a[10],
            need(a, 7)? != 0.0,
            fused,
        )?]),
        17 | 18 => {
            if a.len() < 7 {
                return Err("missing sparse action arguments".into());
            }
            let order = need(a, 0)?;
            if !order.is_finite() || order.fract() != 0.0 || order < 1.0 || order > i32::MAX as f64
            {
                return Err("invalid sparse order".into());
            }
            let n = order as usize;
            let tr = need(a, 6)? != 0.0;
            if a[5].len() != n || !need(a, 4)?.is_finite() || a[5].iter().any(|x| !x.is_finite()) {
                return Err("invalid sparse action vector/time".into());
            }
            let mut key = vec![op as u64, n as u64, tr as u64];
            for v in &a[1..4] {
                key.push(v.len() as u64);
                key.extend(v.iter().map(|x| x.to_bits()));
            }
            SPARSE_ACTION_CACHE.with(|cell| {
                let mut cache = cell.borrow_mut();
                let index = if let Some(i) = cache.iter().position(|(k, _)| *k == key) {
                    i
                } else {
                    let mut rows = Vec::new();
                    let row_slice = if op == 17 {
                        let ptr = a[1];
                        if ptr.len() < n + 1
                            || a[2].len() != a[3].len()
                            || ptr[0] != 1.0
                            || ptr[n] != (a[3].len() + 1) as f64
                            || ptr[n + 1..].iter().any(|x| *x != (a[3].len() + 1) as f64)
                            || ptr.iter().any(|x| {
                                !x.is_finite()
                                    || x.fract() != 0.0
                                    || *x < 1.0
                                    || *x > (a[3].len() + 1) as f64
                            })
                            || ptr.windows(2).any(|w| w[0] > w[1])
                        {
                            return Err("invalid CSR row pointers".into());
                        }
                        rows.reserve(a[3].len());
                        for i in 0..n {
                            for _ in ptr[i] as usize..ptr[i + 1] as usize {
                                rows.push((i + 1) as f64);
                            }
                        }
                        rows.as_slice()
                    } else {
                        a[1]
                    };
                    let plan = ActionPlan::from_triplets(n, row_slice, a[2], a[3], tr)?;
                    if cache.len() >= 8 {
                        cache.remove(0);
                    }
                    cache.push((key, plan));
                    cache.len() - 1
                };
                Ok(vec![cache[index].1.apply(need(a, 4)?, a[5])?])
            })
        }
        19 => {
            let c = CONTEXT_CREATED.load(std::sync::atomic::Ordering::Relaxed);
            let d = CONTEXT_DROPPED.load(std::sync::atomic::Ordering::Relaxed);
            Ok(vec![vec![c as f64, d as f64, (c - d) as f64]])
        }
        20 => {
            let now = std::time::Instant::now();
            let p = ActionPlan::from_triplets(
                need_usize(a, 0, "matrix order")?,
                a[1],
                a[2],
                a[3],
                false,
            )?;
            let elapsed = now.elapsed().as_secs_f64();
            let mut diag = p.diagnostics(need(a, 4)?)?;
            diag.push(elapsed);
            Ok(vec![diag])
        }
        22 => {
            let m = Matrix::new(need_usize(a, 1, "matrix order")?, a[0].to_vec())?;
            let p = ActionPlan::new(&m);
            let n = p.positive_multiply_count(need(a, 2)?, a[3])?;
            Ok(vec![vec![n as f64]])
        }
        23 => {
            use std::sync::atomic::Ordering::Relaxed;
            let created = STRAT_CREATED.load(Relaxed);
            let dropped = STRAT_DROPPED.load(Relaxed);
            Ok(vec![vec![
                created as f64,
                dropped as f64,
                (created - dropped) as f64,
            ]])
        }
        24 => bsm::bridges(a),
        25 => bsm::uniforms(a),
        26 => bsm::whole_tree(a),
        _ => Err(format!("unknown kernel {op}")),
    }
}
/// # Safety
/// `ptrs` and `lens` contain n valid arrays, live until this call returns.
#[no_mangle]
pub unsafe extern "C" fn bgb_kernel(
    op: u32,
    ptrs: *const *const f64,
    lens: *const usize,
    n: usize,
    precision: u32,
) -> *mut Reply {
    let result = catch_unwind(AssertUnwindSafe(|| {
        let ptrs = unsafe { std::slice::from_raw_parts(ptrs, n) };
        let lens = unsafe { std::slice::from_raw_parts(lens, n) };
        let args: Vec<_> = ptrs
            .iter()
            .zip(lens)
            .map(|(&p, &l)| unsafe { std::slice::from_raw_parts(p, l) })
            .collect();
        run(
            op,
            &args,
            if precision != 1 {
                PrecisionMode::LegacyF32
            } else {
                PrecisionMode::StrictF64
            },
            precision == 2,
        )
    }));
    let reply = match result {
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
            error: Some(CString::new("Rust kernel rejected malformed input or panicked").unwrap()),
        },
    };
    Box::into_raw(Box::new(reply))
}
/// # Safety
/// p is a live pointer returned by bgb_kernel.
#[no_mangle]
pub unsafe extern "C" fn bgb_result_error(p: *const Reply) -> *const c_char {
    unsafe { (*p).error.as_ref().map_or(std::ptr::null(), |x| x.as_ptr()) }
}
/// # Safety
/// p is a live pointer returned by bgb_kernel.
#[no_mangle]
pub unsafe extern "C" fn bgb_result_count(p: *const Reply) -> usize {
    unsafe { (*p).buffers.len() }
}
/// # Safety
/// p is live and i is a valid buffer index.
#[no_mangle]
pub unsafe extern "C" fn bgb_result_len(p: *const Reply, i: usize) -> usize {
    unsafe { (&(*p).buffers)[i].len() }
}
/// # Safety
/// p is live and i is a valid buffer index.
#[no_mangle]
pub unsafe extern "C" fn bgb_result_data(p: *const Reply, i: usize) -> *const f64 {
    unsafe { (&(*p).buffers)[i].as_ptr() }
}
/// # Safety
/// p is a uniquely owned live pointer returned by bgb_kernel, freed exactly once.
#[no_mangle]
pub unsafe extern "C" fn bgb_result_free(p: *mut Reply) {
    if !p.is_null() {
        drop(unsafe { Box::from_raw(p) });
    }
}

/// Owned context handle. Creation errors are read and released by the C bridge.
enum Context {
    Single(bgb_step::CompiledModelContext),
    Stratified(bgb_step::CompiledStratifiedContext),
}
pub struct ContextHandle {
    context: Option<Context>,
    error: Option<CString>,
}
unsafe fn buffers<'a>(ptrs: *const *const f64, lens: *const usize, n: usize) -> Vec<&'a [f64]> {
    let p = unsafe { std::slice::from_raw_parts(ptrs, n) };
    let l = unsafe { std::slice::from_raw_parts(lens, n) };
    p.iter()
        .zip(l)
        .map(|(&p, &n)| unsafe { std::slice::from_raw_parts(p, n) })
        .collect()
}
#[no_mangle]
pub unsafe extern "C" fn bgb_context_create(
    ptrs: *const *const f64,
    lens: *const usize,
    n: usize,
    threads: usize,
) -> *mut ContextHandle {
    make_context(|| {
        bgb_step::CompiledModelContext::new(&unsafe { buffers(ptrs, lens, n) }, threads)
            .map(Context::Single)
    })
}
#[no_mangle]
pub unsafe extern "C" fn bgb_strat_context_create(
    ptrs: *const *const f64,
    lens: *const usize,
    n: usize,
    threads: usize,
) -> *mut ContextHandle {
    make_context(|| {
        bgb_step::CompiledStratifiedContext::new(&unsafe { buffers(ptrs, lens, n) }, threads)
            .map(Context::Stratified)
    })
}
fn make_context(create: impl FnOnce() -> Result<Context, String>) -> *mut ContextHandle {
    use std::sync::atomic::Ordering::Relaxed;
    let result = catch_unwind(AssertUnwindSafe(create));
    let handle = match result {
        Ok(Ok(context)) => {
            match &context {
                Context::Single(_) => CONTEXT_CREATED.fetch_add(1, Relaxed),
                Context::Stratified(_) => STRAT_CREATED.fetch_add(1, Relaxed),
            };
            ContextHandle {
                context: Some(context),
                error: None,
            }
        }
        Ok(Err(error)) => ContextHandle {
            context: None,
            error: Some(CString::new(error.replace('\0', " ")).unwrap()),
        },
        Err(_) => ContextHandle {
            context: None,
            error: Some(CString::new("tree context construction panicked").unwrap()),
        },
    };
    Box::into_raw(Box::new(handle))
}
#[no_mangle]
pub unsafe extern "C" fn bgb_context_error(p: *const ContextHandle) -> *const c_char {
    unsafe { (*p).error.as_ref().map_or(std::ptr::null(), |x| x.as_ptr()) }
}
#[no_mangle]
pub unsafe extern "C" fn bgb_context_free(p: *mut ContextHandle) {
    if !p.is_null() {
        let context = unsafe { Box::from_raw(p) };
        use std::sync::atomic::Ordering::Relaxed;
        match &context.context {
            Some(Context::Single(_)) => {
                CONTEXT_DROPPED.fetch_add(1, Relaxed);
            }
            Some(Context::Stratified(_)) => {
                STRAT_DROPPED.fetch_add(1, Relaxed);
            }
            None => {}
        }
        drop(context);
    }
}
#[no_mangle]
pub unsafe extern "C" fn bgb_context_evaluate(
    p: *mut ContextHandle,
    ptrs: *const *const f64,
    lens: *const usize,
    n: usize,
    precision: u32,
) -> *mut Reply {
    let result = catch_unwind(AssertUnwindSafe(|| {
        let context = unsafe { (*p).context.as_mut() }.ok_or("invalid tree context")?;
        let args = unsafe { buffers(ptrs, lens, n) };
        let mode = if precision == 1 {
            PrecisionMode::StrictF64
        } else {
            PrecisionMode::LegacyF32
        };
        match context {
            Context::Single(tree) => tree.evaluate(&args, mode, precision == 2),
            Context::Stratified(tree) => tree.evaluate(&args, mode, precision == 2),
        }
    }));
    let reply = match result {
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
            error: Some(CString::new("tree likelihood panicked").unwrap()),
        },
    };
    Box::into_raw(Box::new(reply))
}

/// Write single-created, single-dropped, stratified-created and stratified-dropped.
#[no_mangle]
pub unsafe extern "C" fn bgb_context_stats(out: *mut usize) {
    use std::sync::atomic::Ordering::Relaxed;
    let counts = [
        CONTEXT_CREATED.load(Relaxed),
        CONTEXT_DROPPED.load(Relaxed),
        STRAT_CREATED.load(Relaxed),
        STRAT_DROPPED.load(Relaxed),
    ];
    unsafe { std::ptr::copy_nonoverlapping(counts.as_ptr(), out, counts.len()) };
}

#[cfg(not(feature = "r-blas"))]
#[no_mangle]
pub extern "C" fn bgb_register_dgemm(function: Option<bgb_numerics::runtime_blas::Dgemm>) {
    bgb_numerics::runtime_blas::register(function);
}

#[cfg(not(feature = "r-blas"))]
#[no_mangle]
pub extern "C" fn bgb_dgemm_registered() -> bool {
    bgb_numerics::runtime_blas::registered()
}

#[cfg(not(feature = "r-blas"))]
#[no_mangle]
pub extern "C" fn bgb_register_dgesv(function: Option<bgb_numerics::runtime_blas::Dgesv>) {
    bgb_numerics::runtime_blas::register_solve(function);
}
