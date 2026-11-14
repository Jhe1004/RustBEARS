//! Column-major linear algebra and original DEmat arithmetic. GPL >= 2.
// Indexed loops intentionally preserve the reference reduction order.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]
use bgb_core::PrecisionMode;
mod positive;
mod templates;
use positive::PositivePlan;
pub use templates::CompiledQTemplate;
#[derive(Clone, Debug)]
pub struct Matrix {
    n: usize,
    a: Vec<f64>,
}
impl Matrix {
    pub fn new(n: usize, a: Vec<f64>) -> Result<Self, String> {
        if n == 0
            || n > i32::MAX as usize
            || n.checked_mul(n) != Some(a.len())
            || a.iter().any(|x| !x.is_finite())
        {
            Err("invalid finite square matrix".into())
        } else {
            Ok(Self { n, a })
        }
    }
    #[inline]
    pub fn dim(&self) -> usize {
        self.n
    }
    #[inline]
    pub fn as_slice(&self) -> &[f64] {
        &self.a
    }
    #[inline]
    pub fn into_vec(self) -> Vec<f64> {
        self.a
    }
    pub fn transpose(&self) -> Self {
        let n = self.n;
        let len = n.checked_mul(n).expect("matrix dimension overflow");
        assert_eq!(self.a.len(), len, "matrix buffer length mismatch");
        let mut a = vec![0.0; len];
        for j in 0..n {
            for i in 0..n {
                a[i + n * j] = self.a[j + n * i];
            }
        }
        Self { n, a }
    }
    pub fn eye(n: usize) -> Self {
        assert!(
            n > 0 && n <= i32::MAX as usize,
            "invalid identity matrix dimension"
        );
        let len = n.checked_mul(n).expect("matrix dimension overflow");
        let mut a = vec![0.0; len];
        for i in 0..n {
            a[i + n * i] = 1.0;
        }
        Self { n, a }
    }
    pub fn scaled(&self, s: f64) -> Self {
        Self {
            n: self.n,
            a: self.a.iter().map(|x| x * s).collect(),
        }
    }
    pub fn add(&self, b: &Self, s: f64) -> Self {
        assert_eq!(self.n, b.n, "matrix dimension mismatch in add");
        let len = self.n.checked_mul(self.n).expect("matrix dimension overflow");
        assert_eq!(self.a.len(), len, "left matrix buffer length mismatch");
        assert_eq!(b.a.len(), len, "right matrix buffer length mismatch");
        Self {
            n: self.n,
            a: self.a.iter().zip(&b.a).map(|(x, y)| x + s * y).collect(),
        }
    }
    pub fn norm1(&self) -> f64 {
        self.a
            .chunks(self.n)
            .map(|c| c.iter().map(|x| x.abs()).sum::<f64>())
            .fold(0.0, f64::max)
    }
    pub fn mul(&self, b: &Self) -> Self {
        assert_eq!(self.n, b.n, "matrix dimension mismatch in mul");
        let n = self.n;
        let len = n.checked_mul(n).expect("matrix dimension overflow");
        assert_eq!(self.a.len(), len, "left matrix buffer length mismatch");
        assert_eq!(b.a.len(), len, "right matrix buffer length mismatch");
        let mut c = vec![0.0; len];
        #[cfg(feature = "r-blas")]
        unsafe {
            unsafe extern "C" {
                fn bgb_dgemm(n: i32, a: *const f64, b: *const f64, c: *mut f64);
            }
            bgb_dgemm(n as i32, self.a.as_ptr(), b.a.as_ptr(), c.as_mut_ptr());
        }
        #[cfg(not(feature = "r-blas"))]
        if !runtime_blas::gemm(n, n, 1.0, &self.a, &b.a, &mut c) {
            for j in 0..n {
                for k in 0..n {
                    let v = b.a[k + n * j];
                    for i in 0..n {
                        c[i + n * j] += self.a[i + n * k] * v;
                    }
                }
            }
        }
        Self { n, a: c }
    }
    pub fn mv(&self, v: &[f64]) -> Vec<f64> {
        assert_eq!(v.len(), self.n, "matrix-vector dimension mismatch in mv");
        let len = self.n.checked_mul(self.n).expect("matrix dimension overflow");
        assert_eq!(self.a.len(), len, "matrix buffer length mismatch");
        let mut w = vec![0.0; self.n];
        for (j, &x) in v.iter().enumerate() {
            for (i, y) in w.iter_mut().enumerate() {
                *y += self.a[i + self.n * j] * x;
            }
        }
        w
    }
    pub fn solve(&self, rhs: &Self) -> Result<Self, String> {
        if self.n != rhs.n {
            return Err("matrix dimension mismatch in solve".into());
        }
        let n = self.n;
        let Some(len) = n.checked_mul(n) else {
            return Err("matrix dimension overflow".into());
        };
        if self.a.len() != len || rhs.a.len() != len {
            return Err("matrix buffer length mismatch in solve".into());
        }
        let mut a = self.a.clone();
        let mut b = rhs.a.clone();
        for k in 0..n {
            let mut p = k;
            for i in k + 1..n {
                if a[i + n * k].abs() > a[p + n * k].abs() {
                    p = i;
                }
            }
            if a[p + n * k] == 0.0 || !a[p + n * k].is_finite() {
                return Err("singular Pade denominator".into());
            }
            if p != k {
                for j in 0..n {
                    a.swap(k + n * j, p + n * j);
                    b.swap(k + n * j, p + n * j);
                }
            }
            for i in k + 1..n {
                let f = a[i + n * k] / a[k + n * k];
                a[i + n * k] = 0.0;
                for j in k + 1..n {
                    a[i + n * j] -= f * a[k + n * j];
                }
                for j in 0..n {
                    b[i + n * j] -= f * b[k + n * j];
                }
            }
        }
        for col in 0..n {
            for i in (0..n).rev() {
                let mut v = b[i + n * col];
                for j in i + 1..n {
                    v -= a[i + n * j] * b[j + n * col];
                }
                b[i + n * col] = v / a[i + n * i];
            }
        }
        Ok(Self { n, a: b })
    }
    /// Conservative norm-based scaling and squaring with [13/13] Pade.
    pub fn exp(&self, t: f64) -> Result<Self, String> {
        if !t.is_finite() {
            return Err("non-finite time".into());
        }
        if t == 0.0 {
            return Ok(Self::eye(self.n));
        }
        let x = self.scaled(t);
        let norm = x.norm1();
        if norm == 0.0 {
            return Ok(Self::eye(self.n));
        }
        let s = (norm / 5.371920351148152).log2().ceil().max(0.0) as i32;
        if s > 1023 {
            return Err("exponential scaling overflow".into());
        }
        let a = x.scaled(2.0f64.powi(-s));
        let a2 = a.mul(&a);
        let a4 = a2.mul(&a2);
        let a6 = a4.mul(&a2);
        let id = Self::eye(self.n);
        let b = [
            64764752532480000.0,
            32382376266240000.0,
            7771770303897600.0,
            1187353796428800.0,
            129060195264000.0,
            10559470521600.0,
            670442572800.0,
            33522128640.0,
            1323241920.0,
            40840800.0,
            960960.0,
            16380.0,
            182.0,
            1.0,
        ];
        let u = a.mul(
            &a6.mul(&a6.scaled(b[13]).add(&a4, b[11]).add(&a2, b[9]))
                .add(&a6, b[7])
                .add(&a4, b[5])
                .add(&a2, b[3])
                .add(&id, b[1]),
        );
        let v = a6
            .mul(&a6.scaled(b[12]).add(&a4, b[10]).add(&a2, b[8]))
            .add(&a6, b[6])
            .add(&a4, b[4])
            .add(&a2, b[2])
            .add(&id, b[0]);
        let mut r = v.add(&u, -1.0).solve(&v.add(&u, 1.0))?;
        for _ in 0..s {
            r = r.mul(&r);
        }
        if r.a.iter().any(|x| !x.is_finite()) {
            return Err("matrix exponential overflow".into());
        }
        Ok(r)
    }
}
/// Shared trace shift and sparse operator preanalysis, no probability projection.
/// Owned sparse row storage. Construction validates indices and combines duplicates.
#[derive(Clone, Debug)]
pub struct CsrMatrix {
    pub n: usize,
    pub row_ptr: Vec<usize>,
    pub col_idx: Vec<usize>,
    pub values: Vec<f64>,
}
impl CsrMatrix {
    pub fn from_triplets(
        n: usize,
        rows: &[f64],
        cols: &[f64],
        vals: &[f64],
        transpose: bool,
    ) -> Result<Self, String> {
        if n == 0 || n > i32::MAX as usize || rows.len() != cols.len() || rows.len() != vals.len() {
            return Err("invalid sparse dimensions".into());
        }
        let mut data = std::collections::BTreeMap::<(usize, usize), f64>::new();
        for k in 0..vals.len() {
            let r = rows[k];
            let c = cols[k];
            if !r.is_finite()
                || !c.is_finite()
                || r.fract() != 0.0
                || c.fract() != 0.0
                || r < 1.0
                || c < 1.0
                || r > n as f64
                || c > n as f64
                || !vals[k].is_finite()
            {
                return Err("invalid sparse entry".into());
            }
            let (i, j) = if transpose {
                (c as usize - 1, r as usize - 1)
            } else {
                (r as usize - 1, c as usize - 1)
            };
            *data.entry((i, j)).or_insert(0.0) += vals[k];
        }
        let mut out = Self {
            n,
            row_ptr: vec![0; n + 1],
            col_idx: Vec::with_capacity(data.len()),
            values: Vec::with_capacity(data.len()),
        };
        for ((i, j), x) in data {
            if !x.is_finite() {
                return Err("sparse entry overflow".into());
            }
            if x != 0.0 {
                out.row_ptr[i + 1] += 1;
                out.col_idx.push(j);
                out.values.push(x);
            }
        }
        for i in 0..n {
            out.row_ptr[i + 1] += out.row_ptr[i];
        }
        Ok(out)
    }
    pub fn transpose(&self) -> Self {
        let mut ptr = vec![0; self.n + 1];
        for &j in &self.col_idx {
            ptr[j + 1] += 1;
        }
        for i in 0..self.n {
            ptr[i + 1] += ptr[i];
        }
        let mut next = ptr.clone();
        let mut cols = vec![0; self.values.len()];
        let mut vals = vec![0.0; cols.len()];
        for i in 0..self.n {
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                let j = self.col_idx[k];
                let at = next[j];
                cols[at] = i;
                vals[at] = self.values[k];
                next[j] += 1;
            }
        }
        Self {
            n: self.n,
            row_ptr: ptr,
            col_idx: cols,
            values: vals,
        }
    }
    pub fn dense(&self) -> Matrix {
        let mut a = vec![0.0; self.n * self.n];
        for i in 0..self.n {
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                a[i + self.n * self.col_idx[k]] = self.values[k];
            }
        }
        Matrix { n: self.n, a }
    }
}
/// Shared trace shift and O(9*(n+nnz)) nonnegative-envelope power bounds.
/// No dense powers, sampling estimate or global-norm early stopping.
pub struct ActionPlan {
    n: usize,
    entries: Vec<(usize, usize, f64)>,
    mu: f64,
    norm: f64,
    alphas: Vec<(usize, f64)>,
    positive: Option<PositivePlan>,
}
impl ActionPlan {
    pub fn new(a: &Matrix) -> Self {
        let mut entries = Vec::new();
        for j in 0..a.n {
            for i in 0..a.n {
                let x = a.a[i + a.n * j];
                if x != 0.0 {
                    entries.push((i, j, x));
                }
            }
        }
        Self::from_entries(a.n, entries)
    }
    pub fn from_triplets(
        n: usize,
        rows: &[f64],
        cols: &[f64],
        vals: &[f64],
        transpose: bool,
    ) -> Result<Self, String> {
        Self::from_csr(CsrMatrix::from_triplets(n, rows, cols, vals, transpose)?)
    }
    pub fn from_csr(csr: CsrMatrix) -> Result<Self, String> {
        let n = csr.n;
        if n == 0
            || csr.row_ptr.len() != n + 1
            || csr.row_ptr[0] != 0
            || csr.row_ptr[n] != csr.values.len()
            || csr.col_idx.len() != csr.values.len()
            || csr.row_ptr.windows(2).any(|w| w[0] > w[1])
            || csr.col_idx.iter().any(|&j| j >= n)
            || csr.values.iter().any(|x| !x.is_finite())
        {
            return Err("invalid CSR matrix".into());
        }
        for i in 0..n {
            if csr.col_idx[csr.row_ptr[i]..csr.row_ptr[i + 1]]
                .windows(2)
                .any(|w| w[0] >= w[1])
            {
                return Err("CSR must have unique sorted columns".into());
            }
        }
        let columns = csr.transpose();
        let mut entries = Vec::with_capacity(columns.values.len());
        for j in 0..n {
            for k in columns.row_ptr[j]..columns.row_ptr[j + 1] {
                entries.push((columns.col_idx[k], j, columns.values[k]));
            }
        }
        Ok(Self::from_entries(n, entries))
    }
    fn from_entries(n: usize, raw: Vec<(usize, usize, f64)>) -> Self {
        fn diagonal_transform(
            raw: &[(usize, usize, f64)],
            n: usize,
            divisor: f64,
            diag_add: f64,
        ) -> Vec<(usize, usize, f64)> {
            let mut out = Vec::with_capacity(raw.len() + n);
            let mut k = 0;
            for j in 0..n {
                let mut inserted = false;
                while k < raw.len() && raw[k].1 == j {
                    let (i, _, x) = raw[k];
                    if !inserted && i >= j {
                        if i > j {
                            if diag_add != 0.0 {
                                out.push((j, j, diag_add));
                            }
                        }
                        inserted = true;
                    }
                    let y = x / divisor + if i == j { diag_add } else { 0.0 };
                    if y != 0.0 {
                        out.push((i, j, y));
                    }
                    k += 1;
                }
                if !inserted && diag_add != 0.0 {
                    out.push((j, j, diag_add));
                }
            }
            out
        }
        let mut diag = vec![0.0; n];
        let mut row = vec![0.0; n];
        let mut col = vec![0.0; n];
        let mut colsabs = vec![0.0; n];
        let mut nonnegative = true;
        for &(i, j, x) in &raw {
            row[i] += x;
            col[j] += x;
            colsabs[j] += x.abs();
            if i == j {
                diag[i] = x;
                if x > 0.0 {
                    nonnegative = false;
                }
            } else if x < 0.0 {
                nonnegative = false;
            }
        }
        let rate = diag.iter().map(|x| -x).fold(0.0, f64::max);
        let qnorm = colsabs.iter().copied().fold(0.0, f64::max);
        let residual = row
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
            .min(col.iter().copied().fold(f64::NEG_INFINITY, f64::max));
        let positive = if nonnegative
            && qnorm.is_finite()
            && raw
                .iter()
                .all(|&(i, j, x)| i == j || x == 0.0 || (x / rate).abs() >= f64::MIN_POSITIVE)
            && residual <= 32.0 * f64::EPSILON * qnorm
            && (rate > 0.0 || raw.is_empty())
        {
            let b = diagonal_transform(&raw, n, rate.max(f64::MIN_POSITIVE), 1.0);
            PositivePlan::new(n, b, rate)
        } else {
            None
        };
        let mu = diag.iter().sum::<f64>() / n as f64;
        let entries = diagonal_transform(&raw, n, 1.0, -mu);
        let mut sums = vec![0.0; n];
        for &(_, j, x) in &entries {
            sums[j] += x.abs();
        }
        let norm = sums.into_iter().fold(0.0, f64::max);
        // Recursively bound ||B^k||1 by ||(|B|^T)^k 1||inf. Normalize
        // each stage to avoid power overflow; retain a conservative rounding margin.
        let mut powers = [0.0f64; 10];
        powers[1] = norm;
        let mut v = vec![1.0; n];
        let mut next = vec![0.0; n];
        let mut log_scale = 0.0;
        let mut alphas = Vec::new();
        if norm.is_finite() && norm > 0.0 {
            for k in 1..=9 {
                next.fill(0.0);
                for &(i, j, x) in &entries {
                    next[j] += x.abs() * v[i];
                }
                let scale = next.iter().copied().fold(0.0, f64::max);
                if scale == 0.0 {
                    powers[k] = 0.0;
                    v.fill(0.0);
                } else if scale.is_finite() {
                    log_scale += scale.ln();
                    powers[k] = (log_scale / k as f64).exp()
                        * (1.0 + 64.0 * (n as f64 + 1.0) * k as f64 * f64::EPSILON);
                    for x in &mut next {
                        *x /= scale;
                    }
                    std::mem::swap(&mut v, &mut next);
                } else {
                    powers[k] = norm;
                    v.fill(1.0);
                    log_scale = k as f64 * norm.ln();
                }
            }
            for p in 2..=8 {
                let upper = powers[p].max(powers[p + 1]);
                if upper.is_finite() {
                    alphas.push((p, upper));
                }
            }
        }
        Self {
            n,
            entries,
            mu,
            norm,
            alphas,
            positive,
        }
    }
    pub fn diagnostics(&self, t: f64) -> Result<Vec<f64>, String> {
        let (m, s) = self.parameters(t)?;
        let mut d = vec![
            self.n as f64,
            self.entries.len() as f64,
            self.norm,
            m as f64,
            s as f64,
            self.positive.is_some() as u8 as f64,
        ];
        d.extend(self.alphas.iter().map(|x| x.1));
        Ok(d)
    }
    fn positive_apply(&self, p: &PositivePlan, t: f64, v: &[f64]) -> Result<Vec<f64>, String> {
        p.apply(t, v).map(|x| x.0)
    }
    pub fn positive_multiply_count(&self, t: f64, v: &[f64]) -> Result<usize, String> {
        self.positive
            .as_ref()
            .ok_or("positive plan unavailable")?
            .apply(t, v)
            .map(|x| x.1)
    }
    pub fn parameters(&self, t: f64) -> Result<(usize, usize), String> {
        // Published double-precision theta_m values (Al-Mohy/Higham 2011).
        const THETA: &[(usize, f64)] = &[
            (1, 2.29e-16),
            (2, 2.58e-8),
            (3, 1.39e-5),
            (4, 3.40e-4),
            (5, 2.40e-3),
            (6, 9.07e-3),
            (7, 2.38e-2),
            (8, 5.00e-2),
            (9, 8.96e-2),
            (10, 1.44e-1),
            (11, 2.14e-1),
            (12, 3.00e-1),
            (13, 4.00e-1),
            (14, 5.14e-1),
            (15, 6.41e-1),
            (16, 7.81e-1),
            (17, 9.31e-1),
            (18, 1.09),
            (19, 1.26),
            (20, 1.44),
            (21, 1.62),
            (22, 1.82),
            (23, 2.01),
            (24, 2.22),
            (25, 2.43),
            (26, 2.64),
            (27, 2.86),
            (28, 3.08),
            (29, 3.31),
            (30, 3.54),
            (35, 4.7),
            (40, 6.0),
            (45, 7.2),
            (50, 8.5),
            (55, 9.9),
        ];
        let mut best = (usize::MAX, 1, 1);
        for &(degree, theta) in THETA {
            let scales = (t.abs() * self.norm / theta).ceil().max(1.0) as usize;
            let cost = degree.saturating_mul(scales);
            if cost < best.0 {
                best = (cost, degree, scales)
            }
            for &(p, alpha) in &self.alphas {
                if degree < p * (p - 1) - 1 {
                    continue;
                }
                let scales = (t.abs() * alpha / theta).ceil().max(1.0) as usize;
                let cost = degree.saturating_mul(scales);
                if cost < best.0 {
                    best = (cost, degree, scales)
                }
            }
        }
        if best.2 > 1_000_000 {
            return Err("Taylor scaling limit exceeded".into());
        }
        Ok((best.1, best.2))
    }
    pub fn apply(&self, t: f64, v: &[f64]) -> Result<Vec<f64>, String> {
        if v.len() != self.n || !t.is_finite() || v.iter().any(|x| !x.is_finite()) {
            return Err("invalid action input".into());
        }
        if t == 0.0 {
            return Ok(v.to_vec());
        }
        if t > 0.0 {
            if let Some(p) = &self.positive {
                if p.rate == 0.0 || p.rate * t > 0.0 {
                    return self.positive_apply(p, t, v);
                }
            }
        }
        let (degree, s) = self.parameters(t)?;
        let degree = degree.max(64);
        let h = t / s as f64;
        let eta = (h * self.mu).exp();
        let mut out = v.to_vec();
        for _ in 0..s {
            let mut term = out.clone();
            let mut sum = out.clone();
            let mut next = vec![0.0; self.n];
            for k in 1..=degree {
                next.fill(0.0);
                for &(i, j, a) in &self.entries {
                    next[i] += a * term[j];
                }
                let f = h / k as f64;
                for x in &mut next {
                    *x *= f;
                }
                for (i, x) in next.iter().enumerate() {
                    sum[i] += x;
                }
                std::mem::swap(&mut term, &mut next);
            }
            for x in &mut sum {
                *x *= eta;
            }
            out = sum;
        }
        if out.iter().any(|x| !x.is_finite()) {
            return Err("exponential action overflow".into());
        }
        Ok(out)
    }
}
/// Optional real eigensystem accelerator. Complex spectra are conservatively
/// delegated to Taylor. The residual and both induced condition numbers gate use.
pub struct EigenPlan {
    v: Matrix,
    inverse: Matrix,
    lambda: Vec<f64>,
    pub condition: f64,
    pub residual: f64,
}
impl EigenPlan {
    #[cfg(feature = "r-blas")]
    pub fn new(a: &Matrix) -> Option<Self> {
        let n = a.n;
        let mut wr = vec![0.0; n];
        let mut wi = vec![0.0; n];
        let mut vr = vec![0.0; n * n];
        unsafe {
            unsafe extern "C" {
                fn bgb_eigen(
                    n: i32,
                    a: *const f64,
                    wr: *mut f64,
                    wi: *mut f64,
                    vr: *mut f64,
                ) -> i32;
            }
            if bgb_eigen(
                n as i32,
                a.a.as_ptr(),
                wr.as_mut_ptr(),
                wi.as_mut_ptr(),
                vr.as_mut_ptr(),
            ) != 0
            {
                return None;
            }
        }
        if wi.iter().any(|x| *x != 0.0) {
            return None;
        }
        let v = Matrix::new(n, vr).ok()?;
        let inverse = v.solve(&Matrix::eye(n)).ok()?;
        let norminf = |m: &Matrix| {
            let mut max = 0.0f64;
            for i in 0..n {
                let sum = (0..n).map(|j| m.a[i + n * j].abs()).sum::<f64>();
                max = max.max(sum);
            }
            max
        };
        let condition = (v.norm1() * inverse.norm1()).max(norminf(&v) * norminf(&inverse));
        if !condition.is_finite() || condition > 1e4 {
            return None;
        }
        let mut vd = v.clone();
        for j in 0..n {
            for i in 0..n {
                vd.a[i + n * j] *= wr[j];
            }
        }
        let residual = vd.mul(&inverse).add(a, -1.0).norm1() / a.norm1().max(f64::MIN_POSITIVE);
        if !residual.is_finite() || residual > 1e-13 {
            return None;
        }
        Some(Self {
            v,
            inverse,
            lambda: wr,
            condition,
            residual,
        })
    }
    #[cfg(not(feature = "r-blas"))]
    pub fn new(_a: &Matrix) -> Option<Self> {
        None
    }
    pub fn apply(&self, t: f64, x: &[f64]) -> Vec<f64> {
        let mut y = self.inverse.mv(x);
        for (j, v) in y.iter_mut().enumerate() {
            *v *= (t * self.lambda[j]).exp();
        }
        self.v.mv(&y)
    }
    pub fn apply_batch(&self, times: &[f64], vectors: &[f64]) -> Vec<f64> {
        let n = self.v.n;
        let cols = times.len();
        let mut y = vec![0.0; n * cols];
        let mut out = vec![0.0; n * cols];
        #[cfg(feature = "r-blas")]
        unsafe {
            unsafe extern "C" {
                fn bgb_dgemm_rect(n: i32, cols: i32, a: *const f64, b: *const f64, c: *mut f64);
            }
            bgb_dgemm_rect(
                n as i32,
                cols as i32,
                self.inverse.a.as_ptr(),
                vectors.as_ptr(),
                y.as_mut_ptr(),
            );
        }
        #[cfg(not(feature = "r-blas"))]
        for j in 0..cols {
            y[j * n..(j + 1) * n].copy_from_slice(&self.inverse.mv(&vectors[j * n..(j + 1) * n]));
        }
        for j in 0..cols {
            for i in 0..n {
                y[i + n * j] *= (times[j] * self.lambda[i]).exp();
            }
        }
        #[cfg(feature = "r-blas")]
        unsafe {
            unsafe extern "C" {
                fn bgb_dgemm_rect(n: i32, cols: i32, a: *const f64, b: *const f64, c: *mut f64);
            }
            bgb_dgemm_rect(
                n as i32,
                cols as i32,
                self.v.a.as_ptr(),
                y.as_ptr(),
                out.as_mut_ptr(),
            );
        }
        #[cfg(not(feature = "r-blas"))]
        for j in 0..cols {
            out[j * n..(j + 1) * n].copy_from_slice(&self.v.mv(&y[j * n..(j + 1) * n]));
        }
        out
    }
}
pub fn demat(
    states: &[Vec<i32>],
    na: usize,
    d: &[f64],
    e: &[f64],
    amat: &[f64],
    normalize: bool,
    coo: bool,
    min: f64,
    mode: PrecisionMode,
) -> Result<Vec<Vec<f64>>, String> {
    if na == 0 || d.len() != na * na || amat.len() != na * na || e.len() < na {
        return Err("invalid DEmat dimensions".into());
    }
    let n = states.len();
    for s in states {
        if s.is_empty() || s.iter().any(|&i| i >= na as i32) {
            return Err("invalid state index".into());
        }
    }
    let mut q = vec![0.0; n * n];
    let mut ia = Vec::new();
    let mut ja = Vec::new();
    let mut vals = Vec::new();
    let threshold = mode.round(if coo { min } else { 1e-15 });
    for (i, a) in states.iter().enumerate() {
        for (j, b) in states.iter().enumerate() {
            let mut x = 0.0;
            let mut legal = false;
            if a.len() == 1 && b.len() == 1 && a[0] >= 0 && b[0] >= 0 && a[0] != b[0] {
                legal = true;
                x = mode.round(amat[a[0] as usize + na * b[0] as usize]);
            } else if a.len() + 1 == b.len() && a.iter().all(|x| b.contains(x)) {
                legal = true;
                let extra = *b.iter().find(|x| !a.contains(x)).ok_or("missing area")?;
                for &from in a {
                    if from < 0 {
                        return Err("invalid empty-range expansion".into());
                    }
                    x = mode.round(x + mode.round(d[from as usize + na * extra as usize]));
                }
            } else if b[0] < 0 && a[0] >= 0 && a.len() == 1 {
                legal = true;
                x = mode.round(e[a[0] as usize]);
            } else if a.len() == b.len() + 1 && b.iter().all(|x| a.contains(x)) {
                legal = true;
                let lost = *a.iter().find(|x| !b.contains(x)).ok_or("missing area")?;
                x = mode.round(e[lost as usize]);
            }
            q[i + n * j] = x;
            if coo && legal && x >= threshold && i != j {
                ia.push((j + 1) as f64);
                ja.push((i + 1) as f64);
                vals.push(x);
            }
        }
    }
    if coo {
        let mut legacy_count = vals.len().saturating_sub(1);
        if normalize {
            for i in 0..n {
                let mut sum = 0.0;
                let end = if mode == PrecisionMode::LegacyF32 {
                    legacy_count
                } else {
                    vals.len()
                };
                for k in 0..end {
                    if ja[k] == (i + 1) as f64 && ia[k] != (i + 1) as f64 {
                        sum = mode.round(sum + vals[k]);
                    }
                }
                if sum >= threshold {
                    ia.push((i + 1) as f64);
                    ja.push((i + 1) as f64);
                    vals.push(-sum);
                    legacy_count += 1;
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
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sparse_and_dense_action_keep_long_path_probability() {
        let n = 30;
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for i in 0..n - 1 {
            rows.extend([(i + 1) as f64, (i + 1) as f64]);
            cols.extend([(i + 1) as f64, (i + 2) as f64]);
            vals.extend([-1.0, 1.0]);
        }
        let csr = CsrMatrix::from_triplets(n, &rows, &cols, &vals, false).unwrap();
        let q = csr.dense();
        let truth = (-0.5f64).exp() * 0.5f64.powi(20) / (1..=20).map(|k| k as f64).product::<f64>();
        let mut v = vec![0.0; n];
        v[20] = 1.0;
        for plan in [
            ActionPlan::new(&q),
            ActionPlan::from_csr(csr).unwrap(),
            ActionPlan::from_triplets(n, &rows, &cols, &vals, false).unwrap(),
        ] {
            let x = plan.apply(0.5, &v).unwrap()[0];
            assert!(x > 0.0);
            assert!((x / truth - 1.0).abs() < 1e-13);
            assert_eq!(plan.alphas.len(), 7);
        }
    }
    #[test]
    fn identities_and_jordan() {
        let q = Matrix::new(2, vec![-1.0, 0.0, 1.0, -1.0]).unwrap();
        let p = q.exp(1.0).unwrap();
        let z = (-1.0f64).exp();
        for (i, v) in [z, 0.0, z, z].iter().enumerate() {
            assert!((p.a[i] - v).abs() < 1e-15)
        }
        let v = vec![2.0, 3.0];
        assert_eq!(ActionPlan::new(&q).apply(0.0, &v).unwrap(), v);
        let w = ActionPlan::new(&q).apply(1.0, &[0.0, 1.0]).unwrap();
        assert!((w[0] - z).abs() < 1e-15);
    }
    #[test]
    fn matrix_dimension_guards_reject_mismatched_orders() {
        let m2 = Matrix::new(2, vec![1.0, 0.0, 0.0, 1.0]).unwrap();
        let m1 = Matrix::new(1, vec![2.0]).unwrap();
        assert!(m2.solve(&m1).is_err());
        assert!(std::panic::catch_unwind(|| m2.mul(&m1)).is_err());
        assert!(std::panic::catch_unwind(|| m2.add(&m1, 1.0)).is_err());
        assert!(std::panic::catch_unwind(|| m2.mv(&[1.0])).is_err());
        let mut out = vec![0.0; 4];
        assert!(!runtime_blas::gemm(2, 2, 1.0, &[1.0; 4], &[1.0; 1], &mut out));
        assert!(runtime_blas::solve(2, &[1.0; 4], &mut [1.0; 1]).is_err());
    }
}

impl Matrix {
    // EXPOKIT's DGEXX accumulation order, ported from the reference BLAS loop.
    // Copyright notices for the original algorithm are in inst/NOTICE-Rust-port.md.
    fn legacy_mul(&self, b: &Self, alpha: f64, fused: bool) -> Self {
        assert_eq!(self.n, b.n, "matrix dimension mismatch in legacy_mul");
        let n = self.n;
        let len = n.checked_mul(n).expect("matrix dimension overflow");
        assert_eq!(self.a.len(), len, "left matrix buffer length mismatch");
        assert_eq!(b.a.len(), len, "right matrix buffer length mismatch");
        let mut out = vec![0.0; len];
        #[cfg(not(feature = "r-blas"))]
        if n >= 64 && runtime_blas::gemm(n, n, alpha, &self.a, &b.a, &mut out) {
            return Self { n, a: out };
        }
        for j in 0..n {
            for k in 0..n {
                if b.a[k + n * j] != 0.0 {
                    let temp = alpha * b.a[k + n * j];
                    for i in 0..n {
                        out[i + n * j] = if fused {
                            temp.mul_add(self.a[i + n * k], out[i + n * j])
                        } else {
                            out[i + n * j] + temp * self.a[i + n * k]
                        };
                    }
                }
            }
        }
        Self { n, a: out }
    }
    /// Historical EXPOKIT DGPADM(ideg=6). In compatibility mode even tiny
    /// rounding changes can reorder equal probabilities in R's weighted sampler.
    pub fn legacy_exp(&self, t: f64, fused: bool) -> Result<Self, String> {
        if !t.is_finite() {
            return Err("non-finite time".into());
        }
        let n = self.n;
        let mut rows = vec![0.0; n];
        for j in 0..n {
            for i in 0..n {
                rows[i] += self.a[i + n * j].abs();
            }
        }
        let norm = (t * rows.into_iter().fold(0.0, f64::max)).abs();
        if norm == 0.0 {
            return Ok(Self::eye(n));
        }
        let ns = ((norm.ln() / 2.0f64.ln()) as i32 + 2).max(0);
        if ns > 1023 {
            return Err("legacy scaling overflow".into());
        }
        let scale = t / 2.0f64.powi(ns);
        let scale2 = scale * scale;
        let mut coefficients = [0.0; 7];
        coefficients[0] = 1.0;
        for k in 1..=6 {
            coefficients[k] = (coefficients[k - 1] * (7 - k) as f64) / (k * (13 - k)) as f64;
        }
        let h2 = self.legacy_mul(self, scale2, fused);
        let mut p = Self::eye(n).scaled(coefficients[5]);
        let mut q = Self::eye(n).scaled(coefficients[6]);
        let mut odd = true;
        for k in (1..=5).rev() {
            let mut tmp = if odd {
                q.legacy_mul(&h2, 1.0, fused)
            } else {
                p.legacy_mul(&h2, 1.0, fused)
            };
            for j in 0..n {
                tmp.a[j + n * j] += coefficients[k - 1];
            }
            if odd {
                q = tmp
            } else {
                p = tmp
            }
            odd = !odd;
        }
        if odd {
            q = q.legacy_mul(self, scale, fused);
        } else {
            p = p.legacy_mul(self, scale, fused);
        }
        for i in 0..q.a.len() {
            q.a[i] += -p.a[i];
        }
        #[cfg(feature = "r-blas")]
        {
            unsafe extern "C" {
                fn bgb_dgesv(n: i32, a: *mut f64, b: *mut f64) -> i32;
            }
            let info = unsafe { bgb_dgesv(n as i32, q.a.as_mut_ptr(), p.a.as_mut_ptr()) };
            if info != 0 {
                return Err("legacy Pade solve failed".into());
            }
        }
        #[cfg(not(feature = "r-blas"))]
        {
            if n < 64 || !runtime_blas::solve(n, &q.a, &mut p.a)? {
                p = q.solve(&p)?;
            }
        }
        for x in &mut p.a {
            *x *= 2.0;
        }
        for j in 0..n {
            p.a[j + n * j] += 1.0;
        }
        if ns == 0 && odd {
            for x in &mut p.a {
                *x = -*x;
            }
        }
        for _ in 0..ns {
            p = p.legacy_mul(&p, 1.0, fused);
        }
        Ok(p)
    }
}

impl Matrix {
    /// Positive uniformization for conservative or killed generators (either orientation).
    /// This is a series for exp(tQ), not a projection or renormalization.
    /// Other finite matrices retain the general Pade path.
    pub fn strict_exp(&self, t: f64) -> Result<Self, String> {
        if !t.is_finite() {
            return Err("non-finite time".into());
        }
        if t == 0.0 {
            return Ok(Self::eye(self.n));
        }
        let n = self.n;
        let mut rows = vec![0.0; n];
        let mut cols = vec![0.0; n];
        let mut rate = 0.0_f64;
        for j in 0..n {
            for i in 0..n {
                let x = self.a[i + n * j];
                if i != j && x < 0.0 {
                    return self.exp(t);
                }
                rows[i] += x;
                cols[j] += x;
                if i == j {
                    rate = rate.max(-x);
                    if x > 0.0 {
                        return self.exp(t);
                    }
                }
            }
        }
        let residual = rows
            .iter()
            .fold(0.0_f64, |m, x| m.max(*x))
            .min(cols.iter().fold(0.0_f64, |m, x| m.max(*x)));
        if t < 0.0 || residual > 32.0 * f64::EPSILON * self.norm1() {
            return self.exp(t);
        }
        if rate == 0.0 {
            return Ok(Self::eye(n));
        }
        let z = rate * t;
        if !z.is_finite() {
            return Err("uniformization scaling overflow".into());
        }
        let power = (z / 0.5).log2().ceil().max(0.0) as i32;
        if power > 1023 {
            return Err("uniformization scaling overflow".into());
        }
        let h = t / 2.0_f64.powi(power);
        let z = rate * h;
        let plan = ActionPlan::new(self);
        if let Some(p) = &plan.positive {
            let mut a = vec![0.0; n * n];
            let mut basis = vec![0.0; n];
            let mut supported = true;
            for j in 0..n {
                basis.fill(0.0);
                basis[j] = 1.0;
                match p.apply(h, &basis) {
                    Ok((col, _)) => a[j * n..(j + 1) * n].copy_from_slice(&col),
                    Err(_) => {
                        supported = false;
                        break;
                    }
                }
            }
            if supported {
                let mut out = Self { n, a };
                for _ in 0..power {
                    out = out.mul(&out);
                }
                return Ok(out);
            }
        }
        let b = Self::eye(n).add(
            &Self {
                n,
                a: self.a.iter().map(|x| x / rate).collect(),
            },
            1.0,
        );
        let mut term = Self::eye(n);
        let mut sum = term.clone();
        // Scaled Poisson mean <= 0.5; 64-term tail is below 1e-110.
        for k in 1..=64 {
            term = term.mul(&b).scaled(z / k as f64);
            sum = sum.add(&term, 1.0);
        }
        sum = sum.scaled((-z).exp());
        for _ in 0..power {
            sum = sum.mul(&sum);
        }
        Ok(sum)
    }
}

/// Historical EXPOKIT Arnoldi action for the rexpokit compatibility API.
/// Its public R shell fixes m=n-1 and tol=0.01; these approximations are
/// observable behavior, distinct from the accurate Taylor action interface.
#[allow(clippy::approx_constant)] // Preserve the historical 3.14 stepsize heuristic.
pub fn legacy_arnoldi(
    n: usize,
    m: usize,
    t: f64,
    input: &[f64],
    tol: f64,
    anorm: f64,
    rows: &[f64],
    cols: &[f64],
    vals: &[f64],
    markov: bool,
    fused: bool,
) -> Result<Vec<f64>, String> {
    if input.len() != n
        || m == 0
        || m >= n
        || rows.len() != cols.len()
        || rows.len() != vals.len()
        || !t.is_finite()
        || !anorm.is_finite()
        || anorm < 0.0
    {
        return Err("invalid legacy Arnoldi input".into());
    }
    let madd = |a: f64, b: f64, c: f64| if fused { a.mul_add(b, c) } else { c + a * b };
    let norm = |x: &[f64]| {
        let mut scale = 0.0_f64;
        let mut ssq = 1.0_f64;
        for &v in x {
            if v != 0.0 {
                let a = v.abs();
                if scale < a {
                    ssq = 1.0 + ssq * (scale / a).powi(2);
                    scale = a;
                } else {
                    ssq += (a / scale).powi(2);
                }
            }
        }
        scale * ssq.sqrt()
    };
    let mv = |v: &[f64]| -> Result<Vec<f64>, String> {
        let mut out = vec![0.0; n];
        for k in 0..vals.len() {
            let i = rows[k] as usize;
            let j = cols[k] as usize;
            if i == 0 || i > n || j == 0 || j > n {
                return Err("legacy COO index out of range".into());
            }
            out[i - 1] = madd(vals[k], v[j - 1], out[i - 1]);
        }
        Ok(out)
    };
    let rounded = |x: f64| {
        let exponent = (x.log10() - 0.1_f64.sqrt()).round() as i32 - 1;
        let unit = 10.0_f64.powi(exponent);
        (x / unit + 0.55).trunc() * unit
    };
    let mut column_sums = vec![0.0; n];
    let mut nonnegative_offdiag = true;
    for k in 0..vals.len() {
        let i = rows[k] as usize;
        let j = cols[k] as usize;
        if i == 0 || i > n || j == 0 || j > n {
            return Err("legacy COO index out of range".into());
        }
        column_sums[j - 1] += vals[k];
        if i != j && vals[k] < 0.0 {
            nonnegative_offdiag = false;
        }
    }
    let input_mass = input.iter().sum::<f64>();
    let probability_path = markov
        && t >= 0.0
        && input.iter().all(|x| *x >= 0.0)
        && (input_mass - 1.0).abs() <= 32.0 * f64::EPSILON
        && nonnegative_offdiag
        && column_sums
            .iter()
            .all(|x| x.abs() <= 32.0 * f64::EPSILON * anorm);
    let mut w = input.to_vec();
    let mut beta = norm(&w);
    if t == 0.0 || beta == 0.0 || anorm == 0.0 {
        return Ok(w);
    }
    let tol = if tol <= f64::EPSILON {
        f64::EPSILON.sqrt()
    } else {
        tol
    };
    let mut xm = 1.0 / m as f64;
    let p =
        tol * ((m + 1) as f64 / 2.72).powi((m + 1) as i32) * (2.0 * 3.14 * (m + 1) as f64).sqrt();
    let mut next = rounded((1.0 / anorm) * (p / (4.0 * beta * anorm)).powf(xm));
    let end = t.abs();
    let mut now = 0.0;
    let size = m + 2;
    for _step in 0..if markov { 500 } else { 1000 } {
        if now >= end {
            return Ok(w);
        }
        let mut dt = (end - now).min(next);
        if !dt.is_finite() || dt <= 0.0 {
            return Err("legacy Arnoldi step underflow".into());
        }
        let mut basis = vec![vec![0.0; n]; m + 1];
        for i in 0..n {
            basis[0][i] = (1.0 / beta) * w[i];
        }
        let mut h = Matrix::new(size, vec![0.0; size * size])?;
        let mut correction: usize = 2;
        let mut breakdown = m;
        let mut avnorm = 0.0;
        for j in 0..m {
            let mut v = mv(&basis[j])?;
            for i in 0..=j {
                let mut dot = 0.0;
                for k in 0..n {
                    dot = madd(basis[i][k], v[k], dot);
                }
                h.a[i + size * j] = dot;
                for k in 0..n {
                    v[k] = madd(-dot, basis[i][k], v[k]);
                }
            }
            let length = norm(&v);
            if length <= 1e-7 {
                correction = 0;
                breakdown = j + 1;
                dt = end - now;
                break;
            }
            h.a[j + 1 + size * j] = length;
            for k in 0..n {
                basis[j + 1][k] = (1.0 / length) * v[k];
            }
        }
        if correction != 0 {
            avnorm = norm(&mv(&basis[m])?);
        }
        h.a[m + 1 + size * m] = 1.0;
        let mut accepted = None;
        for _reject in 0..10_000 {
            let nx = breakdown + correction;
            let mut data = vec![0.0; nx * nx];
            for j in 0..nx {
                for i in 0..nx {
                    data[i + nx * j] = h.a[i + size * j];
                }
            }
            let e = Matrix::new(nx, data)?.legacy_exp(t.signum() * dt, fused)?;
            let error = if correction == 0 {
                tol
            } else {
                let p1 = e.a[m].abs() * beta;
                let p2 = e.a[m + 1].abs() * beta * avnorm;
                if p1 > 10.0 * p2 {
                    xm = 1.0 / m as f64;
                    p2
                } else if p1 > p2 {
                    xm = 1.0 / m as f64;
                    p1 * p2 / (p1 - p2)
                } else {
                    xm = 1.0 / (m as f64 - 1.0);
                    p1
                }
            };
            if correction != 0 && error > 1.2 * dt * tol {
                dt = rounded(0.9 * dt * (dt * tol / error).powf(xm));
                if dt <= 0.0 || !dt.is_finite() {
                    return Err("legacy Arnoldi rejected step underflow".into());
                }
                continue;
            }
            accepted = Some((e, error));
            break;
        }
        let (e, error) = accepted.ok_or("legacy Arnoldi rejection limit exceeded")?;
        w.fill(0.0);
        for j in 0..breakdown + correction.saturating_sub(1) {
            let temp = beta * e.a[j];
            for i in 0..n {
                w[i] = madd(temp, basis[j][i], w[i]);
            }
        }
        beta = norm(&w);
        // DMEXPV is specifically the old Markov-probability interface. No
        // clipping is applied in the general DGEXPV/Taylor/StrictF64 paths.
        if probability_path {
            for x in &mut w {
                if *x < 0.0 {
                    *x = 0.0;
                }
            }
        }
        now += dt;
        if now >= end {
            return Ok(w);
        }
        next = rounded(0.9 * dt * (dt * tol / error).powf(xm));
    }
    Ok(w)
}

#[cfg(not(feature = "r-blas"))]
pub mod runtime_blas {
    use std::ffi::c_char;
    use std::sync::atomic::{AtomicPtr, Ordering};

    pub type Dgemm = unsafe extern "C" fn(
        *const c_char,
        *const c_char,
        *const i32,
        *const i32,
        *const i32,
        *const f64,
        *const f64,
        *const i32,
        *const f64,
        *const i32,
        *const f64,
        *mut f64,
        *const i32,
    );
    pub type Dgesv = unsafe extern "C" fn(
        *const i32,
        *const i32,
        *mut f64,
        *const i32,
        *mut i32,
        *mut f64,
        *const i32,
        *mut i32,
    );
    static DGEMM: AtomicPtr<()> = AtomicPtr::new(std::ptr::null_mut());
    static DGESV: AtomicPtr<()> = AtomicPtr::new(std::ptr::null_mut());

    pub fn register_solve(function: Option<Dgesv>) {
        let address = function.map_or(std::ptr::null_mut(), |f| f as *mut ());
        DGESV.store(address, Ordering::Release);
    }

    pub(crate) fn solve(n: usize, a: &[f64], b: &mut [f64]) -> Result<bool, String> {
        let Some(nn) = n.checked_mul(n) else {
            return Err("matrix dimension overflow".into());
        };
        if a.len() != nn || b.len() != nn {
            return Err("matrix buffer length mismatch in solve".into());
        }
        let f = DGESV.load(Ordering::Acquire);
        if f.is_null() || n == 0 || n > i32::MAX as usize {
            return Ok(false);
        }
        let f: Dgesv = unsafe { std::mem::transmute(f) };
        let mut a = a.to_vec();
        let mut pivots = vec![0; n];
        let n = n as i32;
        let mut status = 0;
        unsafe {
            f(
                &n,
                &n,
                a.as_mut_ptr(),
                &n,
                pivots.as_mut_ptr(),
                b.as_mut_ptr(),
                &n,
                &mut status,
            );
        }
        if status != 0 {
            return Err("native LAPACK Pade solve failed".into());
        }
        Ok(true)
    }

    /// The provider must keep this native symbol alive for the process.
    pub fn register(function: Option<Dgemm>) {
        let address = function.map_or(std::ptr::null_mut(), |f| f as *mut ());
        DGEMM.store(address, Ordering::Release);
    }

    pub fn registered() -> bool {
        !DGEMM.load(Ordering::Acquire).is_null()
    }

    pub(crate) fn gemm(
        n: usize,
        cols: usize,
        alpha: f64,
        a: &[f64],
        b: &[f64],
        c: &mut [f64],
    ) -> bool {
        let Some(nn) = n.checked_mul(n) else {
            return false;
        };
        let Some(nc) = n.checked_mul(cols) else {
            return false;
        };
        if a.len() != nn || b.len() != nc || c.len() != nc {
            return false;
        }
        let f = DGEMM.load(Ordering::Acquire);
        if f.is_null() || n == 0 || cols == 0 || n > i32::MAX as usize || cols > i32::MAX as usize {
            return false;
        }
        let f: Dgemm = unsafe { std::mem::transmute(f) };
        let (n, cols) = (n as i32, cols as i32);
        let t = b'N' as c_char;
        unsafe {
            f(
                &t,
                &t,
                &n,
                &cols,
                &n,
                &alpha,
                a.as_ptr(),
                &n,
                b.as_ptr(),
                &n,
                &0.0,
                c.as_mut_ptr(),
                &n,
            );
        }
        true
    }
}
