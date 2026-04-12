//! Safe bindings to `dgelsd` from the macOS Accelerate NEWLAPACK interface.
//!
//! Provides access to the LAPACK routine `dgelsd`, which computes the
//! minimum-norm solution to a real linear least squares problem using the
//! singular value decomposition (SVD) with divide-and-conquer.

use std::os::raw::{c_double, c_int};

// Raw FFI binding to `dgelsd$NEWLAPACK` in the Accelerate framework.
// All pointer arguments follow the standard LAPACK Fortran calling convention
// (pass-by-reference for scalars, column-major layout for matrices).
#[link(name = "Accelerate", kind = "framework")]
extern "C" {
    #[link_name = "dgelsd$NEWLAPACK"]
    pub fn dgelsd_newlapack(
        m: *const c_int,
        n: *const c_int,
        nrhs: *const c_int,
        a: *mut c_double,
        lda: *const c_int,
        b: *mut c_double,
        ldb: *const c_int,
        s: *mut c_double,
        rcond: *const c_double,
        rank: *mut c_int,
        work: *mut c_double,
        lwork: *const c_int,
        iwork: *mut c_int,
        info: *mut c_int,
    );
}

/// Computes the minimum-norm least squares solution via SVD (divide-and-conquer).
///
/// Solves the problem  min || b - A x ||_2  where `A` is an `m`-by-`n` matrix.
/// On exit, the first `n` entries of `b` contain the solution vector `x`, and
/// `s` contains the singular values of `A` in descending order.
///
/// # Arguments
///
/// * `m`     - Number of rows of `A`.
/// * `n`     - Number of columns of `A`.
/// * `nrhs`  - Number of right-hand side columns in `b`.
/// * `a`     - Matrix `A` in column-major order; overwritten on exit.
/// * `lda`   - Leading dimension of `a` (>= `m`).
/// * `b`     - Right-hand side matrix; on exit holds the solution. Length >= max(`m`,`n`) * `nrhs`.
/// * `ldb`   - Leading dimension of `b` (>= max(`m`,`n`)).
/// * `s`     - Output singular values of `A`, length >= min(`m`,`n`).
/// * `rcond` - Singular values with `s[i]/s[0] < rcond` are treated as zero.
/// * `rank`  - Output effective rank of `A`.
/// * `work`  - Workspace array. Pass `lwork = -1` to query the optimal size.
/// * `lwork` - Length of `work`. If `-1`, a workspace query is performed.
/// * `iwork` - Integer workspace array.
/// * `info`  - `0` on success; `< 0` if the `i`-th argument is invalid; `> 0` if convergence failed.
#[allow(clippy::too_many_arguments)]
pub fn dgelsd(
    m: i32,
    n: i32,
    nrhs: i32,
    a: &mut [f64],
    lda: i32,
    b: &mut [f64],
    ldb: i32,
    s: &mut [f64],
    rcond: f64,
    rank: &mut i32,
    work: &mut [f64],
    lwork: i32,
    iwork: &mut [i32],
    info: &mut i32,
) {
    unsafe {
        dgelsd_newlapack(
            &m,
            &n,
            &nrhs,
            a.as_mut_ptr(),
            &lda,
            b.as_mut_ptr(),
            &ldb,
            s.as_mut_ptr(),
            &rcond,
            rank,
            work.as_mut_ptr(),
            &lwork,
            iwork.as_mut_ptr(),
            info,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dgelsd_overdetermined_system() {
        // Solve the overdetermined system A x ≈ b where
        //   A = [[1, 2], [3, 4], [5, 6]],  b = [1, 2, 3].
        // Expected least squares solution: x ≈ [0, 0.5].
        let m: i32 = 3;
        let n: i32 = 2;
        let nrhs: i32 = 1;

        // Column-major storage of A.
        let mut a: Vec<f64> = vec![
            1.0, 3.0, 5.0, // column 0
            2.0, 4.0, 6.0, // column 1
        ];

        let ldb = m.max(n) as usize;
        let mut b: Vec<f64> = vec![0.0; ldb];
        b[0] = 1.0;
        b[1] = 2.0;
        b[2] = 3.0;

        let min_mn = (m as usize).min(n as usize);
        let mut s: Vec<f64> = vec![0.0; min_mn];

        let rcond = f64::EPSILON * (m.max(n) as f64);
        let mut rank: i32 = 0;
        let mut info: i32 = 0;

        // Workspace query (lwork = -1).
        let mut work_query: Vec<f64> = vec![0.0; 1];
        let mut iwork_query: Vec<i32> = vec![0; 1];

        dgelsd(
            m, n, nrhs, &mut a, m, &mut b, ldb as i32, &mut s, rcond,
            &mut rank, &mut work_query, -1, &mut iwork_query, &mut info,
        );
        assert_eq!(info, 0, "workspace query failed");

        // Allocate optimal workspace.
        let lwork = work_query[0] as i32;
        let mut work: Vec<f64> = vec![0.0; lwork as usize];

        let smlsiz: i32 = 25;
        let nlvl = if min_mn > 0 {
            ((min_mn as f64 / (smlsiz + 1) as f64).log2().floor() as i32 + 1).max(0)
        } else {
            0
        };
        let liwork = (3 * min_mn as i32 * nlvl + 11 * min_mn as i32).max(1);
        let mut iwork: Vec<i32> = vec![0; liwork as usize];

        // Reset A and b (overwritten by the workspace query).
        a = vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0];
        b = vec![0.0; ldb];
        b[0] = 1.0;
        b[1] = 2.0;
        b[2] = 3.0;

        dgelsd(
            m, n, nrhs, &mut a, m, &mut b, ldb as i32, &mut s, rcond,
            &mut rank, &mut work, lwork, &mut iwork, &mut info,
        );
        assert_eq!(info, 0, "dgelsd execution failed");

        assert!(
            b[0].abs() < 1e-10,
            "x[0] should be ≈ 0, got {}",
            b[0]
        );
        assert!(
            (b[1] - 0.5).abs() < 1e-10,
            "x[1] should be ≈ 0.5, got {}",
            b[1]
        );
    }
}
