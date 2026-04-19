//! C-ABI entry points for the TROP (Triply Robust Panel) estimator.
//!
//! Provides `extern "C"` functions that a Stata plugin can call to execute
//! leave-one-out cross-validation, point estimation, and bootstrap inference
//! for both the Twostep and Joint estimation methods.
//!
//! All matrix arguments follow column-major (Fortran) layout so that Stata
//! matrices can be passed without transposition.

extern crate lapack_src;

pub mod bootstrap;
pub mod distance;
pub mod error;
pub mod estimation;
pub mod loocv;
#[cfg(target_os = "macos")]
pub mod newlapack;
pub mod weights;

pub use error::{TropError, TropResult};

use ndarray::{Array1, Array2, ArrayView2, ShapeBuilder};
use rayon::prelude::*;
use std::slice;

// ---------------------------------------------------------------------------
// Column-major pointer conversion helpers
// ---------------------------------------------------------------------------

/// Constructs a 2-D array view from a raw `f64` pointer in column-major order.
///
/// # Safety
/// `ptr` must be non-null and point to at least `rows * cols` contiguous `f64` values
/// whose memory remains valid for the lifetime of the returned view.
#[inline]
unsafe fn ptr_to_array2<'a>(ptr: *const f64, rows: usize, cols: usize) -> ArrayView2<'a, f64> {
    let slice = slice::from_raw_parts(ptr, rows * cols);
    ArrayView2::from_shape((rows, cols).f(), slice).unwrap()
}

/// Constructs a 2-D array view from a raw `i64` pointer in column-major order.
///
/// # Safety
/// Same requirements as [`ptr_to_array2`].
#[inline]
unsafe fn ptr_to_array2_i64<'a>(ptr: *const i64, rows: usize, cols: usize) -> ArrayView2<'a, i64> {
    let slice = slice::from_raw_parts(ptr, rows * cols);
    ArrayView2::from_shape((rows, cols).f(), slice).unwrap()
}

/// Constructs a 2-D array view from a raw `u8` pointer in column-major order.
///
/// # Safety
/// Same requirements as [`ptr_to_array2`].
#[inline]
unsafe fn ptr_to_array2_u8<'a>(ptr: *const u8, rows: usize, cols: usize) -> ArrayView2<'a, u8> {
    let slice = slice::from_raw_parts(ptr, rows * cols);
    ArrayView2::from_shape((rows, cols).f(), slice).unwrap()
}

/// Writes a 2-D `f64` array to a raw pointer in column-major order.
///
/// # Safety
/// `ptr` must be non-null and point to a buffer of at least `arr.len()` elements.
#[inline]
unsafe fn array2_to_ptr(arr: &Array2<f64>, ptr: *mut f64) {
    let slice = slice::from_raw_parts_mut(ptr, arr.len());
    for ((t, i), val) in arr.indexed_iter() {
        let idx = i * arr.nrows() + t;
        slice[idx] = *val;
    }
}

/// Catches any unwinding panic inside `$body` and converts it to an error code.
macro_rules! catch_panic {
    ($body:expr) => {
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| $body)) {
            Ok(result) => result,
            Err(_) => TropError::RustPanic.code(),
        }
    };
}

// ---------------------------------------------------------------------------
// Twostep method — C ABI exports
// ---------------------------------------------------------------------------

/// Leave-one-out cross-validation grid search for Twostep tuning parameters
/// (λ_time, λ_unit, λ_nn).
///
/// Searches over the Cartesian product of the three grids and writes the
/// best triple, its LOOCV score, and diagnostic counters to the output
/// pointers.
///
/// # Returns
/// `0` on success; a non-zero `TropError` code otherwise.
///
/// # Safety
/// All pointers must be non-null and point to properly sized buffers.
#[no_mangle]
pub unsafe extern "C" fn stata_loocv_grid_search(
    y_ptr: *const f64,
    d_ptr: *const f64,
    control_mask_ptr: *const u8,
    time_dist_ptr: *const i64,
    n_periods: i32,
    n_units: i32,
    lambda_time_grid_ptr: *const f64,
    lambda_time_grid_len: i32,
    lambda_unit_grid_ptr: *const f64,
    lambda_unit_grid_len: i32,
    lambda_nn_grid_ptr: *const f64,
    lambda_nn_grid_len: i32,
    max_iter: i32,
    tol: f64,
    best_lambda_time_out: *mut f64,
    best_lambda_unit_out: *mut f64,
    best_lambda_nn_out: *mut f64,
    best_score_out: *mut f64,
    n_valid_out: *mut i32,
    n_attempted_out: *mut i32,
    first_failed_t_out: *mut i32,
    first_failed_i_out: *mut i32,
) -> i32 {
    catch_panic!({
        if y_ptr.is_null()
            || d_ptr.is_null()
            || control_mask_ptr.is_null()
            || time_dist_ptr.is_null()
            || lambda_time_grid_ptr.is_null()
            || lambda_unit_grid_ptr.is_null()
            || lambda_nn_grid_ptr.is_null()
            || best_lambda_time_out.is_null()
            || best_lambda_unit_out.is_null()
            || best_lambda_nn_out.is_null()
            || best_score_out.is_null()
        {
            return TropError::NullPointer.code();
        }

        let np = n_periods as usize;
        let nu = n_units as usize;

        let y = ptr_to_array2(y_ptr, np, nu);
        let d = ptr_to_array2(d_ptr, np, nu);
        let control_mask = ptr_to_array2_u8(control_mask_ptr, np, nu);
        let time_dist = ptr_to_array2_i64(time_dist_ptr, np, np);

        let lambda_time_grid =
            slice::from_raw_parts(lambda_time_grid_ptr, lambda_time_grid_len as usize);
        let lambda_unit_grid =
            slice::from_raw_parts(lambda_unit_grid_ptr, lambda_unit_grid_len as usize);
        let lambda_nn_grid = slice::from_raw_parts(lambda_nn_grid_ptr, lambda_nn_grid_len as usize);

        let (
            best_time,
            best_unit,
            best_nn,
            best_score,
            n_valid,
            n_attempted,
            first_failed,
        ) = loocv::loocv_grid_search(
            &y,
            &d,
            &control_mask,
            &time_dist,
            lambda_time_grid,
            lambda_unit_grid,
            lambda_nn_grid,
            max_iter as usize,
            tol,
        );

        *best_lambda_time_out = best_time;
        *best_lambda_unit_out = best_unit;
        *best_lambda_nn_out = best_nn;
        *best_score_out = best_score;

        if !n_valid_out.is_null() {
            *n_valid_out = n_valid as i32;
        }
        if !n_attempted_out.is_null() {
            *n_attempted_out = n_attempted as i32;
        }
        if !first_failed_t_out.is_null() {
            *first_failed_t_out = match first_failed {
                Some((t, _)) => t as i32,
                None => -1,
            };
        }
        if !first_failed_i_out.is_null() {
            *first_failed_i_out = match first_failed {
                Some((_, i)) => i as i32,
                None => -1,
            };
        }

        TropError::Success.code()
    })
}

/// Twostep point estimation with fixed tuning parameters.
///
/// For each treated observation (t, i), computes observation-specific weights,
/// fits the additive model Y = α_i + β_t + L_{ti} + τ_{ti} D_{ti} via
/// penalised SVD, and returns the ATT (average of individual τ values) together
/// with the averaged parameter matrices.
///
/// # Returns
/// `0` on success; a non-zero `TropError` code otherwise.
///
/// # Safety
/// All pointers must be non-null and point to properly sized buffers.
#[no_mangle]
pub unsafe extern "C" fn stata_estimate_twostep(
    y_ptr: *const f64,
    d_ptr: *const f64,
    control_mask_ptr: *const u8,
    time_dist_ptr: *const i64,
    n_periods: i32,
    n_units: i32,
    lambda_time: f64,
    lambda_unit: f64,
    lambda_nn: f64,
    max_iter: i32,
    tol: f64,
    att_out: *mut f64,
    tau_ptr: *mut f64,
    alpha_ptr: *mut f64,
    beta_ptr: *mut f64,
    l_ptr: *mut f64,
    n_treated_out: *mut i32,
    n_iterations_out: *mut i32,
    converged_out: *mut i32,
    // Optional per-observation diagnostics.  Both nullable; if non-null each
    // must be pre-allocated to hold N_treated i32 values.  The ordering
    // matches the outer iteration `for t in 0..T { for i in 0..N { if D=1 } }`
    // so Mata can reconstruct (t,i) indices by iterating D in the same way.
    converged_by_obs_ptr: *mut i32,
    n_iters_by_obs_ptr: *mut i32,
) -> i32 {
    catch_panic!({
        if y_ptr.is_null()
            || d_ptr.is_null()
            || control_mask_ptr.is_null()
            || time_dist_ptr.is_null()
            || att_out.is_null()
        {
            return TropError::NullPointer.code();
        }

        let np = n_periods as usize;
        let nu = n_units as usize;

        let y = ptr_to_array2(y_ptr, np, nu);
        let d = ptr_to_array2(d_ptr, np, nu);
        let control_mask = ptr_to_array2_u8(control_mask_ptr, np, nu);
        let time_dist = ptr_to_array2_i64(time_dist_ptr, np, np);

        // Collect (period, unit) indices of treated observations.
        let mut treated_obs: Vec<(usize, usize)> = Vec::new();
        for t in 0..np {
            for i in 0..nu {
                if d[[t, i]] == 1.0 {
                    treated_obs.push((t, i));
                }
            }
        }

        if treated_obs.is_empty() {
            return TropError::NoTreated.code();
        }

        // Map +Inf → 0 (no penalty) and +Inf for λ_nn → large finite cap.
        let lt_eff = if lambda_time.is_infinite() {
            0.0
        } else {
            lambda_time
        };
        let lu_eff = if lambda_unit.is_infinite() {
            0.0
        } else {
            lambda_unit
        };
        let ln_eff = if lambda_nn.is_infinite() {
            1e10
        } else {
            lambda_nn
        };

        // Per-observation estimation (parallelised over treated cells).
        struct ObsResult {
            tau: f64,
            alpha: Array1<f64>,
            beta: Array1<f64>,
            l: Array2<f64>,
            n_iters: usize,
            converged: bool,
        }

        // Share a single UnitDistanceCache across all treated observations.
        // Even with a modest N_treated (~10–30), avoiding the per-call
        // O(T) pairwise distance computation is worthwhile when T is
        // large (e.g. PWT at T=48).
        let dist_cache = distance::UnitDistanceCache::build(&y, &d);

        let obs_results: Vec<Option<ObsResult>> = treated_obs
            .par_iter()
            .map(|(t, i)| {
                let weight_matrix = weights::compute_weight_matrix_cached(
                    &y, &d, &dist_cache, np, nu, *i, *t, lt_eff, lu_eff, &time_dist,
                );

                match estimation::estimate_model(
                    &y,
                    &control_mask,
                    &weight_matrix.view(),
                    ln_eff,
                    np,
                    nu,
                    max_iter as usize,
                    tol,
                    None,
                ) {
                    Some((alpha, beta, l, n_iters, did_converge)) => {
                        let tau = y[[*t, *i]] - alpha[*i] - beta[*t] - l[[*t, *i]];
                        Some(ObsResult {
                            tau,
                            alpha,
                            beta,
                            l,
                            n_iters,
                            converged: did_converge,
                        })
                    }
                    None => None,
                }
            })
            .collect();

        // Aggregate: average α, β, L across successful observations.  The
        // per-obs diagnostics (`converged_by_obs`, `n_iters_by_obs`) are
        // aligned with `treated_obs`: one entry per (t,i) in the same order,
        // with -1 entries marking observations whose solver returned None.
        let mut tau_values: Vec<f64> = Vec::with_capacity(treated_obs.len());
        let mut alpha_sum = Array1::<f64>::zeros(nu);
        let mut beta_sum = Array1::<f64>::zeros(np);
        let mut l_sum = Array2::<f64>::zeros((np, nu));
        let mut n_successful: usize = 0;
        let mut max_iters: usize = 0;
        let mut all_successful_converged = true;
        let mut converged_by_obs: Vec<i32> = Vec::with_capacity(treated_obs.len());
        let mut n_iters_by_obs: Vec<i32> = Vec::with_capacity(treated_obs.len());

        for result in obs_results {
            match result {
                Some(obs) => {
                    tau_values.push(obs.tau);
                    alpha_sum += &obs.alpha;
                    beta_sum += &obs.beta;
                    l_sum += &obs.l;
                    n_successful += 1;
                    if obs.n_iters > max_iters {
                        max_iters = obs.n_iters;
                    }
                    if !obs.converged {
                        all_successful_converged = false;
                    }
                    converged_by_obs.push(if obs.converged { 1 } else { 0 });
                    n_iters_by_obs.push(obs.n_iters as i32);
                }
                None => {
                    // Solver failed (e.g. SVD error); keep slot alignment with
                    // `treated_obs` and mark as -1 for Mata/Stata side.
                    converged_by_obs.push(-1);
                    n_iters_by_obs.push(-1);
                }
            }
        }

        if tau_values.is_empty() {
            return TropError::Convergence.code();
        }

        let att = tau_values.iter().sum::<f64>() / tau_values.len() as f64;

        let n_succ_f64 = n_successful as f64;
        let all_alpha = alpha_sum / n_succ_f64;
        let all_beta = beta_sum / n_succ_f64;
        let all_l = l_sum / n_succ_f64;

        *att_out = att;
        *n_treated_out = tau_values.len() as i32;
        *n_iterations_out = max_iters as i32;
        *converged_out = if all_successful_converged { 1 } else { 0 };

        // Write per-obs diagnostics when the caller requests them.
        if !converged_by_obs_ptr.is_null() {
            let slot = slice::from_raw_parts_mut(converged_by_obs_ptr, converged_by_obs.len());
            slot.copy_from_slice(&converged_by_obs);
        }
        if !n_iters_by_obs_ptr.is_null() {
            let slot = slice::from_raw_parts_mut(n_iters_by_obs_ptr, n_iters_by_obs.len());
            slot.copy_from_slice(&n_iters_by_obs);
        }

        if !tau_ptr.is_null() {
            let tau_slice = slice::from_raw_parts_mut(tau_ptr, tau_values.len());
            tau_slice.copy_from_slice(&tau_values);
        }

        if !alpha_ptr.is_null() {
            let alpha_slice = slice::from_raw_parts_mut(alpha_ptr, nu);
            alpha_slice.copy_from_slice(all_alpha.as_slice().unwrap());
        }

        if !beta_ptr.is_null() {
            let beta_slice = slice::from_raw_parts_mut(beta_ptr, np);
            beta_slice.copy_from_slice(all_beta.as_slice().unwrap());
        }

        if !l_ptr.is_null() {
            array2_to_ptr(&all_l, l_ptr);
        }

        TropError::Success.code()
    })
}

/// Bootstrap variance estimation for the Twostep method.
///
/// Resamples units with replacement `n_bootstrap` times, re-estimates the ATT
/// in each replicate, and returns the standard error, percentile confidence
/// interval, and the full vector of bootstrap estimates.
///
/// # Returns
/// `0` on success; a non-zero `TropError` code otherwise.
///
/// # Safety
/// All pointers must be non-null and point to properly sized buffers.
#[no_mangle]
pub unsafe extern "C" fn stata_bootstrap_trop_variance(
    y_ptr: *const f64,
    d_ptr: *const f64,
    control_mask_ptr: *const u8,
    time_dist_ptr: *const i64,
    n_periods: i32,
    n_units: i32,
    lambda_time: f64,
    lambda_unit: f64,
    lambda_nn: f64,
    n_bootstrap: i32,
    max_iter: i32,
    tol: f64,
    seed: u64,
    alpha: f64,
    estimates_ptr: *mut f64,
    se_out: *mut f64,
    ci_lower_out: *mut f64,
    ci_upper_out: *mut f64,
    n_valid_out: *mut i32,
) -> i32 {
    catch_panic!({
        if y_ptr.is_null()
            || d_ptr.is_null()
            || control_mask_ptr.is_null()
            || time_dist_ptr.is_null()
            || se_out.is_null()
        {
            return TropError::NullPointer.code();
        }

        let np = n_periods as usize;
        let nu = n_units as usize;

        let y = ptr_to_array2(y_ptr, np, nu);
        let d = ptr_to_array2(d_ptr, np, nu);
        let control_mask = ptr_to_array2_u8(control_mask_ptr, np, nu);
        let time_dist = ptr_to_array2_i64(time_dist_ptr, np, np);

        let alpha_eff = if alpha <= 0.0 || alpha >= 1.0 {
            0.05
        } else {
            alpha
        };

        let result = bootstrap::bootstrap_trop_variance_full(
            &y,
            &d,
            &control_mask,
            &time_dist,
            lambda_time,
            lambda_unit,
            lambda_nn,
            n_bootstrap as usize,
            max_iter as usize,
            tol,
            seed,
            alpha_eff,
        );

        *se_out = result.se;

        if !ci_lower_out.is_null() {
            *ci_lower_out = result.ci_lower;
        }
        if !ci_upper_out.is_null() {
            *ci_upper_out = result.ci_upper;
        }
        if !n_valid_out.is_null() {
            *n_valid_out = result.n_valid as i32;
        }

        if !estimates_ptr.is_null() {
            let est_slice = slice::from_raw_parts_mut(estimates_ptr, result.estimates.len());
            est_slice.copy_from_slice(&result.estimates);
        }

        TropError::Success.code()
    })
}

// ---------------------------------------------------------------------------
// Joint method — C ABI exports
// ---------------------------------------------------------------------------

/// Coordinate-descent LOOCV search for Joint method tuning parameters
/// (λ_time, λ_unit, λ_nn).
///
/// Performs univariate sweeps along each parameter axis, then cycles until
/// the selected triple stabilises or `max_cycles` is reached.
///
/// # Returns
/// `0` on success; a non-zero `TropError` code otherwise.
///
/// # Safety
/// All pointers must be non-null and point to properly sized buffers.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn stata_loocv_cycling_search_joint(
    y_ptr: *const f64,
    d_ptr: *const f64,
    control_mask_ptr: *const u8,
    n_periods: i32,
    n_units: i32,
    lambda_time_grid_ptr: *const f64,
    lambda_time_grid_len: i32,
    lambda_unit_grid_ptr: *const f64,
    lambda_unit_grid_len: i32,
    lambda_nn_grid_ptr: *const f64,
    lambda_nn_grid_len: i32,
    max_iter: i32,
    tol: f64,
    max_cycles: i32,
    best_lambda_time_out: *mut f64,
    best_lambda_unit_out: *mut f64,
    best_lambda_nn_out: *mut f64,
    best_score_out: *mut f64,
    n_valid_out: *mut i32,
    n_attempted_out: *mut i32,
    first_failed_t_out: *mut i32,
    first_failed_i_out: *mut i32,
) -> i32 {
    catch_panic!({
        if y_ptr.is_null()
            || d_ptr.is_null()
            || control_mask_ptr.is_null()
            || lambda_time_grid_ptr.is_null()
            || lambda_unit_grid_ptr.is_null()
            || lambda_nn_grid_ptr.is_null()
            || best_lambda_time_out.is_null()
            || best_lambda_unit_out.is_null()
            || best_lambda_nn_out.is_null()
            || best_score_out.is_null()
        {
            return TropError::NullPointer.code();
        }

        let np = n_periods as usize;
        let nu = n_units as usize;

        let y = ptr_to_array2(y_ptr, np, nu);
        let d = ptr_to_array2(d_ptr, np, nu);
        let control_mask = ptr_to_array2_u8(control_mask_ptr, np, nu);

        let lambda_time_grid =
            slice::from_raw_parts(lambda_time_grid_ptr, lambda_time_grid_len as usize);
        let lambda_unit_grid =
            slice::from_raw_parts(lambda_unit_grid_ptr, lambda_unit_grid_len as usize);
        let lambda_nn_grid = slice::from_raw_parts(lambda_nn_grid_ptr, lambda_nn_grid_len as usize);

        let (
            best_time,
            best_unit,
            best_nn,
            best_score,
            n_valid,
            n_attempted,
            first_failed,
        ) = loocv::loocv_cycling_search_joint(
            &y,
            &d,
            &control_mask,
            lambda_time_grid,
            lambda_unit_grid,
            lambda_nn_grid,
            max_iter as usize,
            tol,
            max_cycles as usize,
        );

        *best_lambda_time_out = best_time;
        *best_lambda_unit_out = best_unit;
        *best_lambda_nn_out = best_nn;
        *best_score_out = best_score;

        if !n_valid_out.is_null() {
            *n_valid_out = n_valid as i32;
        }
        if !n_attempted_out.is_null() {
            *n_attempted_out = n_attempted as i32;
        }
        if !first_failed_t_out.is_null() {
            *first_failed_t_out = match first_failed {
                Some((t, _)) => t as i32,
                None => -1,
            };
        }
        if !first_failed_i_out.is_null() {
            *first_failed_i_out = match first_failed {
                Some((_, i)) => i as i32,
                None => -1,
            };
        }

        TropError::Success.code()
    })
}

/// Exhaustive (Cartesian) LOOCV grid search for Joint method tuning parameters
/// (λ_time, λ_unit, λ_nn).
///
/// Evaluates every combination in the Cartesian product of the three grids in
/// parallel and returns the triple that minimises the LOOCV criterion Q(λ).
/// Complexity is O(|Λ_time| · |Λ_unit| · |Λ_nn|).  Prefer this over the
/// coordinate-descent variant when the grid is small enough to afford it, as
/// it matches the Python reference (`diff_diff.trop_global._fit_global`,
/// v3.1.1) exactly.
///
/// The signature mirrors [`stata_loocv_cycling_search_joint`] with the sole
/// exception that `max_cycles` is absent.
///
/// # Returns
/// `0` on success; a non-zero `TropError` code otherwise.
///
/// # Safety
/// All pointers must be non-null and point to properly sized buffers.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn stata_loocv_grid_search_joint(
    y_ptr: *const f64,
    d_ptr: *const f64,
    control_mask_ptr: *const u8,
    n_periods: i32,
    n_units: i32,
    lambda_time_grid_ptr: *const f64,
    lambda_time_grid_len: i32,
    lambda_unit_grid_ptr: *const f64,
    lambda_unit_grid_len: i32,
    lambda_nn_grid_ptr: *const f64,
    lambda_nn_grid_len: i32,
    max_iter: i32,
    tol: f64,
    best_lambda_time_out: *mut f64,
    best_lambda_unit_out: *mut f64,
    best_lambda_nn_out: *mut f64,
    best_score_out: *mut f64,
    n_valid_out: *mut i32,
    n_attempted_out: *mut i32,
    first_failed_t_out: *mut i32,
    first_failed_i_out: *mut i32,
) -> i32 {
    catch_panic!({
        if y_ptr.is_null()
            || d_ptr.is_null()
            || control_mask_ptr.is_null()
            || lambda_time_grid_ptr.is_null()
            || lambda_unit_grid_ptr.is_null()
            || lambda_nn_grid_ptr.is_null()
            || best_lambda_time_out.is_null()
            || best_lambda_unit_out.is_null()
            || best_lambda_nn_out.is_null()
            || best_score_out.is_null()
        {
            return TropError::NullPointer.code();
        }

        let np = n_periods as usize;
        let nu = n_units as usize;

        let y = ptr_to_array2(y_ptr, np, nu);
        let d = ptr_to_array2(d_ptr, np, nu);
        let control_mask = ptr_to_array2_u8(control_mask_ptr, np, nu);

        let lambda_time_grid =
            slice::from_raw_parts(lambda_time_grid_ptr, lambda_time_grid_len as usize);
        let lambda_unit_grid =
            slice::from_raw_parts(lambda_unit_grid_ptr, lambda_unit_grid_len as usize);
        let lambda_nn_grid = slice::from_raw_parts(lambda_nn_grid_ptr, lambda_nn_grid_len as usize);

        let (
            best_time,
            best_unit,
            best_nn,
            best_score,
            n_valid,
            n_attempted,
            first_failed,
        ) = loocv::loocv_grid_search_joint(
            &y,
            &d,
            &control_mask,
            lambda_time_grid,
            lambda_unit_grid,
            lambda_nn_grid,
            max_iter as usize,
            tol,
        );

        *best_lambda_time_out = best_time;
        *best_lambda_unit_out = best_unit;
        *best_lambda_nn_out = best_nn;
        *best_score_out = best_score;

        if !n_valid_out.is_null() {
            *n_valid_out = n_valid as i32;
        }
        if !n_attempted_out.is_null() {
            *n_attempted_out = n_attempted as i32;
        }
        if !first_failed_t_out.is_null() {
            *first_failed_t_out = match first_failed {
                Some((t, _)) => t as i32,
                None => -1,
            };
        }
        if !first_failed_i_out.is_null() {
            *first_failed_i_out = match first_failed {
                Some((_, i)) => i as i32,
                None => -1,
            };
        }

        TropError::Success.code()
    })
}

/// Joint point estimation with fixed tuning parameters.
///
/// Computes global weights δ, then solves the weighted least-squares problem
/// Y = μ + α_i + β_t + L_{ti} + τ D_{ti}.  When λ_nn is effectively
/// infinite the low-rank component L is dropped and a direct WLS solve is
/// used; otherwise coordinate descent alternates between (μ, α, β, τ) and
/// the nuclear-norm–penalised L.
///
/// # Returns
/// `0` on success; a non-zero `TropError` code otherwise.
///
/// # Safety
/// All pointers must be non-null and point to properly sized buffers.
#[no_mangle]
pub unsafe extern "C" fn stata_estimate_joint(
    y_ptr: *const f64,
    d_ptr: *const f64,
    n_periods: i32,
    n_units: i32,
    lambda_time: f64,
    lambda_unit: f64,
    lambda_nn: f64,
    max_iter: i32,
    tol: f64,
    tau_out: *mut f64,
    mu_out: *mut f64,
    alpha_ptr: *mut f64,
    beta_ptr: *mut f64,
    l_ptr: *mut f64,
    n_iterations_out: *mut i32,
    converged_out: *mut i32,
    tau_vec_ptr: *mut f64,
    n_treated_out: *mut i32,
) -> i32 {
    catch_panic!({
        if y_ptr.is_null() || d_ptr.is_null() || tau_out.is_null() || mu_out.is_null() {
            return TropError::NullPointer.code();
        }

        let np = n_periods as usize;
        let nu = n_units as usize;

        let y = ptr_to_array2(y_ptr, np, nu);
        let d = ptr_to_array2(d_ptr, np, nu);

        // Map +Inf → 0 (no penalty) and +Inf for λ_nn → large finite cap.
        let lt_eff = if lambda_time.is_infinite() {
            0.0
        } else {
            lambda_time
        };
        let lu_eff = if lambda_unit.is_infinite() {
            0.0
        } else {
            lambda_unit
        };
        let ln_eff = if lambda_nn.is_infinite() {
            1e10
        } else {
            lambda_nn
        };

        // Locate the earliest period containing any treated cell.
        let mut first_treat_period = np;
        for t in 0..np {
            for i in 0..nu {
                if d[[t, i]] == 1.0 {
                    first_treat_period = first_treat_period.min(t);
                    break;
                }
            }
        }
        let treated_periods = np.saturating_sub(first_treat_period);

        let delta = weights::compute_joint_weights(&y, &d, lt_eff, lu_eff, treated_periods);

        // When λ_nn is large enough, skip the low-rank component entirely.
        // τ is post-hoc: mean residual over treated cells (L ≡ 0 here).
        let result = if ln_eff >= 1e10 {
            estimation::solve_joint_no_lowrank(&y, &delta.view()).map(
                |(mu, alpha, beta)| {
                    let mut tau_sum = 0.0_f64;
                    let mut tau_count = 0usize;
                    for t in 0..np {
                        for i in 0..nu {
                            if d[[t, i]] == 1.0 && y[[t, i]].is_finite() {
                                tau_sum += y[[t, i]] - mu - alpha[i] - beta[t];
                                tau_count += 1;
                            }
                        }
                    }
                    let tau = if tau_count > 0 { tau_sum / tau_count as f64 } else { 0.0 };
                    let l = Array2::<f64>::zeros((np, nu));
                    (mu, alpha, beta, l, tau, 1_usize, true)
                },
            )
        } else {
            estimation::solve_joint_with_lowrank(
                &y,
                &d,
                &delta.view(),
                ln_eff,
                max_iter as usize,
                tol,
            )
        };

        match result {
            Some((mu, alpha, beta, l, tau, n_iters, did_converge)) => {
                *tau_out = tau;
                *mu_out = mu;
                *n_iterations_out = n_iters as i32;
                *converged_out = if did_converge { 1 } else { 0 };

                if !alpha_ptr.is_null() {
                    let alpha_slice = slice::from_raw_parts_mut(alpha_ptr, nu);
                    alpha_slice.copy_from_slice(alpha.as_slice().unwrap());
                }

                if !beta_ptr.is_null() {
                    let beta_slice = slice::from_raw_parts_mut(beta_ptr, np);
                    beta_slice.copy_from_slice(beta.as_slice().unwrap());
                }

                if !l_ptr.is_null() {
                    array2_to_ptr(&l, l_ptr);
                }

                // Paper Eq 13: the joint method assembles τ_it for every treated
                // (i,t) cell as the post-hoc residual Y − μ − α_i − β_t − L_{ti}.
                // The scalar tau_out equals the mean of these values (Eq 1 ATT).
                // Exposing the vector lets users inspect cell-level heterogeneity.
                let mut n_treated_cells: i32 = 0;
                if !tau_vec_ptr.is_null() {
                    let mut tau_values: Vec<f64> = Vec::new();
                    for t in 0..np {
                        for i in 0..nu {
                            if d[[t, i]] == 1.0 && y[[t, i]].is_finite() {
                                tau_values.push(
                                    y[[t, i]] - mu - alpha[i] - beta[t] - l[[t, i]],
                                );
                            }
                        }
                    }
                    n_treated_cells = tau_values.len() as i32;
                    let tau_slice = slice::from_raw_parts_mut(tau_vec_ptr, tau_values.len());
                    tau_slice.copy_from_slice(&tau_values);
                } else {
                    for t in 0..np {
                        for i in 0..nu {
                            if d[[t, i]] == 1.0 && y[[t, i]].is_finite() {
                                n_treated_cells += 1;
                            }
                        }
                    }
                }
                if !n_treated_out.is_null() {
                    *n_treated_out = n_treated_cells;
                }

                TropError::Success.code()
            }
            None => TropError::Convergence.code(),
        }
    })
}

/// Bootstrap variance estimation for the Joint method.
///
/// Resamples units with replacement `n_bootstrap` times, re-estimates τ in
/// each replicate, and returns the standard error, percentile confidence
/// interval, and the full vector of bootstrap estimates.
///
/// # Returns
/// `0` on success; a non-zero `TropError` code otherwise.
///
/// # Safety
/// All pointers must be non-null and point to properly sized buffers.
#[no_mangle]
pub unsafe extern "C" fn stata_bootstrap_trop_variance_joint(
    y_ptr: *const f64,
    d_ptr: *const f64,
    n_periods: i32,
    n_units: i32,
    lambda_time: f64,
    lambda_unit: f64,
    lambda_nn: f64,
    n_bootstrap: i32,
    max_iter: i32,
    tol: f64,
    seed: u64,
    alpha: f64,
    estimates_ptr: *mut f64,
    se_out: *mut f64,
    ci_lower_out: *mut f64,
    ci_upper_out: *mut f64,
    n_valid_out: *mut i32,
) -> i32 {
    catch_panic!({
        if y_ptr.is_null() || d_ptr.is_null() || se_out.is_null() {
            return TropError::NullPointer.code();
        }

        let np = n_periods as usize;
        let nu = n_units as usize;

        let y = ptr_to_array2(y_ptr, np, nu);
        let d = ptr_to_array2(d_ptr, np, nu);

        let alpha_eff = if alpha <= 0.0 || alpha >= 1.0 {
            0.05
        } else {
            alpha
        };

        let result = bootstrap::bootstrap_trop_variance_joint_full(
            &y,
            &d,
            lambda_time,
            lambda_unit,
            lambda_nn,
            n_bootstrap as usize,
            max_iter as usize,
            tol,
            seed,
            alpha_eff,
        );

        *se_out = result.se;

        if !ci_lower_out.is_null() {
            *ci_lower_out = result.ci_lower;
        }
        if !ci_upper_out.is_null() {
            *ci_upper_out = result.ci_upper;
        }
        if !n_valid_out.is_null() {
            *n_valid_out = result.n_valid as i32;
        }

        if !estimates_ptr.is_null() {
            let est_slice = slice::from_raw_parts_mut(estimates_ptr, result.estimates.len());
            est_slice.copy_from_slice(&result.estimates);
        }

        TropError::Success.code()
    })
}

// ---------------------------------------------------------------------------
// Utility — C ABI exports
// ---------------------------------------------------------------------------

/// Computes the N×N unit distance matrix based on pre-treatment outcomes.
///
/// Entry (i, j) measures the Euclidean distance between units i and j over
/// their shared control periods.  The result is written in column-major order.
///
/// # Returns
/// `0` on success; a non-zero `TropError` code otherwise.
///
/// # Safety
/// All pointers must be non-null and point to properly sized buffers.
#[no_mangle]
pub unsafe extern "C" fn stata_compute_unit_distance_matrix(
    y_ptr: *const f64,
    d_ptr: *const f64,
    n_periods: i32,
    n_units: i32,
    dist_ptr: *mut f64,
) -> i32 {
    catch_panic!({
        if y_ptr.is_null() || d_ptr.is_null() || dist_ptr.is_null() {
            return TropError::NullPointer.code();
        }

        let np = n_periods as usize;
        let nu = n_units as usize;

        let y = ptr_to_array2(y_ptr, np, nu);
        let d = ptr_to_array2(d_ptr, np, nu);

        let dist_matrix = distance::compute_unit_distance_matrix_internal(&y, &d);

        let dist_slice = slice::from_raw_parts_mut(dist_ptr, nu * nu);
        for i in 0..nu {
            for j in 0..nu {
                let idx = j * nu + i;
                dist_slice[idx] = dist_matrix[[i, j]];
            }
        }

        TropError::Success.code()
    })
}

/// Returns Twostep weight component vectors for a single treated observation.
///
/// Writes the time weight vector θ (T×1) and the unit weight vector ω (N×1)
/// that would be applied when estimating the model centred on
/// (`target_period`, `target_unit`).
///
/// # Returns
/// `0` on success; a non-zero `TropError` code otherwise.
///
/// # Safety
/// All pointers must be non-null and point to properly sized buffers.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn stata_compute_twostep_weight_vectors(
    y_ptr: *const f64,
    d_ptr: *const f64,
    time_dist_ptr: *const i64,
    n_periods: i32,
    n_units: i32,
    target_unit: i32,
    target_period: i32,
    lambda_time: f64,
    lambda_unit: f64,
    theta_out: *mut f64,
    omega_out: *mut f64,
) -> i32 {
    catch_panic!({
        if y_ptr.is_null()
            || d_ptr.is_null()
            || time_dist_ptr.is_null()
            || theta_out.is_null()
            || omega_out.is_null()
        {
            return TropError::NullPointer.code();
        }

        let np = n_periods as usize;
        let nu = n_units as usize;
        let tu = target_unit as usize;
        let tp = target_period as usize;

        let y = ptr_to_array2(y_ptr, np, nu);
        let d = ptr_to_array2(d_ptr, np, nu);
        let time_dist = ptr_to_array2_i64(time_dist_ptr, np, np);

        // Map +Inf → 0 (no penalty), consistent with estimation functions.
        let lt_eff = if lambda_time.is_infinite() { 0.0 } else { lambda_time };
        let lu_eff = if lambda_unit.is_infinite() { 0.0 } else { lambda_unit };

        let (time_weights, unit_weights) = weights::compute_twostep_weight_vectors(
            &y, &d, np, nu, tu, tp, lt_eff, lu_eff, &time_dist,
        );

        // Write output vectors.
        let theta_slice = slice::from_raw_parts_mut(theta_out, np);
        theta_slice.copy_from_slice(time_weights.as_slice().unwrap());

        let omega_slice = slice::from_raw_parts_mut(omega_out, nu);
        omega_slice.copy_from_slice(unit_weights.as_slice().unwrap());

        TropError::Success.code()
    })
}

/// Returns Joint weight component vectors.
///
/// Writes the global time weight vector δ_time (T×1) and the unit weight
/// vector δ_unit (N×1) used in the Joint estimation objective.
///
/// # Returns
/// `0` on success; a non-zero `TropError` code otherwise.
///
/// # Safety
/// All pointers must be non-null and point to properly sized buffers.
#[no_mangle]
pub unsafe extern "C" fn stata_compute_joint_weight_vectors(
    y_ptr: *const f64,
    d_ptr: *const f64,
    n_periods: i32,
    n_units: i32,
    lambda_time: f64,
    lambda_unit: f64,
    delta_time_out: *mut f64,
    delta_unit_out: *mut f64,
) -> i32 {
    catch_panic!({
        if y_ptr.is_null()
            || d_ptr.is_null()
            || delta_time_out.is_null()
            || delta_unit_out.is_null()
        {
            return TropError::NullPointer.code();
        }

        let np = n_periods as usize;
        let nu = n_units as usize;

        let y = ptr_to_array2(y_ptr, np, nu);
        let d = ptr_to_array2(d_ptr, np, nu);

        // Map +Inf → 0 (no penalty).
        let lt_eff = if lambda_time.is_infinite() { 0.0 } else { lambda_time };
        let lu_eff = if lambda_unit.is_infinite() { 0.0 } else { lambda_unit };

        // Locate the earliest period containing any treated cell.
        let mut first_treat_period = np;
        for t in 0..np {
            for i in 0..nu {
                if d[[t, i]] == 1.0 {
                    first_treat_period = first_treat_period.min(t);
                    break;
                }
            }
        }
        let treated_periods = np.saturating_sub(first_treat_period);

        let (delta_time, delta_unit) = weights::compute_joint_weight_vectors(
            &y, &d, lt_eff, lu_eff, treated_periods,
        );

        // Write output vectors.
        let dt_slice = slice::from_raw_parts_mut(delta_time_out, np);
        dt_slice.copy_from_slice(delta_time.as_slice().unwrap());

        let du_slice = slice::from_raw_parts_mut(delta_unit_out, nu);
        du_slice.copy_from_slice(delta_unit.as_slice().unwrap());

        TropError::Success.code()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_codes() {
        assert_eq!(TropError::Success.code(), 0);
        assert_eq!(TropError::NullPointer.code(), 1);
        assert_eq!(TropError::RustPanic.code(), 8);
    }
}
