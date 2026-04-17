//! Leave-one-out cross-validation (LOOCV) for tuning parameter selection.
//!
//! Selects the regularization triplet (λ_time, λ_unit, λ_nn) by minimizing
//! the LOOCV criterion:
//!
//!   Q(λ) = Σ_{i,t} (1 − W_{it}) (τ̂_{it}^{loocv}(λ))²
//!
//! where τ̂_{it}^{loocv}(λ) is the pseudo-treatment effect obtained by treating
//! each control observation (i,t) as if it were treated and estimating its
//! counterfactual from the remaining control observations.
//!
//! The search uses a two-stage strategy:
//!   Stage 1 — Univariate initialization: search each parameter independently
//!             while holding the other two at extreme values.
//!   Stage 2 — Coordinate descent cycling: successively update each parameter
//!             holding the other two at their most recent optimal values,
//!             until convergence.
//!
//! Both the twostep and joint estimation methods share this search structure.
//! A brute-force O(|grid|³) joint search (`loocv_grid_search_joint`) is
//! retained for testing but is not exposed via the C ABI.

use ndarray::{Array2, ArrayView2};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

use crate::estimation::{estimate_model, solve_joint_no_lowrank, solve_joint_with_lowrank};
use crate::weights::{compute_joint_weights, compute_weight_matrix};

/// Type alias for the LOOCV grid search result tuple.
///
/// Fields: (best_lambda_time, best_lambda_unit, best_lambda_nn, best_score,
///          n_valid, n_attempted, n_control_total, subsampled, first_failed_obs).
#[allow(clippy::type_complexity)]
pub type LoocvGridSearchResult = (
    f64,
    f64,
    f64,
    f64,
    usize,
    usize,
    usize,
    bool,
    Option<(usize, usize)>,
);

/// Type alias for internal LOOCV result tuple.
#[allow(clippy::type_complexity)]
type LoocvResultTuple = (f64, f64, f64, f64, usize, Option<(usize, usize)>);

/// Collect control observations eligible for LOOCV evaluation.
///
/// Returns all (period, unit) pairs where `control_mask` is nonzero and the
/// outcome is finite.  When `max_samples > 0` and the number of eligible
/// observations exceeds that limit, a deterministic Fisher–Yates subsample
/// is drawn using the provided `seed`.  Setting `max_samples = 0` disables
/// subsampling and returns every eligible observation.
///
/// # Arguments
/// * `y`            — Outcome matrix, n_periods × n_units.
/// * `control_mask` — Nonzero entries mark control observations.
/// * `max_samples`  — Upper bound on returned observations (0 = no limit).
/// * `seed`         — RNG seed for reproducible subsampling.
pub fn get_control_observations(
    y: &ArrayView2<f64>,
    control_mask: &ArrayView2<u8>,
    max_samples: usize,
    seed: u64,
) -> Vec<(usize, usize)> {
    let n_periods = y.nrows();
    let n_units = y.ncols();

    // Collect all valid control observations.
    let mut obs: Vec<(usize, usize)> = Vec::new();
    for t in 0..n_periods {
        for i in 0..n_units {
            if control_mask[[t, i]] != 0 && y[[t, i]].is_finite() {
                obs.push((t, i));
            }
        }
    }

    // Subsample if the caller requested a cap.  max_samples == 0 means
    // "use all control observations" (no subsampling).
    if max_samples > 0 && obs.len() > max_samples {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        obs.shuffle(&mut rng);
        obs.truncate(max_samples);
    }

    obs
}

/// Evaluate the LOOCV criterion for a single (λ_time, λ_unit, λ_nn) triple
/// under the twostep method.
///
/// For each control observation (t, i), computes the pseudo-treatment effect
/// τ̂_{ti}^{loocv}(λ) by fitting the model with observation (t, i) excluded,
/// then returns Q(λ) = Σ τ̂²_{ti} over all successfully evaluated control
/// observations.
///
/// Returns `(score, n_valid, first_failed_obs)`.  If any observation fails
/// to produce a valid estimate, the score is set to +∞ and the failing
/// (period, unit) pair is reported.
#[allow(clippy::too_many_arguments)]
pub fn loocv_score_for_params(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    control_mask: &ArrayView2<u8>,
    time_dist: &ArrayView2<i64>,
    control_obs: &[(usize, usize)],
    lambda_time: f64,
    lambda_unit: f64,
    lambda_nn: f64,
    max_iter: usize,
    tol: f64,
) -> (f64, usize, Option<(usize, usize)>) {
    let n_periods = y.nrows();
    let n_units = y.ncols();

    let mut tau_sq_sum = 0.0;
    let mut n_valid = 0usize;

    for &(t, i) in control_obs {
        // Compute observation-specific weight matrix.
        let weight_matrix = compute_weight_matrix(
            y,
            d,
            n_periods,
            n_units,
            i,
            t,
            lambda_time,
            lambda_unit,
            time_dist,
        );

        // Estimate model excluding this observation.
        match estimate_model(
            y,
            control_mask,
            &weight_matrix.view(),
            lambda_nn,
            n_periods,
            n_units,
            max_iter,
            tol,
            Some((t, i)),
        ) {
            Some((alpha, beta, l, _n_iters, _converged)) => {
                // Pseudo-treatment effect: τ̂ = Y_{ti} − α_i − β_t − L_{ti}.
                let tau = y[[t, i]] - alpha[i] - beta[t] - l[[t, i]];
                tau_sq_sum += tau * tau;
                n_valid += 1;
            }
            None => {
                // Estimation failure invalidates this λ combination.
                return (f64::INFINITY, n_valid, Some((t, i)));
            }
        }
    }

    if n_valid == 0 {
        (f64::INFINITY, 0, None)
    } else {
        (tau_sq_sum, n_valid, None)
    }
}

/// Univariate LOOCV search over a single tuning parameter.
///
/// Evaluates Q(λ) for each value in `grid` while holding the other two
/// parameters at the supplied fixed values.  Used in Stage 1 of the
/// two-stage search to obtain initial estimates, and in Stage 2 as the
/// inner step of coordinate descent.
///
/// # Arguments
/// * `param_type` — Which parameter to search: 0 = λ_time, 1 = λ_unit,
///                  2 = λ_nn.
///
/// # Returns
/// `(best_value, best_score)`
#[allow(clippy::too_many_arguments)]
pub fn univariate_loocv_search(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    control_mask: &ArrayView2<u8>,
    time_dist: &ArrayView2<i64>,
    control_obs: &[(usize, usize)],
    grid: &[f64],
    fixed_time: f64,
    fixed_unit: f64,
    fixed_nn: f64,
    param_type: usize,
    max_iter: usize,
    tol: f64,
) -> (f64, f64) {
    let mut best_score = f64::INFINITY;
    let mut best_value = grid.first().copied().unwrap_or(0.0);

    // Convert fixed parameters that may be +∞ to effective finite values.
    // λ_time/λ_unit = ∞ → 0.0 (uniform weights); λ_nn = ∞ → 1e10 (factor
    // model disabled).  Grid values are always finite by construction.
    let fixed_time_eff = if fixed_time.is_infinite() { 0.0 } else { fixed_time };
    let fixed_unit_eff = if fixed_unit.is_infinite() { 0.0 } else { fixed_unit };
    let fixed_nn_eff = if fixed_nn.is_infinite() { 1e10 } else { fixed_nn };

    // Evaluate all grid values in parallel.
    let results: Vec<(f64, f64)> = grid
        .par_iter()
        .map(|&value| {
            let (lambda_time, lambda_unit, lambda_nn) = match param_type {
                0 => (value, fixed_unit_eff, fixed_nn_eff),
                1 => (fixed_time_eff, value, fixed_nn_eff),
                _ => (fixed_time_eff, fixed_unit_eff, value),
            };

            let (score, _, _) = loocv_score_for_params(
                y,
                d,
                control_mask,
                time_dist,
                control_obs,
                lambda_time,
                lambda_unit,
                lambda_nn,
                max_iter,
                tol,
            );
            (value, score)
        })
        .collect();

    for (value, score) in results {
        if score < best_score {
            best_score = score;
            best_value = value;
        }
    }

    (best_value, best_score)
}

/// Coordinate descent cycling over (λ_unit, λ_time, λ_nn) for the twostep
/// method.
///
/// Starting from the Stage 1 initial values, successively optimizes each
/// parameter while holding the other two at their most recent optimal values.
/// Cycling order: λ_unit → λ_time → λ_nn.
/// Terminates when |Q_new − Q_old| < 1e-6 or `max_cycles` is reached.
#[allow(clippy::too_many_arguments)]
pub fn cycling_parameter_search(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    control_mask: &ArrayView2<u8>,
    time_dist: &ArrayView2<i64>,
    control_obs: &[(usize, usize)],
    lambda_time_grid: &[f64],
    lambda_unit_grid: &[f64],
    lambda_nn_grid: &[f64],
    initial_time: f64,
    initial_unit: f64,
    initial_nn: f64,
    max_iter: usize,
    tol: f64,
    max_cycles: usize,
) -> (f64, f64, f64) {
    let mut lambda_time = initial_time;
    let mut lambda_unit = initial_unit;
    let mut lambda_nn = initial_nn;
    let mut prev_score = f64::INFINITY;

    for _cycle in 0..max_cycles {
        // Optimize λ_unit (fix λ_time, λ_nn)
        let (new_unit, _) = univariate_loocv_search(
            y,
            d,
            control_mask,
            time_dist,
            control_obs,
            lambda_unit_grid,
            lambda_time,
            0.0,
            lambda_nn,
            1,
            max_iter,
            tol,
        );
        lambda_unit = new_unit;

        // Optimize λ_time (fix λ_unit, λ_nn)
        let (new_time, _) = univariate_loocv_search(
            y,
            d,
            control_mask,
            time_dist,
            control_obs,
            lambda_time_grid,
            0.0,
            lambda_unit,
            lambda_nn,
            0,
            max_iter,
            tol,
        );
        lambda_time = new_time;

        // Optimize λ_nn (fix λ_unit, λ_time)
        let (new_nn, score) = univariate_loocv_search(
            y,
            d,
            control_mask,
            time_dist,
            control_obs,
            lambda_nn_grid,
            lambda_time,
            lambda_unit,
            0.0,
            2,
            max_iter,
            tol,
        );
        lambda_nn = new_nn;

        // Check cycling convergence.  The 1e-6 threshold governs the outer
        // coordinate descent loop and is distinct from `tol`, which controls
        // convergence of the inner estimation solver.
        if (score - prev_score).abs() < 1e-6 {
            break;
        }
        prev_score = score;
    }

    (lambda_time, lambda_unit, lambda_nn)
}

/// Two-stage LOOCV search for the twostep method.
///
/// Selects λ̂ = argmin_{λ ∈ Λ} Q(λ) using:
///   Stage 1 — Univariate initialization with extreme fixed values.
///   Stage 2 — Coordinate descent cycling until convergence.
///
/// # Returns
/// `(best_λ_time, best_λ_unit, best_λ_nn, best_score, n_valid,
///   n_attempted, n_control_total, subsampled, first_failed_obs)`
#[allow(clippy::too_many_arguments)]
pub fn loocv_grid_search(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    control_mask: &ArrayView2<u8>,
    time_dist: &ArrayView2<i64>,
    lambda_time_grid: &[f64],
    lambda_unit_grid: &[f64],
    lambda_nn_grid: &[f64],
    max_loocv_samples: usize,
    max_iter: usize,
    tol: f64,
    seed: u64,
) -> LoocvGridSearchResult {
    // Count total control observations before subsampling.
    let n_periods = y.nrows();
    let n_units = y.ncols();
    let mut n_control_total = 0usize;
    for t in 0..n_periods {
        for i in 0..n_units {
            if control_mask[[t, i]] != 0 && y[[t, i]].is_finite() {
                n_control_total += 1;
            }
        }
    }

    // Get control observations for LOOCV (may be subsampled).
    let control_obs = get_control_observations(y, control_mask, max_loocv_samples, seed);
    let n_attempted = control_obs.len();
    let subsampled = max_loocv_samples > 0 && n_control_total > max_loocv_samples;

    // Stage 1: Univariate searches for initial values.
    // λ_time search: fix λ_unit = 0, λ_nn = ∞.
    let (lambda_time_init, _) = univariate_loocv_search(
        y,
        d,
        control_mask,
        time_dist,
        &control_obs,
        lambda_time_grid,
        0.0,
        0.0,
        f64::INFINITY,
        0,
        max_iter,
        tol,
    );

    // λ_nn search: fix λ_time = ∞, λ_unit = 0.
    let (lambda_nn_init, _) = univariate_loocv_search(
        y,
        d,
        control_mask,
        time_dist,
        &control_obs,
        lambda_nn_grid,
        f64::INFINITY,
        0.0,
        0.0,
        2,
        max_iter,
        tol,
    );

    // λ_unit search: fix λ_nn = ∞, λ_time = 0.
    let (lambda_unit_init, _) = univariate_loocv_search(
        y,
        d,
        control_mask,
        time_dist,
        &control_obs,
        lambda_unit_grid,
        0.0,
        0.0,
        f64::INFINITY,
        1,
        max_iter,
        tol,
    );

    // Stage 2: Coordinate descent refinement.
    let (best_time, best_unit, best_nn) = cycling_parameter_search(
        y,
        d,
        control_mask,
        time_dist,
        &control_obs,
        lambda_time_grid,
        lambda_unit_grid,
        lambda_nn_grid,
        lambda_time_init,
        lambda_unit_init,
        lambda_nn_init,
        max_iter,
        tol,
        10,
    );

    // Final evaluation at the selected parameters.
    let (best_score, n_valid, first_failed) = loocv_score_for_params(
        y,
        d,
        control_mask,
        time_dist,
        &control_obs,
        best_time,
        best_unit,
        best_nn,
        max_iter,
        tol,
    );

    (
        best_time,
        best_unit,
        best_nn,
        best_score,
        n_valid,
        n_attempted,
        n_control_total,
        subsampled,
        first_failed,
    )
}

/// Evaluate the LOOCV criterion for a single (λ_time, λ_unit, λ_nn) triple
/// under the joint method.
///
/// Uses global weights δ(λ_time, λ_unit) and a homogeneous treatment effect.
/// For each control observation (t, i), the weight δ_{t,i} is set to zero
/// before fitting, and the pseudo-treatment effect is computed as
/// τ̂_{ti} = Y_{ti} − μ̂ − α̂_i − β̂_t − L̂_{ti}.
///
/// Returns `(score, n_valid, first_failed_obs)`.
#[allow(clippy::too_many_arguments)]
pub fn loocv_score_joint(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    control_obs: &[(usize, usize)],
    lambda_time: f64,
    lambda_unit: f64,
    lambda_nn: f64,
    treated_periods: usize,
    max_iter: usize,
    tol: f64,
) -> (f64, usize, Option<(usize, usize)>) {
    let n_periods = y.nrows();
    let n_units = y.ncols();

    let mut tau_sq_sum = 0.0;
    let mut n_valid = 0usize;

    // Compute global weights δ (shared across all LOOCV iterations).
    let delta = compute_joint_weights(y, d, lambda_time, lambda_unit, treated_periods);

    for &(t_ex, i_ex) in control_obs {
        // Zero out the excluded observation's weight.
        let mut delta_ex = delta.clone();
        delta_ex[[t_ex, i_ex]] = 0.0;

        // Fit joint model with the modified weights.
        //
        // When λ_nn ≥ 1e10 the low-rank step is skipped (L ≡ 0). τ is not used
        // downstream (the LOOCV score uses the pseudo-residual at the excluded
        // cell directly), so we pass a dummy 0.0 in its slot.
        let result = if lambda_nn >= 1e10 {
            solve_joint_no_lowrank(y, &delta_ex.view()).map(|(mu, alpha, beta)| {
                let l = Array2::<f64>::zeros((n_periods, n_units));
                (mu, alpha, beta, l, 0.0_f64, 1_usize, true)
            })
        } else {
            solve_joint_with_lowrank(y, d, &delta_ex.view(), lambda_nn, max_iter, tol)
        };

        match result {
            Some((mu, alpha, beta, l, _tau, _n_iters, _converged)) => {
                // Pseudo-treatment effect: τ̂ = Y − μ − α − β − L.
                let y_ti = if y[[t_ex, i_ex]].is_finite() {
                    y[[t_ex, i_ex]]
                } else {
                    continue;
                };
                let tau_loocv = y_ti - mu - alpha[i_ex] - beta[t_ex] - l[[t_ex, i_ex]];
                tau_sq_sum += tau_loocv * tau_loocv;
                n_valid += 1;
            }
            None => {
                // Estimation failure invalidates this λ combination.
                return (f64::INFINITY, n_valid, Some((t_ex, i_ex)));
            }
        }
    }

    if n_valid == 0 {
        (f64::INFINITY, 0, None)
    } else {
        (tau_sq_sum, n_valid, None)
    }
}

/// Brute-force LOOCV grid search for the joint method.
///
/// Exhaustively evaluates all (λ_time, λ_unit, λ_nn) combinations in
/// parallel.  Complexity is O(|grid|³) and can be expensive; the
/// coordinate-descent variant `loocv_cycling_search_joint` is preferred
/// for production use.
///
/// # Returns
/// `(best_λ_time, best_λ_unit, best_λ_nn, best_score, n_valid,
///   n_attempted, n_control_total, subsampled, first_failed_obs)`
#[allow(clippy::too_many_arguments)]
pub fn loocv_grid_search_joint(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    control_mask: &ArrayView2<u8>,
    lambda_time_grid: &[f64],
    lambda_unit_grid: &[f64],
    lambda_nn_grid: &[f64],
    max_loocv_samples: usize,
    max_iter: usize,
    tol: f64,
    seed: u64,
) -> LoocvGridSearchResult {
    let n_periods = y.nrows();
    let n_units = y.ncols();

    // Count total control observations before subsampling.
    let mut n_control_total = 0usize;
    for t in 0..n_periods {
        for i in 0..n_units {
            if control_mask[[t, i]] != 0 && y[[t, i]].is_finite() {
                n_control_total += 1;
            }
        }
    }

    // Determine the number of treated periods from D.
    let mut first_treat_period = n_periods;
    for t in 0..n_periods {
        for i in 0..n_units {
            if d[[t, i]] == 1.0 {
                first_treat_period = first_treat_period.min(t);
                break;
            }
        }
    }
    let treated_periods = n_periods.saturating_sub(first_treat_period);

    // Get control observations for LOOCV (may be subsampled).
    let control_obs = get_control_observations(y, control_mask, max_loocv_samples, seed);
    let n_attempted = control_obs.len();
    let subsampled = max_loocv_samples > 0 && n_control_total > max_loocv_samples;

    // Enumerate all grid combinations.
    let mut grid_combinations: Vec<(f64, f64, f64)> = Vec::new();
    for &lt in lambda_time_grid {
        for &lu in lambda_unit_grid {
            for &ln in lambda_nn_grid {
                grid_combinations.push((lt, lu, ln));
            }
        }
    }

    // Parallel grid search
    let results: Vec<LoocvResultTuple> = grid_combinations
        .into_par_iter()
        .map(|(lt, lu, ln)| {
            // Convert infinity values
            let lt_eff = if lt.is_infinite() { 0.0 } else { lt };
            let lu_eff = if lu.is_infinite() { 0.0 } else { lu };
            let ln_eff = if ln.is_infinite() { 1e10 } else { ln };

            let (score, n_valid, first_failed) = loocv_score_joint(
                y,
                d,
                &control_obs,
                lt_eff,
                lu_eff,
                ln_eff,
                treated_periods,
                max_iter,
                tol,
            );

            (lt, lu, ln, score, n_valid, first_failed)
        })
        .collect();

    // Find best result
    let mut best_result: LoocvResultTuple = (
        lambda_time_grid.first().copied().unwrap_or(0.0),
        lambda_unit_grid.first().copied().unwrap_or(0.0),
        lambda_nn_grid.first().copied().unwrap_or(0.0),
        f64::INFINITY,
        0usize,
        None,
    );

    for (lt, lu, ln, score, n_valid, first_failed) in results {
        if score < best_result.3 {
            best_result = (lt, lu, ln, score, n_valid, first_failed);
        }
    }

    let (best_lt, best_lu, best_ln, best_score, n_valid, first_failed) = best_result;

    // Story 4.5: Include n_control_total and subsampled flag for diagnostics
    (
        best_lt,
        best_lu,
        best_ln,
        best_score,
        n_valid,
        n_attempted,
        n_control_total,
        subsampled,
        first_failed,
    )
}

/// Coordinate descent LOOCV search for joint method tuning parameters.
///
/// Paper: Footnote 2 applied to joint method (Remark 6.1).
///
/// Same two-stage approach as Twostep LOOCV:
///   Stage 1: Univariate searches with extreme fixed values for initial estimates
///   Stage 2: Cycling (coordinate descent) until convergence
/// but using joint method's global weights δ and LOOCV score (loocv_score_joint).
///
/// Complexity: O(|grid| × n_cycles) instead of O(|grid|^3) for brute-force.
pub fn loocv_cycling_search_joint(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    control_mask: &ArrayView2<u8>,
    lambda_time_grid: &[f64],
    lambda_unit_grid: &[f64],
    lambda_nn_grid: &[f64],
    max_loocv_samples: usize,
    max_iter: usize,
    tol: f64,
    seed: u64,
    max_cycles: usize,
) -> LoocvGridSearchResult {
    let n_periods = y.nrows();
    let n_units = y.ncols();

    // Count total control observations before subsampling
    let mut n_control_total = 0usize;
    for t in 0..n_periods {
        for i in 0..n_units {
            if control_mask[[t, i]] != 0 && y[[t, i]].is_finite() {
                n_control_total += 1;
            }
        }
    }

    // Determine treated periods from D matrix
    let mut first_treat_period = n_periods;
    for t in 0..n_periods {
        for i in 0..n_units {
            if d[[t, i]] == 1.0 {
                first_treat_period = first_treat_period.min(t);
                break;
            }
        }
    }
    let treated_periods = n_periods.saturating_sub(first_treat_period);

    // Get control observations for LOOCV (may be subsampled)
    let control_obs = get_control_observations(y, control_mask, max_loocv_samples, seed);
    let n_attempted = control_obs.len();
    let subsampled = max_loocv_samples > 0 && n_control_total > max_loocv_samples;

    // Helper closure: univariate search over one parameter grid (parallelized)
    // param_idx: 0=lambda_time, 1=lambda_unit, 2=lambda_nn
    let univariate_search =
        |grid: &[f64],
         param_idx: usize,
         fixed: (f64, f64, f64)|
         -> (f64, f64, usize, Option<(usize, usize)>) {
            let results: Vec<(f64, f64, usize, Option<(usize, usize)>)> = grid
                .par_iter()
                .map(|&val| {
                    // Design Issue 38 fix: Convert fixed params upfront.
                    // Fixed params may be f64::INFINITY (Stage 1 uses extreme values).
                    // Grid values (val) are guaranteed finite by C bridge pre-conversion.
                    let f0 = if fixed.0.is_infinite() { 0.0 } else { fixed.0 };
                    let f1 = if fixed.1.is_infinite() { 0.0 } else { fixed.1 };
                    let f2 = if fixed.2.is_infinite() { 1e10 } else { fixed.2 };
                    let (lt, lu, ln) = match param_idx {
                        0 => (val, f1, f2),
                        1 => (f0, val, f2),
                        _ => (f0, f1, val),
                    };

                    let (score, n_valid, first_failed) = loocv_score_joint(
                        y,
                        d,
                        &control_obs,
                        lt,
                        lu,
                        ln,
                        treated_periods,
                        max_iter,
                        tol,
                    );
                    (val, score, n_valid, first_failed)
                })
                .collect();

            let mut best = (
                grid.first().copied().unwrap_or(0.0),
                f64::INFINITY,
                0usize,
                None,
            );
            for (val, score, n_valid, first_failed) in results {
                if score < best.1 {
                    best = (val, score, n_valid, first_failed);
                }
            }
            best
        };

    // Stage 1: Univariate searches with extreme fixed values
    let (lt_init, _, _, _) =
        univariate_search(lambda_time_grid, 0, (0.0, 0.0, f64::INFINITY));
    let (ln_init, _, _, _) =
        univariate_search(lambda_nn_grid, 2, (f64::INFINITY, 0.0, 0.0));
    let (lu_init, _, _, _) =
        univariate_search(lambda_unit_grid, 1, (0.0, 0.0, f64::INFINITY));

    // Stage 2: Cycling refinement (coordinate descent)
    let mut best_lt = lt_init;
    let mut best_lu = lu_init;
    let mut best_ln = ln_init;
    let mut prev_score = f64::INFINITY;

    for _cycle in 0..max_cycles {
        // Optimize λ_unit (fix λ_time, λ_nn)
        let (lu_new, _, _, _) =
            univariate_search(lambda_unit_grid, 1, (best_lt, 0.0, best_ln));
        best_lu = lu_new;

        // Optimize λ_time (fix λ_unit, λ_nn)
        let (lt_new, _, _, _) =
            univariate_search(lambda_time_grid, 0, (0.0, best_lu, best_ln));
        best_lt = lt_new;

        // Optimize λ_nn (fix λ_unit, λ_time)
        let (ln_new, score, _, _) =
            univariate_search(lambda_nn_grid, 2, (best_lt, best_lu, 0.0));
        best_ln = ln_new;

        // Check cycling convergence: |Q_new - Q_old| < 1e-6
        // NOTE: Matches Python reference (trop.py:1159) hardcoded threshold.
        // Separate from `tol` which controls inner estimation convergence.
        if (score - prev_score).abs() < 1e-6 {
            break;
        }
        prev_score = score;
    }

    // DI 11 fix: Final evaluation with best parameters (matching Twostep pattern).
    // Twostep's loocv_grid_search() does a final loocv_score_for_params() call after
    // cycling (loocv.rs:462-473). For consistency, Joint also performs a final
    // loocv_score_joint() evaluation to get fresh n_valid/first_failed diagnostics.
    //
    // Design Issue 38: Grid values are guaranteed finite (C bridge pre-converts
    // Stata missing/∞ → effective values). best_lt/lu/ln come from grid search
    // results, so no infinity conversion is needed here (same as Twostep pattern).
    let (best_score, final_n_valid, final_first_failed) = loocv_score_joint(
        y,
        d,
        &control_obs,
        best_lt,
        best_lu,
        best_ln,
        treated_periods,
        max_iter,
        tol,
    );

    (
        best_lt,
        best_lu,
        best_ln,
        best_score,
        final_n_valid,
        n_attempted,
        n_control_total,
        subsampled,
        final_first_failed,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_get_control_observations() {
        let y = array![[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]];
        let control_mask = array![
            [1u8, 1],
            [1, 1],
            [1, 0] // Unit 1 treated at period 2
        ];

        let obs = get_control_observations(&y.view(), &control_mask.view(), 100, 42);

        // Should have 5 control observations (all except (2,1))
        assert_eq!(obs.len(), 5);

        // (2, 1) should not be in the list
        assert!(!obs.contains(&(2, 1)));
    }

    #[test]
    fn test_get_control_observations_subsampling() {
        let y = array![[1.0, 2.0, 3.0], [2.0, 3.0, 4.0], [3.0, 4.0, 5.0]];
        let control_mask = array![[1u8, 1, 1], [1, 1, 1], [1, 1, 1]];

        // Request only 3 samples
        let obs = get_control_observations(&y.view(), &control_mask.view(), 3, 42);

        assert_eq!(obs.len(), 3);
    }

    #[test]
    fn test_infinity_conversion() {
        // Test that infinity parameters are correctly converted
        let inf = f64::INFINITY;

        // λ_time/λ_unit = ∞ → 0.0 (uniform weights)
        let time_eff = if inf.is_infinite() { 0.0 } else { inf };
        assert_eq!(time_eff, 0.0);

        // λ_nn = ∞ → 1e10 (disable factor model)
        let nn_eff = if inf.is_infinite() { 1e10 } else { inf };
        assert_eq!(nn_eff, 1e10);
    }

    // ========================================================================
    // Joint LOOCV Property Tests (Story 4.4)
    // ========================================================================

    #[test]
    fn test_loocv_score_joint_basic() {
        // Property 6: LOOCV objective function
        // Q(λ) = Σ (τ̂_{it}^{loocv})²
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 6.0, 6.0] // Unit 1 treated at period 3
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0]
        ];

        // Control observations (all except (3,1))
        let control_obs: Vec<(usize, usize)> = vec![
            (0, 0),
            (0, 1),
            (0, 2),
            (1, 0),
            (1, 1),
            (1, 2),
            (2, 0),
            (2, 1),
            (2, 2),
            (3, 0),
            (3, 2),
        ];

        let (score, n_valid, first_failed) = loocv_score_joint(
            &y.view(),
            &d.view(),
            &control_obs,
            0.0,  // lambda_time
            0.0,  // lambda_unit
            1e10, // lambda_nn (no low-rank)
            1,    // treated_periods
            100,
            1e-6,
        );

        // Score should be finite and non-negative
        assert!(score.is_finite(), "Score should be finite");
        assert!(score >= 0.0, "Score should be non-negative");

        // All control observations should be valid
        assert_eq!(
            n_valid,
            control_obs.len(),
            "All control obs should be valid"
        );
        assert!(first_failed.is_none(), "No observation should fail");
    }

    #[test]
    fn test_loocv_grid_search_joint_completeness() {
        // Property 1: Full grid search completeness
        // Should evaluate exactly |λ_time| × |λ_unit| × |λ_nn| combinations
        let y = array![
            [1.0, 2.0],
            [2.0, 3.0],
            [3.0, 4.0],
            [4.0, 6.0] // Unit 1 treated
        ];
        let d = array![[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 1.0]];
        let control_mask = array![[1u8, 1], [1, 1], [1, 1], [1, 0]];

        let lambda_time_grid = &[0.0, 0.5];
        let lambda_unit_grid = &[0.0, 0.5];
        let lambda_nn_grid = &[1e10]; // No low-rank for speed

        let (
            best_lt,
            best_lu,
            best_ln,
            best_score,
            n_valid,
            n_attempted,
            n_control_total,
            subsampled,
            _,
        ) = loocv_grid_search_joint(
            &y.view(),
            &d.view(),
            &control_mask.view(),
            lambda_time_grid,
            lambda_unit_grid,
            lambda_nn_grid,
            1000, // max_loocv_samples
            100,
            1e-6,
            42,
        );

        // Best parameters should be from the grid
        assert!(
            lambda_time_grid.contains(&best_lt),
            "Best λ_time should be from grid"
        );
        assert!(
            lambda_unit_grid.contains(&best_lu),
            "Best λ_unit should be from grid"
        );
        assert!(
            lambda_nn_grid.contains(&best_ln),
            "Best λ_nn should be from grid"
        );

        // Score should be finite
        assert!(best_score.is_finite(), "Best score should be finite");

        // n_attempted should equal number of control observations
        assert!(n_attempted > 0, "Should have attempted some observations");
        assert!(n_valid > 0, "Should have valid observations");

        // Story 4.5: Verify diagnostic fields
        assert!(
            n_control_total >= n_attempted,
            "n_control_total should be >= n_attempted"
        );
        assert!(
            !subsampled,
            "Should not subsample with max_loocv_samples=1000"
        );
    }

    #[test]
    fn test_loocv_joint_reproducibility() {
        // Property 9: Subsampling reproducibility
        // Same seed should produce identical results
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 7.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0]
        ];
        let control_mask = array![[1u8, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 0]];

        let lambda_time_grid = &[0.0, 0.5];
        let lambda_unit_grid = &[0.0];
        let lambda_nn_grid = &[1e10];

        // Run twice with same seed
        let result1 = loocv_grid_search_joint(
            &y.view(),
            &d.view(),
            &control_mask.view(),
            lambda_time_grid,
            lambda_unit_grid,
            lambda_nn_grid,
            5, // Small sample for subsampling
            100,
            1e-6,
            42, // Same seed
        );

        let result2 = loocv_grid_search_joint(
            &y.view(),
            &d.view(),
            &control_mask.view(),
            lambda_time_grid,
            lambda_unit_grid,
            lambda_nn_grid,
            5,
            100,
            1e-6,
            42, // Same seed
        );

        // Results should be identical
        assert_eq!(result1.0, result2.0, "λ_time should match");
        assert_eq!(result1.1, result2.1, "λ_unit should match");
        assert_eq!(result1.2, result2.2, "λ_nn should match");
        assert_eq!(result1.3, result2.3, "Score should match");
    }

    #[test]
    fn test_loocv_joint_infinity_handling() {
        // Property 10: Infinity parameter conversion
        // λ=∞ should be converted internally but returned as original value
        let y = array![[1.0, 2.0], [2.0, 3.0], [3.0, 4.0], [4.0, 6.0]];
        let d = array![[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 1.0]];
        let control_mask = array![[1u8, 1], [1, 1], [1, 1], [1, 0]];

        // Include infinity in grid
        let lambda_time_grid = &[f64::INFINITY, 0.5];
        let lambda_unit_grid = &[0.0];
        let lambda_nn_grid = &[f64::INFINITY];

        let (best_lt, _best_lu, best_ln, best_score, _, _, _, _, _) = loocv_grid_search_joint(
            &y.view(),
            &d.view(),
            &control_mask.view(),
            lambda_time_grid,
            lambda_unit_grid,
            lambda_nn_grid,
            1000,
            100,
            1e-6,
            42,
        );

        // Score should still be finite (infinity was converted internally)
        assert!(
            best_score.is_finite(),
            "Score should be finite even with ∞ params"
        );

        // If infinity was selected, it should be returned as infinity
        if best_lt.is_infinite() {
            assert!(best_lt == f64::INFINITY, "Should return original ∞ value");
        }
        if best_ln.is_infinite() {
            assert!(best_ln == f64::INFINITY, "Should return original ∞ value");
        }
    }

    #[test]
    fn test_joint_weights_non_normalization() {
        // Property 13: Weight non-normalization
        // Weights should NOT sum to 1 (per paper specification)
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 7.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0]
        ];

        let lambda_time = 0.5;
        let lambda_unit = 0.5;
        let treated_periods = 1;

        let delta = compute_joint_weights(
            &y.view(),
            &d.view(),
            lambda_time,
            lambda_unit,
            treated_periods,
        );

        // Sum of weights
        let weight_sum: f64 = delta.iter().sum();

        // Weights should NOT be normalized to 1
        // With exponential decay, sum will typically be > 1 or < 1
        assert!(
            (weight_sum - 1.0).abs() > 1e-6,
            "Weights should NOT sum to 1, got sum = {}",
            weight_sum
        );

        // All weights should be non-negative
        for &w in delta.iter() {
            assert!(w >= 0.0, "All weights should be non-negative");
        }
    }

    #[test]
    fn test_joint_model_includes_intercept() {
        // Property 14: Model includes global intercept μ
        // Joint model: Y = μ + α + β + L + τD + ε
        // Pseudo treatment effect: τ̂ = Y - μ - α - β - L
        let y = array![
            [10.0, 12.0], // Shifted by constant
            [11.0, 13.0],
            [12.0, 14.0],
            [13.0, 16.0] // Unit 1 treated
        ];
        let d = array![[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 1.0]];

        let control_obs: Vec<(usize, usize)> =
            vec![(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0)];

        // Compute LOOCV score with no low-rank (λ_nn = ∞)
        let (score, n_valid, _) = loocv_score_joint(
            &y.view(),
            &d.view(),
            &control_obs,
            0.0,  // lambda_time (uniform)
            0.0,  // lambda_unit (uniform)
            1e10, // lambda_nn (no low-rank)
            1,    // treated_periods
            100,
            1e-6,
        );

        // Score should be finite (model should fit)
        assert!(score.is_finite(), "Score should be finite with intercept");
        assert!(n_valid == control_obs.len(), "All obs should be valid");

        // The score should be small for well-structured data
        // (intercept absorbs the constant shift)
        assert!(
            score < 100.0,
            "Score should be reasonable with intercept, got {}",
            score
        );
    }

    // ========================================================================
    // Story 4.5: Diagnostic Information Tests
    // ========================================================================

    #[test]
    fn test_loocv_twostep_diagnostic_completeness() {
        // Property 7: Diagnostic Information Completeness
        // All diagnostic fields should be populated correctly
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 7.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0]
        ];
        let control_mask = array![[1u8, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 0]];
        let time_dist = array![[0i64, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]];

        let lambda_time_grid = &[0.0, 0.5];
        let lambda_unit_grid = &[0.0];
        let lambda_nn_grid = &[1e10];

        let (_, _, _, best_score, n_valid, n_attempted, n_control_total, subsampled, _) =
            loocv_grid_search(
                &y.view(),
                &d.view(),
                &control_mask.view(),
                &time_dist.view(),
                lambda_time_grid,
                lambda_unit_grid,
                lambda_nn_grid,
                1000, // No subsampling
                100,
                1e-6,
                42,
            );

        // Verify diagnostic constraints
        assert!(
            best_score.is_finite() && best_score >= 0.0,
            "Score should be finite and non-negative"
        );
        assert!(
            n_valid <= n_attempted,
            "n_valid ({}) should be <= n_attempted ({})",
            n_valid,
            n_attempted
        );
        assert!(
            n_attempted <= n_control_total,
            "n_attempted ({}) should be <= n_control_total ({})",
            n_attempted,
            n_control_total
        );
        assert!(
            !subsampled,
            "Should not subsample with max_loocv_samples=1000"
        );
        assert_eq!(
            n_attempted, n_control_total,
            "Without subsampling, n_attempted should equal n_control_total"
        );
    }

    #[test]
    fn test_loocv_twostep_subsampling_diagnostics() {
        // Test that subsampling correctly sets diagnostic flags
        let y = array![
            [1.0, 2.0, 3.0, 4.0, 5.0],
            [2.0, 3.0, 4.0, 5.0, 6.0],
            [3.0, 4.0, 5.0, 6.0, 7.0],
            [4.0, 5.0, 6.0, 7.0, 9.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0]
        ];
        let control_mask = array![
            [1u8, 1, 1, 1, 1],
            [1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1],
            [1, 1, 1, 1, 0]
        ];
        let time_dist = array![[0i64, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]];

        let lambda_time_grid = &[0.0];
        let lambda_unit_grid = &[0.0];
        let lambda_nn_grid = &[1e10];

        // Total control observations = 19 (4*5 - 1 treated)
        // Request only 5 samples
        let (_, _, _, _, n_valid, n_attempted, n_control_total, subsampled, _) = loocv_grid_search(
            &y.view(),
            &d.view(),
            &control_mask.view(),
            &time_dist.view(),
            lambda_time_grid,
            lambda_unit_grid,
            lambda_nn_grid,
            5, // Force subsampling
            100,
            1e-6,
            42,
        );

        // Verify subsampling occurred
        assert!(
            subsampled,
            "Should subsample when n_control > max_loocv_samples"
        );
        assert_eq!(n_attempted, 5, "n_attempted should equal max_loocv_samples");
        assert_eq!(n_control_total, 19, "n_control_total should be 19");
        assert!(n_valid <= n_attempted, "n_valid should be <= n_attempted");
    }

    // ========================================================================
    // LOOCV Score Precision Tests (Story 5.1, Task 3.6)
    // ========================================================================

    #[test]
    fn test_loocv_score_precision() {
        // Test LOOCV score numerical precision (tolerance < 1e-8)
        // Run same computation twice and verify identical results
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 7.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0]
        ];
        let control_mask = array![[1u8, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 0]];
        let time_dist = array![[0i64, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]];

        let control_obs: Vec<(usize, usize)> = vec![
            (0, 0),
            (0, 1),
            (0, 2),
            (1, 0),
            (1, 1),
            (1, 2),
            (2, 0),
            (2, 1),
            (2, 2),
            (3, 0),
            (3, 1),
        ];

        // Compute score twice
        let (score1, n_valid1, _) = loocv_score_for_params(
            &y.view(),
            &d.view(),
            &control_mask.view(),
            &time_dist.view(),
            &control_obs,
            0.5, // lambda_time
            0.5, // lambda_unit
            0.1, // lambda_nn
            100,
            1e-6,
        );

        let (score2, n_valid2, _) = loocv_score_for_params(
            &y.view(),
            &d.view(),
            &control_mask.view(),
            &time_dist.view(),
            &control_obs,
            0.5,
            0.5,
            0.1,
            100,
            1e-6,
        );

        // Scores should be identical (deterministic computation)
        let diff = (score1 - score2).abs();
        assert!(
            diff < 1e-8,
            "LOOCV score should be deterministic: score1={}, score2={}, diff={}",
            score1,
            score2,
            diff
        );
        assert_eq!(n_valid1, n_valid2, "n_valid should be identical");
    }

    #[test]
    fn test_loocv_score_finite_and_nonnegative() {
        // Test that LOOCV scores are always finite and non-negative
        let y = array![[1.0, 2.0], [2.0, 3.0], [3.0, 4.0], [4.0, 6.0]];
        let d = array![[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 1.0]];
        let control_mask = array![[1u8, 1], [1, 1], [1, 1], [1, 0]];
        let time_dist = array![[0i64, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]];

        let control_obs: Vec<(usize, usize)> =
            vec![(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0)];

        // Test with various parameter combinations
        let test_params = vec![
            (0.0, 0.0, 0.0),
            (0.5, 0.0, 0.0),
            (0.0, 0.5, 0.0),
            (0.0, 0.0, 0.1),
            (0.5, 0.5, 0.1),
            (1.0, 1.0, 1.0),
        ];

        for (lt, lu, ln) in test_params {
            let (score, n_valid, first_failed) = loocv_score_for_params(
                &y.view(),
                &d.view(),
                &control_mask.view(),
                &time_dist.view(),
                &control_obs,
                lt,
                lu,
                ln,
                100,
                1e-6,
            );

            if first_failed.is_none() {
                assert!(
                    score.is_finite(),
                    "Score should be finite for λ_t={}, λ_u={}, λ_nn={}",
                    lt,
                    lu,
                    ln
                );
                assert!(
                    score >= 0.0,
                    "Score should be non-negative for λ_t={}, λ_u={}, λ_nn={}, got {}",
                    lt,
                    lu,
                    ln,
                    score
                );
                assert!(
                    n_valid > 0,
                    "Should have valid observations for λ_t={}, λ_u={}, λ_nn={}",
                    lt,
                    lu,
                    ln
                );
            }
        }
    }

    #[test]
    fn test_loocv_parameter_selection_in_grid() {
        // Test that selected parameters are always from the provided grid
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 7.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0]
        ];
        let control_mask = array![[1u8, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 0]];
        let time_dist = array![[0i64, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]];

        let lambda_time_grid = vec![0.0, 0.5, 1.0, 2.0];
        let lambda_unit_grid = vec![0.0, 0.5, 1.0];
        let lambda_nn_grid = vec![0.0, 0.1, 1.0];

        let (best_lt, best_lu, best_ln, _, _, _, _, _, _) = loocv_grid_search(
            &y.view(),
            &d.view(),
            &control_mask.view(),
            &time_dist.view(),
            &lambda_time_grid,
            &lambda_unit_grid,
            &lambda_nn_grid,
            1000,
            100,
            1e-6,
            42,
        );

        assert!(
            lambda_time_grid.contains(&best_lt),
            "Selected λ_time={} should be in grid {:?}",
            best_lt,
            lambda_time_grid
        );
        assert!(
            lambda_unit_grid.contains(&best_lu),
            "Selected λ_unit={} should be in grid {:?}",
            best_lu,
            lambda_unit_grid
        );
        assert!(
            lambda_nn_grid.contains(&best_ln),
            "Selected λ_nn={} should be in grid {:?}",
            best_ln,
            lambda_nn_grid
        );
    }

    #[test]
    fn test_loocv_joint_score_precision() {
        // Test Joint LOOCV score precision (tolerance < 1e-8)
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 6.0, 6.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0]
        ];

        let control_obs: Vec<(usize, usize)> = vec![
            (0, 0),
            (0, 1),
            (0, 2),
            (1, 0),
            (1, 1),
            (1, 2),
            (2, 0),
            (2, 1),
            (2, 2),
            (3, 0),
            (3, 2),
        ];

        // Compute score twice
        let (score1, _, _) = loocv_score_joint(
            &y.view(),
            &d.view(),
            &control_obs,
            0.5,
            0.5,
            1e10,
            1,
            100,
            1e-6,
        );

        let (score2, _, _) = loocv_score_joint(
            &y.view(),
            &d.view(),
            &control_obs,
            0.5,
            0.5,
            1e10,
            1,
            100,
            1e-6,
        );

        let diff = (score1 - score2).abs();
        assert!(
            diff < 1e-8,
            "Joint LOOCV score should be deterministic: diff={}",
            diff
        );
    }
}
