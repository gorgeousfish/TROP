//! Bootstrap variance estimation for the TROP estimator.
//!
//! Implements unit-level block bootstrap with stratified resampling to preserve
//! the ratio of treated to control units.
//!
//! Procedure:
//!   1. For b = 1, ..., B: draw N_0 control units and N_1 treated units with replacement.
//!   2. Re-estimate the treatment effect on each bootstrap sample.
//!   3. Compute the population variance of the B bootstrap estimates.
//!
//! The PRNG is Xoshiro256PlusPlus, seeded deterministically per iteration.

use ndarray::{Array1, Array2, ArrayView2};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

use crate::estimation::{estimate_model, solve_joint_no_lowrank, solve_joint_with_lowrank};
use crate::weights::{compute_joint_weights, compute_weight_matrix};

// ============================================================================
// Unit Classification
// ============================================================================

/// Unit classification for stratified bootstrap sampling.
#[derive(Debug, Clone)]
pub struct UnitClassification {
    /// Column indices of control units (D[t,i] = 0 for all t).
    pub control_units: Vec<usize>,
    /// Column indices of treated units (D[t,i] = 1 for some t).
    pub treated_units: Vec<usize>,
    /// Number of control units (N_0).
    pub n_control: usize,
    /// Number of treated units (N_1).
    pub n_treated: usize,
}

/// Partition units into control and treated groups.
///
/// A unit is classified as treated if D[t,i] = 1 for any period t,
/// and as control otherwise.
///
/// # Arguments
/// * `d` - Treatment indicator matrix (T × N).
///
/// # Returns
/// A [`UnitClassification`] with separated index vectors.
pub fn classify_units(d: &ArrayView2<f64>) -> UnitClassification {
    let n_periods = d.nrows();
    let n_units = d.ncols();

    let mut control_units: Vec<usize> = Vec::new();
    let mut treated_units: Vec<usize> = Vec::new();

    for i in 0..n_units {
        let is_ever_treated = (0..n_periods).any(|t| d[[t, i]] == 1.0);
        if is_ever_treated {
            treated_units.push(i);
        } else {
            control_units.push(i);
        }
    }

    let n_control = control_units.len();
    let n_treated = treated_units.len();

    UnitClassification {
        control_units,
        treated_units,
        n_control,
        n_treated,
    }
}

// ============================================================================
// Stratified Sampling
// ============================================================================

/// Draw a stratified bootstrap sample for one iteration.
///
/// Independently resamples N_0 control units and N_1 treated units
/// with replacement, preserving the group ratio.
///
/// # Arguments
/// * `classification` - Pre-computed unit partition.
/// * `seed` - Deterministic seed for this iteration.
///
/// # Returns
/// Vector of sampled column indices (length N_0 + N_1).
pub fn stratified_sample(classification: &UnitClassification, seed: u64) -> Vec<usize> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut sampled_units: Vec<usize> =
        Vec::with_capacity(classification.n_control + classification.n_treated);

    // Sample control units with replacement
    for _ in 0..classification.n_control {
        if classification.n_control > 0 {
            let idx = rng.gen_range(0..classification.n_control);
            sampled_units.push(classification.control_units[idx]);
        }
    }

    // Sample treated units with replacement
    for _ in 0..classification.n_treated {
        if classification.n_treated > 0 {
            let idx = rng.gen_range(0..classification.n_treated);
            sampled_units.push(classification.treated_units[idx]);
        }
    }

    sampled_units
}

// ============================================================================
// Bootstrap Matrix Construction
// ============================================================================

/// Construct bootstrap outcome and treatment matrices by column selection.
///
/// Selects columns from Y and D according to `sampled_units`, preserving
/// all T rows per unit (block bootstrap at the unit level).
///
/// # Arguments
/// * `y` - Outcome matrix (T × N).
/// * `d` - Treatment indicator matrix (T × N).
/// * `sampled_units` - Column indices drawn by stratified sampling.
///
/// # Returns
/// `(Y_boot, D_boot)` — bootstrap matrices of dimension T × len(sampled_units).
pub fn build_bootstrap_matrices(
    y: &Array2<f64>,
    d: &Array2<f64>,
    sampled_units: &[usize],
) -> (Array2<f64>, Array2<f64>) {
    let n_periods = y.nrows();
    let n_units = sampled_units.len();

    let mut y_boot = Array2::<f64>::zeros((n_periods, n_units));
    let mut d_boot = Array2::<f64>::zeros((n_periods, n_units));

    for (new_idx, &old_idx) in sampled_units.iter().enumerate() {
        for t in 0..n_periods {
            y_boot[[t, new_idx]] = y[[t, old_idx]];
            d_boot[[t, new_idx]] = d[[t, old_idx]];
        }
    }

    (y_boot, d_boot)
}

/// Construct bootstrap matrices including the control-observation mask.
///
/// Same column-selection logic as [`build_bootstrap_matrices`], with an
/// additional binary mask indicating control observations (used by the
/// twostep estimator).
///
/// # Arguments
/// * `y` - Outcome matrix (T × N).
/// * `d` - Treatment indicator matrix (T × N).
/// * `control_mask` - Binary mask for control observations (T × N).
/// * `sampled_units` - Column indices drawn by stratified sampling.
///
/// # Returns
/// `(Y_boot, D_boot, control_mask_boot)`.
pub fn build_bootstrap_matrices_with_mask(
    y: &Array2<f64>,
    d: &Array2<f64>,
    control_mask: &Array2<u8>,
    sampled_units: &[usize],
) -> (Array2<f64>, Array2<f64>, Array2<u8>) {
    let n_periods = y.nrows();
    let n_units = sampled_units.len();

    let mut y_boot = Array2::<f64>::zeros((n_periods, n_units));
    let mut d_boot = Array2::<f64>::zeros((n_periods, n_units));
    let mut control_mask_boot = Array2::<u8>::zeros((n_periods, n_units));

    for (new_idx, &old_idx) in sampled_units.iter().enumerate() {
        for t in 0..n_periods {
            y_boot[[t, new_idx]] = y[[t, old_idx]];
            d_boot[[t, new_idx]] = d[[t, old_idx]];
            control_mask_boot[[t, new_idx]] = control_mask[[t, old_idx]];
        }
    }

    (y_boot, d_boot, control_mask_boot)
}

// ============================================================================
// Variance and CI Calculation
// ============================================================================

/// Compute the population variance and standard error of bootstrap estimates.
///
/// Uses the 1/B denominator (population variance), since the B replications
/// form the empirical sampling distribution of the estimator.
///
///   V = (1/B) Σ (τ_b − τ̄)²,   SE = √V
///
/// Non-finite values are silently discarded before computation.
///
/// # Arguments
/// * `estimates` - Bootstrap ATT estimates (may contain NaN/Inf).
///
/// # Returns
/// `(mean, variance, se)`. Returns `(0, 0, 0)` when fewer than two
/// finite values are available.
pub fn compute_bootstrap_variance(estimates: &[f64]) -> (f64, f64, f64) {
    // Filter non-finite values before computing statistics.
    let finite: Vec<f64> = estimates.iter().copied().filter(|x| x.is_finite()).collect();

    if finite.len() < 2 {
        let m = if finite.is_empty() { 0.0 } else { finite[0] };
        return (m, 0.0, 0.0);
    }

    let n = finite.len() as f64;
    let mean = finite.iter().sum::<f64>() / n;
    // Population variance (1/B denominator).
    let variance = finite.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    let se = variance.sqrt();

    (mean, variance, se)
}

/// Compute a percentile confidence interval from bootstrap estimates.
///
/// Returns the α/2 and 1−α/2 quantiles of the finite bootstrap values.
/// Fractional indices are resolved by linear interpolation:
/// for index f, the quantile is `sorted[⌊f⌋]·(1−frac) + sorted[⌈f⌉]·frac`.
///
/// # Arguments
/// * `estimates` - Bootstrap ATT estimates (may contain NaN/Inf).
/// * `alpha` - Significance level (e.g. 0.05 for a 95 % CI).
///
/// # Returns
/// `(ci_lower, ci_upper)`. Returns `(NaN, NaN)` when no finite values exist.
pub fn compute_percentile_ci(estimates: &[f64], alpha: f64) -> (f64, f64) {
    // Filter non-finite values before computing percentiles.
    let mut sorted: Vec<f64> = estimates.iter().copied().filter(|x| x.is_finite()).collect();
    if sorted.is_empty() {
        return (f64::NAN, f64::NAN);
    }

    // partial_cmp is total after NaN removal.
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("non-finite values already removed"));

    let n = sorted.len();

    // Percentile index: p-th quantile at position (n−1)·p.
    let lower_p = alpha / 2.0;
    let upper_p = 1.0 - alpha / 2.0;

    let lower_idx_f = (n - 1) as f64 * lower_p;
    let upper_idx_f = (n - 1) as f64 * upper_p;

    // Linear interpolation for fractional indices.
    let ci_lower = interpolate_percentile(&sorted, lower_idx_f);
    let ci_upper = interpolate_percentile(&sorted, upper_idx_f);

    (ci_lower, ci_upper)
}

/// Linearly interpolate a quantile from a sorted slice.
///
/// For fractional index f, returns
/// `sorted[⌊f⌋] * (1 − frac) + sorted[⌈f⌉] * frac`.
fn interpolate_percentile(sorted: &[f64], idx_f: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 {
        return sorted[0];
    }

    let idx_low = idx_f.floor() as usize;
    let idx_high = idx_f.ceil() as usize;

    let idx_low = idx_low.min(n - 1);
    let idx_high = idx_high.min(n - 1);

    if idx_low == idx_high {
        sorted[idx_low]
    } else {
        let frac = idx_f - idx_low as f64;
        sorted[idx_low] * (1.0 - frac) + sorted[idx_high] * frac
    }
}

// ============================================================================
// Bootstrap Result Structure
// ============================================================================

/// Aggregated bootstrap inference results.
#[derive(Debug, Clone)]
pub struct BootstrapResult {
    /// Finite bootstrap ATT estimates retained after filtering.
    pub estimates: Vec<f64>,
    /// Standard error (√ of population variance).
    pub se: f64,
    /// Mean of the bootstrap distribution.
    pub mean: f64,
    /// Lower bound of the percentile confidence interval.
    pub ci_lower: f64,
    /// Upper bound of the percentile confidence interval.
    pub ci_upper: f64,
    /// Number of iterations that produced a finite estimate.
    pub n_valid: usize,
    /// Total number of bootstrap iterations attempted.
    pub n_total: usize,
    /// Nominal confidence level (1 − α).
    pub level: f64,
}

// ============================================================================
// Main Bootstrap Functions
// ============================================================================

/// Bootstrap variance for the twostep estimator (convenience wrapper).
///
/// Delegates to [`bootstrap_trop_variance_full`] and returns only the
/// vector of bootstrap estimates and the standard error.
///
/// # Arguments
/// * `y` - Outcome matrix (T × N).
/// * `d` - Treatment indicator matrix (T × N).
/// * `control_mask` - Binary mask for control observations (T × N).
/// * `time_dist` - Pre-computed time distance matrix (T × T).
/// * `lambda_time` - Time kernel bandwidth (Inf → uniform weights).
/// * `lambda_unit` - Unit kernel bandwidth (Inf → uniform weights).
/// * `lambda_nn` - Nuclear norm penalty (Inf → no low-rank component).
/// * `n_bootstrap` - Number of bootstrap replications B.
/// * `max_iter` - Maximum ADMM iterations per estimation.
/// * `tol` - Convergence tolerance for ADMM.
/// * `seed` - Base random seed.
/// * `alpha` - Significance level for the percentile CI.
///
/// # Returns
/// `(estimates, se)` where `estimates` has length ≤ B.

#[allow(clippy::too_many_arguments)]
pub fn bootstrap_trop_variance(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    control_mask: &ArrayView2<u8>,
    time_dist: &ArrayView2<i64>,
    lambda_time: f64,
    lambda_unit: f64,
    lambda_nn: f64,
    n_bootstrap: usize,
    max_iter: usize,
    tol: f64,
    seed: u64,
    alpha: f64,
) -> (Array1<f64>, f64) {
    let result = bootstrap_trop_variance_full(
        y,
        d,
        control_mask,
        time_dist,
        lambda_time,
        lambda_unit,
        lambda_nn,
        n_bootstrap,
        max_iter,
        tol,
        seed,
        alpha,
    );
    (Array1::from_vec(result.estimates), result.se)
}

/// Bootstrap variance estimation for the twostep method with full output.
///
/// For each replication b = 1, …, B:
///   1. Draw a stratified sample of N_0 control + N_1 treated units.
///   2. Estimate per-observation treatment effects on the bootstrap data.
///   3. Average over treated observations to obtain τ̂_b.
///
/// Returns the complete [`BootstrapResult`] including percentile CI.
///
/// # Arguments
/// See [`bootstrap_trop_variance`] — all parameters are identical.
///
/// # Returns
/// A [`BootstrapResult`] containing estimates, SE, mean, and CI.
#[allow(clippy::too_many_arguments)]
pub fn bootstrap_trop_variance_full(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    control_mask: &ArrayView2<u8>,
    time_dist: &ArrayView2<i64>,
    lambda_time: f64,
    lambda_unit: f64,
    lambda_nn: f64,
    n_bootstrap: usize,
    max_iter: usize,
    tol: f64,
    seed: u64,
    alpha: f64,
) -> BootstrapResult {
    let y_arr = y.to_owned();
    let d_arr = d.to_owned();
    let control_mask_arr = control_mask.to_owned();
    let time_dist_arr = time_dist.to_owned();

    let n_periods = y_arr.nrows();

    // Map infinite lambda values to effective computation values:
    //   lambda_time/unit = Inf  →  0.0  (uniform kernel weights)
    //   lambda_nn        = Inf  →  1e10 (effectively no low-rank penalty)
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

    let classification = classify_units(&d_arr.view());

    // Parallel bootstrap: failed iterations yield None and are discarded.
    let bootstrap_estimates: Vec<f64> = (0..n_bootstrap)
        .into_par_iter()
        .filter_map(|b| {
            let iteration_seed = seed.wrapping_add(b as u64);
            let sampled_units = stratified_sample(&classification, iteration_seed);

            let (y_boot, d_boot, control_mask_boot) = build_bootstrap_matrices_with_mask(
                &y_arr,
                &d_arr,
                &control_mask_arr,
                &sampled_units,
            );

            let n_boot_units = d_boot.ncols();

            // Identify treated (t, i) pairs in the bootstrap sample.
            let mut boot_treated: Vec<(usize, usize)> = Vec::new();
            for t in 0..n_periods {
                for i in 0..n_boot_units {
                    if d_boot[[t, i]] == 1.0 {
                        boot_treated.push((t, i));
                    }
                }
            }

            if boot_treated.is_empty() {
                return None;
            }

            // Identify control units (never treated in bootstrap sample).
            let mut boot_control_units: Vec<usize> = Vec::new();
            for i in 0..n_boot_units {
                let is_control = (0..n_periods).all(|t| d_boot[[t, i]] == 0.0);
                if is_control {
                    boot_control_units.push(i);
                }
            }

            if boot_control_units.is_empty() {
                return None;
            }

            // Estimate τ(t,i) for each treated observation.
            let mut tau_values = Vec::with_capacity(boot_treated.len());

            for (t, i) in boot_treated {
                let weight_matrix = compute_weight_matrix(
                    &y_boot.view(),
                    &d_boot.view(),
                    n_periods,
                    n_boot_units,
                    i,
                    t,
                    lt_eff,
                    lu_eff,
                    &time_dist_arr.view(),
                );

                if let Some((alpha_est, beta, l, _n_iters, _converged)) = estimate_model(
                    &y_boot.view(),
                    &control_mask_boot.view(),
                    &weight_matrix.view(),
                    ln_eff,
                    n_periods,
                    n_boot_units,
                    max_iter,
                    tol,
                    None,
                ) {
                    let tau = y_boot[[t, i]] - alpha_est[i] - beta[t] - l[[t, i]];
                    if tau.is_finite() {
                        tau_values.push(tau);
                    }
                }
            }

            if tau_values.is_empty() {
                None
            } else {
                let att = tau_values.iter().sum::<f64>() / tau_values.len() as f64;
                if att.is_finite() { Some(att) } else { None }
            }
        })
        .collect();

    let n_valid = bootstrap_estimates.len();

    let (mean, _, se) = compute_bootstrap_variance(&bootstrap_estimates);
    let (ci_lower, ci_upper) = compute_percentile_ci(&bootstrap_estimates, alpha);

    BootstrapResult {
        estimates: bootstrap_estimates,
        se,
        mean,
        ci_lower,
        ci_upper,
        n_valid,
        n_total: n_bootstrap,
        level: 1.0 - alpha,
    }
}

/// Bootstrap variance estimation for the joint method (convenience wrapper).
///
/// Delegates to [`bootstrap_trop_variance_joint_full`] and returns only the
/// estimate vector and standard error.
///
/// # Arguments
/// * `y` - Outcome matrix (T × N).
/// * `d` - Treatment indicator matrix (T × N).
/// * `lambda_time` - Time kernel bandwidth (Inf → uniform weights).
/// * `lambda_unit` - Unit kernel bandwidth (Inf → uniform weights).
/// * `lambda_nn` - Nuclear norm penalty (Inf → no low-rank component).
/// * `n_bootstrap` - Number of bootstrap replications B.
/// * `max_iter` - Maximum ADMM iterations per estimation.
/// * `tol` - Convergence tolerance for ADMM.
/// * `seed` - Base random seed.
/// * `alpha` - Significance level for the percentile CI.
///
/// # Returns
/// `(estimates, se)` where `estimates` has length ≤ B.

#[allow(clippy::too_many_arguments)]
pub fn bootstrap_trop_variance_joint(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    lambda_time: f64,
    lambda_unit: f64,
    lambda_nn: f64,
    n_bootstrap: usize,
    max_iter: usize,
    tol: f64,
    seed: u64,
    alpha: f64,
) -> (Array1<f64>, f64) {
    let result = bootstrap_trop_variance_joint_full(
        y,
        d,
        lambda_time,
        lambda_unit,
        lambda_nn,
        n_bootstrap,
        max_iter,
        tol,
        seed,
        alpha,
    );
    (Array1::from_vec(result.estimates), result.se)
}

/// Bootstrap variance estimation for the joint method with full output.
///
/// For each replication b = 1, …, B:
///   1. Draw a stratified sample of N_0 control + N_1 treated units.
///   2. Compute global joint weights on the bootstrap data.
///   3. Solve the joint WLS problem to obtain τ̂_b.
///
/// Returns the complete [`BootstrapResult`] including percentile CI.
///
/// # Arguments
/// See [`bootstrap_trop_variance_joint`] — all parameters are identical.
///
/// # Returns
/// A [`BootstrapResult`] containing estimates, SE, mean, and CI.
#[allow(clippy::too_many_arguments)]
pub fn bootstrap_trop_variance_joint_full(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    lambda_time: f64,
    lambda_unit: f64,
    lambda_nn: f64,
    n_bootstrap: usize,
    max_iter: usize,
    tol: f64,
    seed: u64,
    alpha: f64,
) -> BootstrapResult {
    let y_arr = y.to_owned();
    let d_arr = d.to_owned();

    let n_units = y_arr.ncols();
    let n_periods = y_arr.nrows();

    let classification = classify_units(&d_arr.view());

    // Determine the number of treated periods from the original D matrix.
    // This count is fixed across all bootstrap draws.
    let mut first_treat_period = n_periods;
    for t in 0..n_periods {
        for i in 0..n_units {
            if d_arr[[t, i]] == 1.0 {
                first_treat_period = first_treat_period.min(t);
                break;
            }
        }
    }
    let treated_periods = n_periods.saturating_sub(first_treat_period);

    // Map infinite lambda values to effective computation values.
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

    // Parallel bootstrap: failed iterations yield None and are discarded.
    let bootstrap_estimates: Vec<f64> = (0..n_bootstrap)
        .into_par_iter()
        .filter_map(|b| {
            let iteration_seed = seed.wrapping_add(b as u64);
            let sampled_units = stratified_sample(&classification, iteration_seed);

            let (y_boot, d_boot) = build_bootstrap_matrices(&y_arr, &d_arr, &sampled_units);

            // Compute joint weights using the original treated-period count.
            let delta = compute_joint_weights(
                &y_boot.view(),
                &d_boot.view(),
                lt_eff,
                lu_eff,
                treated_periods,
            );

            // Solve the joint model; branch on whether low-rank is active.
            //
            // τ is post-hoc: mean residual (Y − μ − α − β − L) over treated cells.
            // When λ_nn ≥ 1e10 we skip the low-rank fit and L ≡ 0.
            let result = if ln_eff >= 1e10 {
                solve_joint_no_lowrank(&y_boot.view(), &delta.view()).map(
                    |(mu, alpha_est, beta)| {
                        let mut tau_sum = 0.0_f64;
                        let mut tau_count = 0usize;
                        let nu_b = y_boot.ncols();
                        for t in 0..n_periods {
                            for i in 0..nu_b {
                                if d_boot[[t, i]] == 1.0 && y_boot[[t, i]].is_finite() {
                                    tau_sum += y_boot[[t, i]] - mu - alpha_est[i] - beta[t];
                                    tau_count += 1;
                                }
                            }
                        }
                        if tau_count > 0 {
                            tau_sum / tau_count as f64
                        } else {
                            f64::NAN
                        }
                    },
                )
            } else {
                solve_joint_with_lowrank(
                    &y_boot.view(),
                    &d_boot.view(),
                    &delta.view(),
                    ln_eff,
                    max_iter,
                    tol,
                )
                .map(|(_, _, _, _, tau, _, _)| tau)
            };

            result.filter(|tau| tau.is_finite())
        })
        .collect();

    let n_valid = bootstrap_estimates.len();

    let (mean, _, se) = compute_bootstrap_variance(&bootstrap_estimates);
    let (ci_lower, ci_upper) = compute_percentile_ci(&bootstrap_estimates, alpha);

    BootstrapResult {
        estimates: bootstrap_estimates,
        se,
        mean,
        ci_lower,
        ci_upper,
        n_valid,
        n_total: n_bootstrap,
        level: 1.0 - alpha,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_classify_units() {
        // Create a simple D matrix: 3 periods, 4 units
        // Unit 0, 1: never treated (control)
        // Unit 2, 3: treated at some point
        let d = Array2::from_shape_vec(
            (3, 4),
            vec![
                0.0, 0.0, 0.0, 0.0, // t=0
                0.0, 0.0, 1.0, 0.0, // t=1
                0.0, 0.0, 1.0, 1.0, // t=2
            ],
        )
        .unwrap();

        let classification = classify_units(&d.view());

        assert_eq!(classification.n_control, 2);
        assert_eq!(classification.n_treated, 2);
        assert_eq!(classification.control_units, vec![0, 1]);
        assert_eq!(classification.treated_units, vec![2, 3]);
    }

    #[test]
    fn test_stratified_sample_preserves_counts() {
        let classification = UnitClassification {
            control_units: vec![0, 1, 2, 3, 4],
            treated_units: vec![5, 6, 7],
            n_control: 5,
            n_treated: 3,
        };

        let seed = 42u64;
        let sampled = stratified_sample(&classification, seed);

        // Should sample exactly n_control + n_treated units
        assert_eq!(sampled.len(), 8);
    }

    #[test]
    fn test_stratified_sample_deterministic() {
        let classification = UnitClassification {
            control_units: vec![0, 1, 2, 3, 4],
            treated_units: vec![5, 6, 7],
            n_control: 5,
            n_treated: 3,
        };

        let seed = 42u64;
        let sampled1 = stratified_sample(&classification, seed);
        let sampled2 = stratified_sample(&classification, seed);

        // Same seed should produce same result
        assert_eq!(sampled1, sampled2);
    }

    #[test]
    fn test_stratified_sample_different_seeds() {
        let classification = UnitClassification {
            control_units: vec![0, 1, 2, 3, 4],
            treated_units: vec![5, 6, 7],
            n_control: 5,
            n_treated: 3,
        };

        let sampled1 = stratified_sample(&classification, 42);
        let sampled2 = stratified_sample(&classification, 43);

        // Different seeds should produce different results
        assert_ne!(sampled1, sampled2);
    }

    #[test]
    fn test_build_bootstrap_matrices() {
        let y = Array2::from_shape_vec(
            (2, 3),
            vec![
                1.0, 2.0, 3.0, // t=0
                4.0, 5.0, 6.0, // t=1
            ],
        )
        .unwrap();

        let d = Array2::from_shape_vec(
            (2, 3),
            vec![
                0.0, 0.0, 0.0, // t=0
                0.0, 1.0, 1.0, // t=1
            ],
        )
        .unwrap();

        // Sample units [1, 0, 2] (with replacement)
        let sampled_units = vec![1, 0, 2];
        let (y_boot, d_boot) = build_bootstrap_matrices(&y, &d, &sampled_units);

        // Check dimensions
        assert_eq!(y_boot.shape(), &[2, 3]);
        assert_eq!(d_boot.shape(), &[2, 3]);

        // Check values - unit 1 should be in position 0
        assert_eq!(y_boot[[0, 0]], 2.0);
        assert_eq!(y_boot[[1, 0]], 5.0);

        // Check values - unit 0 should be in position 1
        assert_eq!(y_boot[[0, 1]], 1.0);
        assert_eq!(y_boot[[1, 1]], 4.0);
    }

    #[test]
    fn test_compute_bootstrap_variance_population() {
        // variance uses 1/B (population variance) matching standard algorithm
        let estimates = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (mean, variance, se) = compute_bootstrap_variance(&estimates);

        // Mean should be 3.0
        assert!((mean - 3.0).abs() < 1e-10);

        // Variance with 1/B: sum((x-3)^2) / 5 = 10/5 = 2.0
        assert!((variance - 2.0).abs() < 1e-10);

        // SE = sqrt(2.0)
        assert!((se - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_compute_bootstrap_variance_empty() {
        let estimates: Vec<f64> = vec![];
        let (mean, variance, se) = compute_bootstrap_variance(&estimates);

        assert_eq!(mean, 0.0);
        assert_eq!(variance, 0.0);
        assert_eq!(se, 0.0);
    }

    #[test]
    fn test_compute_bootstrap_variance_single() {
        let estimates = vec![5.0];
        let (mean, variance, se) = compute_bootstrap_variance(&estimates);

        // Single value -> mean = that value, variance/se = 0
        assert_eq!(mean, 5.0);
        assert_eq!(variance, 0.0);
        assert_eq!(se, 0.0);
    }

    #[test]
    fn test_compute_percentile_ci() {
        // Create sorted estimates: 1, 2, 3, ..., 100
        let estimates: Vec<f64> = (1..=100).map(|x| x as f64).collect();

        let (ci_lower, ci_upper) = compute_percentile_ci(&estimates, 0.05);

        // For 95% CI with 100 values:
        // Lower: index = (100-1) * 0.025 = 2.475 → interpolate between sorted[2] and sorted[3]
        // Upper: index = (100-1) * 0.975 = 96.525 → interpolate between sorted[96] and sorted[97]
        assert!(
            (ci_lower - 3.475).abs() < 1e-10,
            "ci_lower={}, expected 3.475",
            ci_lower
        );
        assert!(
            (ci_upper - 97.525).abs() < 1e-10,
            "ci_upper={}, expected 97.525",
            ci_upper
        );
    }

    #[test]
    fn test_compute_percentile_ci_empty() {
        let estimates: Vec<f64> = vec![];
        let (ci_lower, ci_upper) = compute_percentile_ci(&estimates, 0.05);

        assert!(ci_lower.is_nan());
        assert!(ci_upper.is_nan());
    }

    /// NaN values must be filtered before CI computation.
    #[test]
    fn test_compute_percentile_ci_with_nan() {
        // Mix of valid values and NaN
        let estimates = vec![1.0, 2.0, f64::NAN, 3.0, f64::NAN, 4.0, 5.0];
        let (ci_lower, ci_upper) = compute_percentile_ci(&estimates, 0.05);

        // CI should be computed from [1, 2, 3, 4, 5] only
        assert!(ci_lower.is_finite(), "ci_lower should be finite, got {}", ci_lower);
        assert!(ci_upper.is_finite(), "ci_upper should be finite, got {}", ci_upper);
        assert!(ci_lower <= ci_upper, "ci_lower ({}) should be <= ci_upper ({})", ci_lower, ci_upper);

        // Compare with NaN-free version
        let clean = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (ci_lower_clean, ci_upper_clean) = compute_percentile_ci(&clean, 0.05);
        assert!((ci_lower - ci_lower_clean).abs() < 1e-12,
            "NaN-filtered CI lower ({}) should match clean ({})", ci_lower, ci_lower_clean);
        assert!((ci_upper - ci_upper_clean).abs() < 1e-12,
            "NaN-filtered CI upper ({}) should match clean ({})", ci_upper, ci_upper_clean);
    }

    /// All-NaN input yields (NaN, NaN).
    #[test]
    fn test_compute_percentile_ci_all_nan() {
        let estimates = vec![f64::NAN, f64::NAN, f64::NAN];
        let (ci_lower, ci_upper) = compute_percentile_ci(&estimates, 0.05);
        assert!(ci_lower.is_nan());
        assert!(ci_upper.is_nan());
    }

    /// Inf values must be filtered identically to NaN.
    #[test]
    fn test_compute_percentile_ci_with_inf() {
        let estimates = vec![1.0, 2.0, f64::INFINITY, 3.0, f64::NEG_INFINITY, 4.0, 5.0];
        let (ci_lower, ci_upper) = compute_percentile_ci(&estimates, 0.05);

        let clean = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (ci_lower_clean, ci_upper_clean) = compute_percentile_ci(&clean, 0.05);
        assert!((ci_lower - ci_lower_clean).abs() < 1e-12);
        assert!((ci_upper - ci_upper_clean).abs() < 1e-12);
    }

    /// Variance computation filters NaN values transparently.
    #[test]
    fn test_compute_bootstrap_variance_with_nan() {
        let estimates = vec![1.0, 2.0, f64::NAN, 3.0, 4.0, 5.0];
        let (mean, variance, se) = compute_bootstrap_variance(&estimates);

        let clean = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (mean_clean, var_clean, se_clean) = compute_bootstrap_variance(&clean);

        assert!((mean - mean_clean).abs() < 1e-12);
        assert!((variance - var_clean).abs() < 1e-12);
        assert!((se - se_clean).abs() < 1e-12);
    }

    #[test]
    fn test_seed_increment() {
        // Test that different bootstrap iterations use different seeds
        let seed = 42u64;

        let mut rng1 = Xoshiro256PlusPlus::seed_from_u64(seed.wrapping_add(0));
        let mut rng2 = Xoshiro256PlusPlus::seed_from_u64(seed.wrapping_add(1));

        let val1: u64 = rng1.gen();
        let val2: u64 = rng2.gen();

        // Different seeds should produce different values
        assert_ne!(val1, val2);
    }

    #[test]
    fn test_population_variance_vs_bessel() {
        // Verify our function uses 1/B (population), not 1/(B-1) (Bessel)
        let estimates = [1.0, 2.0, 3.0, 4.0, 5.0];
        let n = estimates.len() as f64;
        let mean = estimates.iter().sum::<f64>() / n;

        // Population variance (1/B)
        let variance_population = estimates.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;

        // Bessel-corrected variance (1/(B-1))
        let variance_bessel = estimates.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);

        // Our function should match population variance
        let (_, variance_func, _) = compute_bootstrap_variance(&estimates.to_vec());
        assert!((variance_func - variance_population).abs() < 1e-10);

        // And should NOT match Bessel variance
        assert!((variance_func - variance_bessel).abs() > 1e-3);
    }

    // ========================================================================
    // SE precision, CI reasonableness, seed determinism
    // ========================================================================

    /// SE is non-negative for any input.
    #[test]
    fn test_se_positive_definiteness() {
        // Test with various estimate distributions
        let test_cases: Vec<Vec<f64>> = vec![
            vec![1.0, 2.0, 3.0, 4.0, 5.0],           // Normal spread
            vec![0.0, 0.0, 0.0, 0.0, 0.0],           // All zeros
            vec![1.0, 1.0, 1.0, 1.0, 1.0],           // All same
            vec![-5.0, -3.0, 0.0, 3.0, 5.0],         // Symmetric around 0
            vec![100.0, 100.1, 100.2, 100.3, 100.4], // Small variance
            vec![-1000.0, 0.0, 1000.0],              // Large variance
        ];

        for estimates in test_cases {
            let (_, variance, se) = compute_bootstrap_variance(&estimates);

            // SE must be non-negative
            assert!(
                se >= 0.0,
                "SE must be non-negative, got {} for estimates {:?}",
                se,
                estimates
            );

            // Variance must be non-negative
            assert!(
                variance >= 0.0,
                "Variance must be non-negative, got {} for estimates {:?}",
                variance,
                estimates
            );

            // SE should equal sqrt(variance)
            assert!(
                (se - variance.sqrt()).abs() < 1e-12,
                "SE should equal sqrt(variance)"
            );
        }
    }

    /// CI bounds are ordered and contain the mean for a symmetric distribution.
    #[test]
    fn test_ci_reasonableness() {
        // Test with sorted estimates for predictable CI
        let estimates: Vec<f64> = (1..=100).map(|x| x as f64).collect();
        let (mean, _, _) = compute_bootstrap_variance(&estimates);
        let (ci_lower, ci_upper) = compute_percentile_ci(&estimates, 0.05);

        // CI bounds should be ordered correctly
        assert!(
            ci_lower <= ci_upper,
            "CI lower ({}) should be <= CI upper ({})",
            ci_lower,
            ci_upper
        );

        // For symmetric distribution, mean should be within CI
        assert!(
            ci_lower <= mean && mean <= ci_upper,
            "Mean ({}) should be within CI [{}, {}]",
            mean,
            ci_lower,
            ci_upper
        );

        // CI should be within data range
        let min_val = estimates.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_val = estimates.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            ci_lower >= min_val && ci_upper <= max_val,
            "CI [{}, {}] should be within data range [{}, {}]",
            ci_lower,
            ci_upper,
            min_val,
            max_val
        );
    }

    /// CI width increases with the confidence level.
    #[test]
    fn test_ci_width_vs_alpha() {
        let estimates: Vec<f64> = (1..=100).map(|x| x as f64).collect();

        // 90% CI (alpha=0.10)
        let (ci_lower_90, ci_upper_90) = compute_percentile_ci(&estimates, 0.10);
        let width_90 = ci_upper_90 - ci_lower_90;

        // 95% CI (alpha=0.05)
        let (ci_lower_95, ci_upper_95) = compute_percentile_ci(&estimates, 0.05);
        let width_95 = ci_upper_95 - ci_lower_95;

        // 99% CI (alpha=0.01)
        let (ci_lower_99, ci_upper_99) = compute_percentile_ci(&estimates, 0.01);
        let width_99 = ci_upper_99 - ci_lower_99;

        // Wider confidence level should give wider CI
        assert!(
            width_90 < width_95,
            "90% CI width ({}) should be < 95% CI width ({})",
            width_90,
            width_95
        );
        assert!(
            width_95 < width_99,
            "95% CI width ({}) should be < 99% CI width ({})",
            width_95,
            width_99
        );
    }

    /// Identical seeds produce identical sampling sequences.
    #[test]
    fn test_bootstrap_seed_determinism() {
        let classification = UnitClassification {
            control_units: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            treated_units: vec![10, 11, 12, 13, 14],
            n_control: 10,
            n_treated: 5,
        };

        let base_seed = 12345u64;
        let n_iterations = 50;

        // Run two complete bootstrap sampling sequences with same seed
        let mut samples1: Vec<Vec<usize>> = Vec::new();
        let mut samples2: Vec<Vec<usize>> = Vec::new();

        for b in 0..n_iterations {
            let iter_seed = base_seed.wrapping_add(b as u64);
            samples1.push(stratified_sample(&classification, iter_seed));
            samples2.push(stratified_sample(&classification, iter_seed));
        }

        // All samples should be identical
        for b in 0..n_iterations {
            assert_eq!(
                samples1[b], samples2[b],
                "Iteration {} should produce identical samples with same seed",
                b
            );
        }
    }

    /// Distinct base seeds produce distinct sampling sequences.
    #[test]
    fn test_bootstrap_different_seeds_differ() {
        let classification = UnitClassification {
            control_units: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            treated_units: vec![10, 11, 12, 13, 14],
            n_control: 10,
            n_treated: 5,
        };

        let n_iterations = 20;

        // Run with two different base seeds
        let mut samples_seed1: Vec<Vec<usize>> = Vec::new();
        let mut samples_seed2: Vec<Vec<usize>> = Vec::new();

        for b in 0..n_iterations {
            samples_seed1.push(stratified_sample(
                &classification,
                100u64.wrapping_add(b as u64),
            ));
            samples_seed2.push(stratified_sample(
                &classification,
                200u64.wrapping_add(b as u64),
            ));
        }

        // At least some samples should differ
        let mut any_different = false;
        for b in 0..n_iterations {
            if samples_seed1[b] != samples_seed2[b] {
                any_different = true;
                break;
            }
        }
        assert!(
            any_different,
            "Different base seeds should produce different sample sequences"
        );
    }

    /// SE converges to the theoretical value for a known distribution.
    #[test]
    fn test_bootstrap_se_precision() {
        // Uniform[0, 12]: theoretical variance = (12-0)^2 / 12 = 12.
        let n_samples = 10000;
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let estimates: Vec<f64> = (0..n_samples).map(|_| rng.gen_range(0.0..12.0)).collect();

        let (_, variance, se) = compute_bootstrap_variance(&estimates);

        // Theoretical variance for uniform [0, 12] is 12
        let theoretical_variance: f64 = 12.0;
        let theoretical_se = theoretical_variance.sqrt();

        // With 10000 samples, sample variance should be close to theoretical
        // Allow 5% relative error for statistical variation
        let relative_error = (variance - theoretical_variance).abs() / theoretical_variance;
        assert!(
            relative_error < 0.05,
            "Variance {} should be close to theoretical {} (relative error: {})",
            variance,
            theoretical_variance,
            relative_error
        );

        // SE precision check
        let se_error = (se - theoretical_se).abs();
        assert!(
            se_error < 0.1, // Allow 0.1 absolute error for SE ≈ 3.46
            "SE {} should be close to theoretical {} (error: {})",
            se,
            theoretical_se,
            se_error
        );
    }

    /// Variance formula: V = Σ(τ_b − τ̄)² / B.
    #[test]
    fn test_bootstrap_variance_formula() {
        let estimates = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        let n = estimates.len() as f64;

        // Manual calculation using population variance (1/B).
        let mean_manual = estimates.iter().sum::<f64>() / n;
        let variance_manual = estimates
            .iter()
            .map(|x| (x - mean_manual).powi(2))
            .sum::<f64>()
            / n;
        let se_manual = variance_manual.sqrt();

        // Function calculation
        let (mean_func, variance_func, se_func) = compute_bootstrap_variance(&estimates);

        // Verify with high precision (< 1e-12)
        assert!(
            (mean_func - mean_manual).abs() < 1e-12,
            "Mean mismatch: {} vs {}",
            mean_func,
            mean_manual
        );
        assert!(
            (variance_func - variance_manual).abs() < 1e-12,
            "Variance mismatch: {} vs {}",
            variance_func,
            variance_manual
        );
        assert!(
            (se_func - se_manual).abs() < 1e-12,
            "SE mismatch: {} vs {}",
            se_func,
            se_manual
        );

        // Verify expected values
        // Mean = (2+4+6+8+10)/5 = 6
        assert!((mean_func - 6.0).abs() < 1e-12);
        // Variance = ((2-6)^2 + (4-6)^2 + (6-6)^2 + (8-6)^2 + (10-6)^2) / 5
        //          = (16 + 4 + 0 + 4 + 16) / 5 = 40/5 = 8
        assert!((variance_func - 8.0).abs() < 1e-12);
        // SE = sqrt(8) ≈ 2.828
        assert!((se_func - 8.0_f64.sqrt()).abs() < 1e-12);
    }

    /// Percentile CI with linear interpolation on 10 values.
    #[test]
    fn test_percentile_ci_interpolation() {
        let estimates: Vec<f64> = (1..=10).map(|x| x as f64).collect();

        // 95% CI (alpha = 0.05):
        // lower_idx = (10-1) * 0.025 = 0.225 → interpolate between index 0 and 1
        // upper_idx = (10-1) * 0.975 = 8.775 → interpolate between index 8 and 9
        let (ci_lower, ci_upper) = compute_percentile_ci(&estimates, 0.05);

        // Lower: 1 * (1 - 0.225) + 2 * 0.225 = 1.225
        let expected_lower = 1.0 * (1.0 - 0.225) + 2.0 * 0.225;
        assert!(
            (ci_lower - expected_lower).abs() < 1e-10,
            "CI lower {} should be {}",
            ci_lower,
            expected_lower
        );

        // Upper: 9 * (1 - 0.775) + 10 * 0.775 = 9.775
        let expected_upper = 9.0 * (1.0 - 0.775) + 10.0 * 0.775;
        assert!(
            (ci_upper - expected_upper).abs() < 1e-10,
            "CI upper {} should be {}",
            ci_upper,
            expected_upper
        );
    }

    // ========================================================================
    // Regression tests: non-finite value filtering
    // ========================================================================

    #[test]
    fn test_compute_bootstrap_variance_filters_nan() {
        // NaN values should be filtered before computing statistics
        let estimates = vec![1.0, 2.0, f64::NAN, 4.0, 5.0];
        let (mean, variance, se) = compute_bootstrap_variance(&estimates);

        // Should compute stats from [1.0, 2.0, 4.0, 5.0] only
        assert!(mean.is_finite(), "Mean should be finite after NaN filtering");
        assert!((mean - 3.0).abs() < 1e-10, "Mean of [1,2,4,5] = 3.0, got {}", mean);
        assert!(variance.is_finite(), "Variance should be finite");
        assert!(se.is_finite(), "SE should be finite");
    }

    #[test]
    fn test_compute_bootstrap_variance_filters_inf() {
        // Inf values should be filtered before computing statistics
        let estimates = vec![1.0, 2.0, f64::INFINITY, 4.0, f64::NEG_INFINITY];
        let (mean, variance, se) = compute_bootstrap_variance(&estimates);

        // Should compute stats from [1.0, 2.0, 4.0] only
        assert!(mean.is_finite(), "Mean should be finite after Inf filtering");
        let expected_mean = (1.0 + 2.0 + 4.0) / 3.0;
        assert!((mean - expected_mean).abs() < 1e-10);
        assert!(variance.is_finite());
        assert!(se.is_finite());
    }

    #[test]
    fn test_compute_bootstrap_variance_all_nan() {
        // All NaN → treated as empty
        let estimates = vec![f64::NAN, f64::NAN, f64::NAN];
        let (mean, variance, se) = compute_bootstrap_variance(&estimates);
        assert_eq!(mean, 0.0);
        assert_eq!(variance, 0.0);
        assert_eq!(se, 0.0);
    }

    #[test]
    fn test_compute_percentile_ci_filters_nan() {
        // NaN values should be filtered before computing percentiles
        let estimates = vec![1.0, 2.0, f64::NAN, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let (ci_lower, ci_upper) = compute_percentile_ci(&estimates, 0.05);

        assert!(ci_lower.is_finite(), "CI lower should be finite after NaN filtering");
        assert!(ci_upper.is_finite(), "CI upper should be finite after NaN filtering");
        assert!(ci_lower < ci_upper, "CI lower < upper");
    }

    #[test]
    fn test_compute_percentile_ci_all_nan_two_elements() {
        let estimates = vec![f64::NAN, f64::NAN];
        let (ci_lower, ci_upper) = compute_percentile_ci(&estimates, 0.05);
        assert!(ci_lower.is_nan());
        assert!(ci_upper.is_nan());
    }
}
