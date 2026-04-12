//! Exponential-decay weight computation for the TROP estimator.
//!
//! Implements the distance-based weighting scheme defined in Equation (3):
//!
//!   θ_s^{i,t}(λ) = exp(−λ_time · dist_time(s, t))
//!   ω_j^{i,t}(λ) = exp(−λ_unit · dist_unit_{−t}(j, i))
//!
//! where dist_time(s, t) = |t − s| and dist_unit_{−t}(j, i) is the RMSE of
//! outcome differences over mutually untreated periods excluding t.
//!
//! Two modes are provided:
//! - Per-observation weights (Algorithm 2, Twostep): each treated (i, t) pair
//!   receives its own weight matrix W with W[s, j] = θ_s^{i,t} · ω_j^{i,t}.
//! - Global weights (Joint): a single weight matrix δ shared across all treated
//!   observations, with time center at the midpoint of the treated block and
//!   unit distance measured against the average treated trajectory.

use crate::distance::compute_unit_distance_for_obs;
use ndarray::{Array1, Array2, ArrayView2};

/// Compute the per-observation weight matrix for the Twostep method.
///
/// For a treated observation (i, t), the weight matrix is the outer product
/// of time weights and unit weights (Equation 3):
///
///   W[s, j] = θ_s · ω_j
///
/// where
///   θ_s = exp(−λ_time · |t − s|)
///   ω_j = exp(−λ_unit · dist_{−t}(j, i))
///
/// Units treated at the target period receive zero weight. The target unit
/// itself receives unit weight 1.0. Weights are unnormalized.
///
/// # Arguments
/// * `y` - Outcome matrix (T × N), column-major
/// * `d` - Treatment indicator matrix (T × N)
/// * `n_periods` - Number of time periods T
/// * `n_units` - Number of units N
/// * `target_unit` - Index of the treated unit i
/// * `target_period` - Index of the treated period t
/// * `lambda_time` - Decay parameter for time weights
/// * `lambda_unit` - Decay parameter for unit weights
/// * `time_dist` - Pre-computed time distance matrix (T × T)
///
/// # Returns
/// Weight matrix W of dimension (T × N)
#[allow(clippy::too_many_arguments)]
pub fn compute_weight_matrix(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    n_periods: usize,
    n_units: usize,
    target_unit: usize,
    target_period: usize,
    lambda_time: f64,
    lambda_unit: f64,
    time_dist: &ArrayView2<i64>,
) -> Array2<f64> {
    // Decompose into time and unit weight vectors, then form the outer product.
    let (time_weights, unit_weights) = compute_twostep_weight_vectors(
        y, d, n_periods, n_units, target_unit, target_period,
        lambda_time, lambda_unit, time_dist,
    );

    // W[s, j] = θ_s · ω_j
    let mut weight_matrix = Array2::<f64>::zeros((n_periods, n_units));
    for t in 0..n_periods {
        for i in 0..n_units {
            weight_matrix[[t, i]] = time_weights[t] * unit_weights[i];
        }
    }

    weight_matrix
}

/// Compute the global weight matrix for the Joint method.
///
/// Unlike the Twostep method, the Joint method uses a single weight matrix
/// shared across all treated observations. The weight matrix is the outer
/// product δ[s, j] = δ_time[s] · δ_unit[j], where:
///
///   δ_time[s] = exp(−λ_time · |s − center|)
///     with center = T − T_post / 2 (midpoint of the treated block)
///
///   δ_unit[j] = exp(−λ_unit · RMSE_j)
///     where RMSE_j is the root-mean-square deviation of unit j's
///     pre-treatment outcomes from the average treated trajectory.
///
/// Weights are unnormalized.
///
/// # Arguments
/// * `y` - Outcome matrix (T × N)
/// * `d` - Treatment indicator matrix (T × N)
/// * `lambda_time` - Decay parameter for time weights
/// * `lambda_unit` - Decay parameter for unit weights
/// * `treated_periods` - Number of post-treatment periods T_post
///
/// # Returns
/// Weight matrix δ of dimension (T × N)
pub fn compute_joint_weights(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    lambda_time: f64,
    lambda_unit: f64,
    treated_periods: usize,
) -> Array2<f64> {
    let n_periods = y.nrows();
    let n_units = y.ncols();

    let (delta_time, delta_unit) = compute_joint_weight_vectors(
        y, d, lambda_time, lambda_unit, treated_periods,
    );

    // δ[s, j] = δ_time[s] · δ_unit[j]
    let mut delta = Array2::<f64>::zeros((n_periods, n_units));
    for t in 0..n_periods {
        for i in 0..n_units {
            delta[[t, i]] = delta_time[t] * delta_unit[i];
        }
    }

    delta
}

/// Decompose Twostep weights into separate time and unit vectors.
///
/// Returns (θ, ω) where:
///   θ_s = exp(−λ_time · |t − s|)           (Equation 3, time component)
///   ω_j = exp(−λ_unit · dist_{−t}(j, i))   (Equation 3, unit component)
///
/// The target unit receives ω_i = 1.0. Units treated at the target period
/// and units with non-finite distance receive ω_j = 0.0. When λ_unit = 0,
/// all eligible control units receive ω_j = 1.0 (uniform weighting).
///
/// # Arguments
/// Same as [`compute_weight_matrix`].
///
/// # Returns
/// `(time_weights, unit_weights)` — both `Array1<f64>` of length T and N
#[allow(clippy::too_many_arguments)]
pub fn compute_twostep_weight_vectors(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    n_periods: usize,
    n_units: usize,
    target_unit: usize,
    target_period: usize,
    lambda_time: f64,
    lambda_unit: f64,
    time_dist: &ArrayView2<i64>,
) -> (Array1<f64>, Array1<f64>) {
    // θ_s = exp(−λ_time · |t − s|)
    let time_weights: Array1<f64> = Array1::from_shape_fn(n_periods, |s| {
        let dist = time_dist[[target_period, s]] as f64;
        (-lambda_time * dist).exp()
    });

    // ω_j = exp(−λ_unit · dist_{−t}(j, i))
    let mut unit_weights = Array1::<f64>::zeros(n_units);

    if lambda_unit == 0.0 {
        // Uniform weighting: all untreated units at the target period get ω_j = 1.
        for j in 0..n_units {
            if d[[target_period, j]] == 0.0 && j != target_unit {
                unit_weights[j] = 1.0;
            }
        }
    } else {
        for j in 0..n_units {
            if d[[target_period, j]] == 0.0 && j != target_unit {
                let dist = compute_unit_distance_for_obs(y, d, j, target_unit, target_period);
                if dist.is_finite() {
                    unit_weights[j] = (-lambda_unit * dist).exp();
                }
            }
        }
    }

    // The target unit always participates with weight 1.
    unit_weights[target_unit] = 1.0;

    (time_weights, unit_weights)
}

/// Decompose Joint weights into separate time and unit vectors.
///
/// Returns (δ_time, δ_unit) where:
///   δ_time[s] = exp(−λ_time · |s − center|)
///     with center = T − T_post / 2
///   δ_unit[j] = exp(−λ_unit · RMSE_j)
///     where RMSE_j is computed over pre-treatment periods against the
///     average treated trajectory.
///
/// Units with no valid pre-treatment observations receive δ_unit = 0.
/// When there are no pre-treatment periods, all units receive δ_unit = 1.
///
/// # Arguments
/// Same as [`compute_joint_weights`].
///
/// # Returns
/// `(delta_time, delta_unit)` — both `Array1<f64>` of length T and N
pub fn compute_joint_weight_vectors(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    lambda_time: f64,
    lambda_unit: f64,
    treated_periods: usize,
) -> (Array1<f64>, Array1<f64>) {
    let n_periods = y.nrows();
    let n_units = y.ncols();

    // Identify ever-treated units.
    let mut treated_unit_idx: Vec<usize> = Vec::new();
    for i in 0..n_units {
        if (0..n_periods).any(|t| d[[t, i]] == 1.0) {
            treated_unit_idx.push(i);
        }
    }

    // δ_time[s] = exp(−λ_time · |s − center|), center = T − T_post / 2.
    let center = n_periods as f64 - treated_periods as f64 / 2.0;
    let mut delta_time = Array1::<f64>::zeros(n_periods);
    for t in 0..n_periods {
        let dist = (t as f64 - center).abs();
        delta_time[t] = (-lambda_time * dist).exp();
    }

    // δ_unit[j] = exp(−λ_unit · RMSE_j), where RMSE_j is the root-mean-square
    // deviation of unit j from the average treated trajectory over pre-periods.
    let n_pre = n_periods.saturating_sub(treated_periods);

    // Average outcome trajectory across treated units.
    let mut average_treated = Array1::<f64>::from_elem(n_periods, f64::NAN);
    if !treated_unit_idx.is_empty() {
        for t in 0..n_periods {
            let mut sum = 0.0;
            let mut count = 0;
            for &i in &treated_unit_idx {
                if y[[t, i]].is_finite() {
                    sum += y[[t, i]];
                    count += 1;
                }
            }
            if count > 0 {
                average_treated[t] = sum / count as f64;
            }
        }
    }

    let mut delta_unit = Array1::<f64>::zeros(n_units);
    for i in 0..n_units {
        if n_pre > 0 {
            let mut sum_sq = 0.0;
            let mut n_valid = 0;
            for t in 0..n_pre {
                if y[[t, i]].is_finite() && average_treated[t].is_finite() {
                    let diff = average_treated[t] - y[[t, i]];
                    sum_sq += diff * diff;
                    n_valid += 1;
                }
            }
            let dist = if n_valid > 0 {
                (sum_sq / n_valid as f64).sqrt()
            } else {
                f64::INFINITY
            };
            // Guard against NaN from IEEE 754: −0.0 × ∞ = NaN, exp(NaN) = NaN.
            // Units with no valid pre-period data receive zero weight.
            delta_unit[i] = if dist.is_infinite() {
                0.0
            } else {
                (-lambda_unit * dist).exp()
            };
        } else {
            // No pre-treatment periods; assign uniform weight.
            delta_unit[i] = 1.0;
        }
    }

    (delta_time, delta_unit)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_weight_matrix_structure() {
        // 3 periods, 2 units, all control.
        let y = array![[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]];
        let d = array![[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]];
        let time_dist = array![[0i64, 1, 2], [1, 0, 1], [2, 1, 0]];

        let weights = compute_weight_matrix(
            &y.view(),
            &d.view(),
            3,
            2,
            0,
            1,
            0.5,
            0.5,
            &time_dist.view(),
        );

        let total: f64 = weights.sum();
        assert!(total > 0.0, "Weights should be positive, got {}", total);
    }

    #[test]
    fn test_weight_matrix_zero_lambda() {
        // λ_time = λ_unit = 0 ⟹ all weights equal exp(0) = 1 (unnormalized).
        let y = array![[1.0, 2.0, 3.0], [2.0, 3.0, 4.0], [3.0, 4.0, 5.0]];
        let d = array![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let time_dist = array![[0i64, 1, 2], [1, 0, 1], [2, 1, 0]];

        let weights = compute_weight_matrix(
            &y.view(),
            &d.view(),
            3,
            3,
            0,
            1,
            0.0,
            0.0,
            &time_dist.view(),
        );

        let expected = 1.0;
        for t in 0..3 {
            for i in 0..3 {
                assert!(
                    (weights[[t, i]] - expected).abs() < 1e-10,
                    "Expected uniform weight {}, got {} at [{}, {}]",
                    expected,
                    weights[[t, i]],
                    t,
                    i
                );
            }
        }

        let total: f64 = weights.sum();
        assert!(
            (total - 9.0).abs() < 1e-10,
            "Weights should sum to 9, got {}",
            total
        );
    }

    #[test]
    fn test_joint_weights_unnormalized() {
        // Joint weights are unnormalized; verify they do not sum to 1.
        let y = array![[1.0, 2.0], [2.0, 3.0], [3.0, 4.0], [4.0, 5.0]];
        let d = array![
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 1.0],
            [0.0, 1.0]
        ];

        let weights = compute_joint_weights(&y.view(), &d.view(), 0.5, 0.5, 2);

        let total: f64 = weights.sum();
        assert!(total > 0.0, "Weights should be positive");
    }

    #[test]
    fn test_joint_weights_time_center() {
        // λ_time = 0 ⟹ uniform time weights; all rows should have equal sums.
        let y = array![[1.0, 2.0], [2.0, 3.0], [3.0, 4.0], [4.0, 5.0]];
        let d = array![[0.0, 0.0], [0.0, 0.0], [0.0, 1.0], [0.0, 1.0]];

        let weights = compute_joint_weights(&y.view(), &d.view(), 0.0, 0.0, 2);

        let row0_sum: f64 = weights.row(0).sum();
        let row1_sum: f64 = weights.row(1).sum();
        let row2_sum: f64 = weights.row(2).sum();
        let row3_sum: f64 = weights.row(3).sum();

        assert!((row0_sum - row1_sum).abs() < 1e-10);
        assert!((row1_sum - row2_sum).abs() < 1e-10);
        assert!((row2_sum - row3_sum).abs() < 1e-10);
    }

    // ========================================================================
    // Joint weight property tests
    // ========================================================================

    #[test]
    fn test_joint_weights_outer_product_structure() {
        // δ[t,i] = δ_time[t] · δ_unit[i] ⟹ δ[t,i]/δ[t,j] = δ[s,i]/δ[s,j].
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 6.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0]
        ];

        let delta = compute_joint_weights(&y.view(), &d.view(), 0.5, 0.5, 2);

        // Verify rank-1 structure: δ[t,i]/δ[t,j] = δ[s,i]/δ[s,j]
        // for all (t, s, i, j) where denominators are non-zero.
        let (n_periods, n_units) = delta.dim();

        for t in 0..n_periods {
            for s in 0..n_periods {
                for i in 0..n_units {
                    for j in 0..n_units {
                        if delta[[t, j]] > 1e-10 && delta[[s, j]] > 1e-10 {
                            let ratio_t = delta[[t, i]] / delta[[t, j]];
                            let ratio_s = delta[[s, i]] / delta[[s, j]];
                            assert!(
                                (ratio_t - ratio_s).abs() < 1e-8,
                                "Outer product structure violated at t={}, s={}, i={}, j={}",
                                t,
                                s,
                                i,
                                j
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_joint_weights_time_formula() {
        // δ_time[t] = exp(−λ_time · |t − center|), center = T − T_post / 2.
        let y = array![
            [1.0, 2.0],
            [2.0, 3.0],
            [3.0, 4.0],
            [4.0, 5.0],
            [5.0, 6.0],
            [6.0, 7.0]
        ];
        let d = array![
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 1.0],
            [0.0, 1.0]
        ];

        let lambda_time = 0.3;
        let treated_periods = 2;
        let n_periods = 6;
        let center = n_periods as f64 - treated_periods as f64 / 2.0;

        let delta = compute_joint_weights(&y.view(), &d.view(), lambda_time, 0.0, treated_periods);

        // With λ_unit = 0, δ_unit[j] = 1 for all j, so δ[t, j] = δ_time[t].
        for t in 0..n_periods {
            let expected_time_weight = (-lambda_time * (t as f64 - center).abs()).exp();
            let actual = delta[[t, 0]];
            assert!(
                (actual - expected_time_weight).abs() < 1e-8,
                "Time weight formula incorrect at t={}: expected {}, got {}",
                t,
                expected_time_weight,
                actual
            );
        }
    }

    #[test]
    fn test_joint_weights_unit_formula() {
        // δ_unit[j] = exp(−λ_unit · RMSE_j) over pre-treatment periods.
        // Unit 2 has a trajectory far from the treated average ⟹ lower weight.
        let y = array![
            [1.0, 1.0, 10.0],
            [2.0, 2.0, 20.0],
            [3.0, 3.0, 30.0],
            [4.0, 5.0, 40.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0]
        ];

        let lambda_unit = 0.1;
        let delta = compute_joint_weights(&y.view(), &d.view(), 0.0, lambda_unit, 1);

        assert!(
            delta[[0, 2]] < delta[[0, 0]],
            "Unit far from treated average should have lower weight"
        );
    }

    #[test]
    fn test_joint_weights_non_negative() {
        // exp(·) ≥ 0 for all real arguments; verify across λ combinations.
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 6.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0]
        ];

        for lambda_time in [0.0, 0.5, 1.0, 2.0] {
            for lambda_unit in [0.0, 0.5, 1.0, 2.0] {
                let delta =
                    compute_joint_weights(&y.view(), &d.view(), lambda_time, lambda_unit, 2);

                for t in 0..4 {
                    for i in 0..3 {
                        assert!(
                            delta[[t, i]] >= 0.0,
                            "Weight should be non-negative at [{}, {}] with λ_t={}, λ_u={}",
                            t,
                            i,
                            lambda_time,
                            lambda_unit
                        );
                    }
                }
            }
        }
    }

    // ========================================================================
    // Twostep weight tests
    // ========================================================================

    #[test]
    fn test_twostep_weight_non_negative() {
        // exp(·) ≥ 0 ⟹ all Twostep weights are non-negative.
        let y = array![
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 3.0, 4.0, 5.0],
            [3.0, 4.0, 5.0, 6.0],
            [4.0, 5.0, 6.0, 7.0],
            [5.0, 6.0, 7.0, 8.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 1.0]
        ];
        let time_dist = array![
            [0i64, 1, 2, 3, 4],
            [1, 0, 1, 2, 3],
            [2, 1, 0, 1, 2],
            [3, 2, 1, 0, 1],
            [4, 3, 2, 1, 0]
        ];

        // Sweep over λ combinations.
        for lambda_time in [0.0, 0.5, 1.0, 2.0] {
            for lambda_unit in [0.0, 0.5, 1.0, 2.0] {
                let weights = compute_weight_matrix(
                    &y.view(),
                    &d.view(),
                    5,
                    4,
                    3,
                    4,
                    lambda_time,
                    lambda_unit,
                    &time_dist.view(),
                );

                for t in 0..5 {
                    for i in 0..4 {
                        assert!(
                            weights[[t, i]] >= 0.0,
                            "Twostep weight should be non-negative at [{}, {}] with λ_t={}, λ_u={}",
                            t,
                            i,
                            lambda_time,
                            lambda_unit
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_twostep_weight_positive() {
        // Total weight sum should be strictly positive for any λ ≥ 0.
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 6.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0]
        ];
        let time_dist = array![[0i64, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]];

        for lambda_time in [0.0, 0.5, 1.0, 2.0] {
            for lambda_unit in [0.0, 0.5, 1.0, 2.0] {
                let weights = compute_weight_matrix(
                    &y.view(),
                    &d.view(),
                    4,
                    3,
                    2,
                    3,
                    lambda_time,
                    lambda_unit,
                    &time_dist.view(),
                );

                let total: f64 = weights.sum();
                assert!(
                    total > 0.0,
                    "Twostep weights should be positive, got {} with λ_t={}, λ_u={}",
                    total,
                    lambda_time,
                    lambda_unit
                );
            }
        }
    }

    #[test]
    fn test_twostep_time_weight_formula() {
        // θ_s = exp(−λ_time · |t − s|), unnormalized.
        let y = array![[1.0, 2.0], [2.0, 3.0], [3.0, 4.0], [4.0, 5.0], [5.0, 6.0]];
        let d = array![[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 1.0], [0.0, 1.0]];
        let time_dist = array![
            [0i64, 1, 2, 3, 4],
            [1, 0, 1, 2, 3],
            [2, 1, 0, 1, 2],
            [3, 2, 1, 0, 1],
            [4, 3, 2, 1, 0]
        ];

        let lambda_time = 0.5;
        let target_period = 3;

        // Expected θ_s = exp(−0.5 · |3 − s|).
        let mut expected_time_weights = [0.0f64; 5];
        for (s, weight) in expected_time_weights.iter_mut().enumerate() {
            let dist = (target_period as i64 - s as i64).abs() as f64;
            *weight = (-lambda_time * dist).exp();
        }

        // With λ_unit = 0, ω_j = 1 for all j, so row sum = θ_s · N.
        let weights = compute_weight_matrix(
            &y.view(),
            &d.view(),
            5,
            2,
            1,
            target_period,
            lambda_time,
            0.0,
            &time_dist.view(),
        );

        for (t, &expected_weight) in expected_time_weights.iter().enumerate() {
            let row_sum: f64 = weights.row(t).sum();
            let expected_row_sum = expected_weight * 2.0;
            assert!(
                (row_sum - expected_row_sum).abs() < 1e-10,
                "Time weight mismatch at t={}: expected {}, got {}",
                t,
                expected_row_sum,
                row_sum
            );
        }
    }

    #[test]
    fn test_twostep_lambda_decay_behavior() {
        // Larger λ_time ⟹ faster exponential decay ⟹ higher concentration
        // of weight near the target period.
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 6.0],
            [5.0, 6.0, 7.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0]
        ];
        let time_dist = array![
            [0i64, 1, 2, 3, 4],
            [1, 0, 1, 2, 3],
            [2, 1, 0, 1, 2],
            [3, 2, 1, 0, 1],
            [4, 3, 2, 1, 0]
        ];

        let target_period = 4;

        // λ_time = 0.1 (slow decay)
        let weights_low = compute_weight_matrix(
            &y.view(),
            &d.view(),
            5,
            3,
            2,
            target_period,
            0.1,
            0.0,
            &time_dist.view(),
        );

        // λ_time = 2.0 (fast decay)
        let weights_high = compute_weight_matrix(
            &y.view(),
            &d.view(),
            5,
            3,
            2,
            target_period,
            2.0,
            0.0,
            &time_dist.view(),
        );

        // The ratio W[target_period, ·] / W[distant_period, ·] should grow with λ.
        let weight_at_target_low: f64 = weights_low.row(target_period).sum();
        let weight_at_target_high: f64 = weights_high.row(target_period).sum();

        let distant_period = 0;
        let weight_at_distant_low: f64 = weights_low.row(distant_period).sum();
        let weight_at_distant_high: f64 = weights_high.row(distant_period).sum();

        let ratio_low = weight_at_target_low / weight_at_distant_low;
        let ratio_high = weight_at_target_high / weight_at_distant_high;

        assert!(
            ratio_high > ratio_low,
            "Larger λ should increase target/distant ratio: low={}, high={}",
            ratio_low,
            ratio_high
        );

        assert!(
            weight_at_distant_high < weight_at_distant_low,
            "Larger λ should reduce weight at distant periods: low={}, high={}",
            weight_at_distant_low,
            weight_at_distant_high
        );
    }

    #[test]
    fn test_twostep_weight_outer_product_structure() {
        // W[t,i] = θ_t · ω_i ⟹ W[t,i]/W[t,j] = W[s,i]/W[s,j].
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 6.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0]
        ];
        let time_dist = array![[0i64, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]];

        let weights = compute_weight_matrix(
            &y.view(),
            &d.view(),
            4,
            3,
            2,
            3,
            0.5,
            0.5,
            &time_dist.view(),
        );

        // Verify outer product: W[t,i] / W[t,j] = W[s,i] / W[s,j]
        for t in 0..4 {
            for s in 0..4 {
                for i in 0..3 {
                    for j in 0..3 {
                        if weights[[t, j]] > 1e-10 && weights[[s, j]] > 1e-10 {
                            let ratio_t = weights[[t, i]] / weights[[t, j]];
                            let ratio_s = weights[[s, i]] / weights[[s, j]];
                            assert!(
                                (ratio_t - ratio_s).abs() < 1e-10,
                                "Outer product structure violated at t={}, s={}, i={}, j={}",
                                t,
                                s,
                                i,
                                j
                            );
                        }
                    }
                }
            }
        }
    }

    /// When a unit has no valid pre-period data, dist = ∞. IEEE 754 gives
    /// −0.0 × ∞ = NaN, so exp(NaN) = NaN. The implementation guards against
    /// this by returning δ_unit = 0 for such units.
    #[test]
    fn test_joint_weights_infinite_distance_guard() {
        // Unit 2: all NaN in pre-periods ⟹ dist = ∞ ⟹ δ_unit = 0.
        let y = array![
            [1.0, 2.0, f64::NAN],
            [2.0, 3.0, f64::NAN],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 6.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0]
        ];

        let weights = compute_joint_weights(&y.view(), &d.view(), 0.5, 0.0, 2);

        for t in 0..4 {
            for i in 0..3 {
                assert!(
                    !weights[[t, i]].is_nan(),
                    "NaN weight at [{}, {}] with lambda_unit=0 and NaN unit",
                    t, i
                );
                assert!(
                    weights[[t, i]] >= 0.0,
                    "Weight should be non-negative at [{}, {}]", t, i
                );
            }
        }

        // Unit 2 should have zero weight (no valid pre-period data)
        for t in 0..4 {
            assert!(
                weights[[t, 2]] == 0.0,
                "Unit with no valid pre-period data should have zero weight, got {} at t={}",
                weights[[t, 2]], t
            );
        }

        // Units 0 and 1 should have positive weights (valid pre-period data)
        assert!(weights[[0, 0]] > 0.0, "Unit 0 should have positive weight");
        assert!(weights[[0, 1]] > 0.0, "Unit 1 should have positive weight");
    }

    /// Same guard tested at the vector level: δ_unit[j] = 0 when dist = ∞.
    #[test]
    fn test_joint_weight_vectors_infinite_distance_guard() {
        let y = array![
            [1.0, 2.0, f64::NAN],
            [2.0, 3.0, f64::NAN],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 6.0]
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0]
        ];

        let (delta_time, delta_unit) =
            compute_joint_weight_vectors(&y.view(), &d.view(), 0.5, 0.0, 2);

        // delta_unit[2] should be 0.0, not NaN
        assert!(
            !delta_unit[2].is_nan(),
            "delta_unit NaN for unit with no valid pre-period data"
        );
        assert_eq!(
            delta_unit[2], 0.0,
            "Unit with no valid pre-period data should have zero weight"
        );

        // Other units should have positive weights
        assert!(delta_unit[0] > 0.0);
        assert!(delta_unit[1] > 0.0);

        // Time weights should all be finite and positive
        for t in 0..4 {
            assert!(delta_time[t].is_finite() && delta_time[t] > 0.0);
        }
    }

    /// Numerical cross-validation: compare Twostep weight matrix against
    /// independently computed reference values (T=5, N=4, seed=42,
    /// target_unit=2, target_period=3, λ_time=0.5, λ_unit=0.5).
    #[test]
    fn test_weight_matrix_reference_values() {
        let y = array![
            [ 0.496714153011233,  0.361735698828815,  1.647688538100692,  3.023029856408026],
            [-0.234153374723336,  0.265863043050819,  2.579212815507391,  2.267434729152909],
            [-0.469474385934952,  1.042560043585965,  0.536582307187538,  1.034270246429743],
            [ 0.241962271566034, -1.413280244657798, -0.724917832513033,  0.937712470759027],
            [-1.012831120334424,  0.814247332595274,  0.091975924478789,  0.087696298664708]
        ];
        let d = array![
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0, 0.0]
        ];
        let time_dist = array![
            [0i64, 1, 2, 3, 4],
            [1, 0, 1, 2, 3],
            [2, 1, 0, 1, 2],
            [3, 2, 1, 0, 1],
            [4, 3, 2, 1, 0]
        ];

        let weights = compute_weight_matrix(
            &y.view(), &d.view(), 5, 4, 2, 3, 0.5, 0.5, &time_dist.view(),
        );

        // Reference values (tolerance < 1e-10).
        let reference = [
            [0.088540252840800, 0.102500687323097, 0.223130160148430, 0.144900482487827],
            [0.145978198171795, 0.168995063450973, 0.367879441171442, 0.238900507612393],
            [0.240677360384317, 0.278625755754937, 0.606530659712633, 0.393880348481609],
            [0.396809883441583, 0.459376210078063, 1.000000000000000, 0.649398908652408],
            [0.240677360384317, 0.278625755754937, 0.606530659712633, 0.393880348481609],
        ];

        for t in 0..5 {
            for i in 0..4 {
                assert!(
                    (weights[[t, i]] - reference[t][i]).abs() < 1e-10,
                    "Weight mismatch at [{}, {}]: got={:.15} expected={:.15}",
                    t, i, weights[[t, i]], reference[t][i]
                );
            }
        }
    }

    #[test]
    fn test_weight_numerical_precision() {
        // Verify individual weight entries against closed-form values
        // with tolerance < 1e-10.
        let y = array![[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]];
        let d = array![[0.0, 0.0], [0.0, 0.0], [0.0, 1.0]];
        let time_dist = array![[0i64, 1, 2], [1, 0, 1], [2, 1, 0]];

        let lambda_time = 1.0;
        let target_period = 2;
        let target_unit = 1;

        // θ_s = exp(−λ_time · |target_period − s|)
        let expected_time = [
            (-lambda_time * 2.0_f64).exp(),
            (-lambda_time * 1.0_f64).exp(),
            (-lambda_time * 0.0_f64).exp(),
        ];

        // With λ_unit = 0, ω_j = 1 for both units ⟹ row sum = θ_s · 2.
        let expected_unit_sum = 2.0;

        let weights = compute_weight_matrix(
            &y.view(),
            &d.view(),
            3,
            2,
            target_unit,
            target_period,
            lambda_time,
            0.0,
            &time_dist.view(),
        );

        // Row sums: Σ_j W[t, j] = θ_t · Σ_j ω_j = θ_t · 2.
        for (t, &expected_t) in expected_time.iter().enumerate().take(3) {
            let row_sum: f64 = weights.row(t).sum();
            let expected_row_sum = expected_t * expected_unit_sum;
            let diff = (row_sum - expected_row_sum).abs();
            assert!(
                diff < 1e-10,
                "Time weight precision error at t={}: expected {}, got {}, diff={}",
                t,
                expected_row_sum,
                row_sum,
                diff
            );
        }

        // Individual entries: W[t, j] = θ_t · ω_j.
        for t in 0..3 {
            let expected_w0 = expected_time[t] * 1.0;
            assert!(
                (weights[[t, 0]] - expected_w0).abs() < 1e-10,
                "Weight precision error at [{}, 0]",
                t
            );
            let expected_w1 = expected_time[t] * 1.0;
            assert!(
                (weights[[t, 1]] - expected_w1).abs() < 1e-10,
                "Weight precision error at [{}, 1]",
                t
            );
        }
    }
}

/// Property-based tests using proptest.
#[cfg(test)]
mod proptests {
    use super::*;
    use ndarray::Array2;
    use proptest::prelude::*;

    /// Generate a random outcome matrix Y of dimension (T × N).
    fn y_matrix_strategy(
        n_periods: usize,
        n_units: usize,
    ) -> impl Strategy<Value = Array2<f64>> {
        prop::collection::vec(-100.0..100.0_f64, n_periods * n_units).prop_map(move |v| {
            Array2::from_shape_vec((n_periods, n_units), v).unwrap()
        })
    }

    /// Treatment matrix with a single treated cell at (T−1, 0).
    fn d_single_treated(n_periods: usize, n_units: usize) -> Array2<f64> {
        let mut d = Array2::<f64>::zeros((n_periods, n_units));
        d[[n_periods - 1, 0]] = 1.0;
        d
    }

    /// Absolute-difference time distance matrix: |t − s|.
    fn time_dist_matrix(n_periods: usize) -> Array2<i64> {
        Array2::from_shape_fn((n_periods, n_periods), |(t, s)| (t as i64 - s as i64).abs())
    }

    proptest! {
        /// All weights are non-negative (exp(·) ≥ 0).
        #[test]
        fn prop_weights_non_negative(
            y in y_matrix_strategy(6, 4),
            lambda_time in 0.0..5.0_f64,
            lambda_unit in 0.0..5.0_f64,
            target_unit in 0..4_usize,
        ) {
            let d = d_single_treated(6, 4);
            let td = time_dist_matrix(6);
            let target_period = 5; // last period (treated)

            let w = compute_weight_matrix(
                &y.view(), &d.view(), 6, 4,
                target_unit, target_period,
                lambda_time, lambda_unit, &td.view(),
            );

            for t in 0..6 {
                for j in 0..4 {
                    prop_assert!(w[[t, j]] >= 0.0,
                        "Weight [{}, {}] = {} should be >= 0", t, j, w[[t, j]]);
                }
            }
        }

        /// Rank-1 (outer product) structure: W[t1,j]/W[t2,j] is constant
        /// across all j with non-zero weight.
        #[test]
        fn prop_weights_outer_product_structure(
            y in y_matrix_strategy(6, 4),
            lambda_time in 0.01..3.0_f64,
            lambda_unit in 0.01..3.0_f64,
        ) {
            let d = d_single_treated(6, 4);
            let td = time_dist_matrix(6);
            let target_unit = 0;
            let target_period = 5;

            let w = compute_weight_matrix(
                &y.view(), &d.view(), 6, 4,
                target_unit, target_period,
                lambda_time, lambda_unit, &td.view(),
            );

            let mut ratio: Option<f64> = None;
            for j in 0..4 {
                if w[[0, j]] > 1e-15 && w[[1, j]] > 1e-15 {
                    let r = w[[0, j]] / w[[1, j]];
                    if let Some(prev_r) = ratio {
                        prop_assert!((r - prev_r).abs() < 1e-10,
                            "Rank-1 structure violated: ratio {} vs {}", r, prev_r);
                    }
                    ratio = Some(r);
                }
            }
        }

        /// Exponential decay: θ_{s1}/θ_{s2} = exp(λ_time · (|t−s2| − |t−s1|)).
        #[test]
        fn prop_time_weight_exponential_decay(
            lambda_time in 0.01..5.0_f64,
        ) {
            let y = Array2::from_shape_fn((5, 3), |(t, i)| (t as f64) + (i as f64) * 0.1);
            let d = Array2::<f64>::zeros((5, 3));
            let td = time_dist_matrix(5);
            let target_unit = 0;
            let target_period = 4;

            let w = compute_weight_matrix(
                &y.view(), &d.view(), 5, 3,
                target_unit, target_period,
                lambda_time, 0.0, &td.view(),
            );

            // With λ_unit = 0, ω_j = 1 ⟹ W[t, j] = θ_t.
            // Ratio W[3,j]/W[2,j] = exp(−λ·1)/exp(−λ·2) = exp(λ).
            for j in 1..3 {
                if w[[3, j]] > 1e-15 && w[[2, j]] > 1e-15 {
                    let ratio = w[[3, j]] / w[[2, j]];
                    let expected_ratio = lambda_time.exp();
                    prop_assert!((ratio - expected_ratio).abs() < 1e-8,
                        "Time decay ratio {} vs expected {}", ratio, expected_ratio);
                }
            }
        }

        /// W[target_period, target_unit] = θ_t(0) · ω_i = 1.0 · 1.0 = 1.0.
        #[test]
        fn prop_target_unit_weight(
            y in y_matrix_strategy(5, 3),
            lambda_time in 0.0..3.0_f64,
            lambda_unit in 0.0..3.0_f64,
        ) {
            let d = Array2::<f64>::zeros((5, 3));
            let td = time_dist_matrix(5);
            let target_unit = 1;
            let target_period = 4;

            let w = compute_weight_matrix(
                &y.view(), &d.view(), 5, 3,
                target_unit, target_period,
                lambda_time, lambda_unit, &td.view(),
            );

            // θ_{target_period} = exp(0) = 1, ω_{target_unit} = 1 ⟹ W = 1.
            prop_assert!((w[[target_period, target_unit]] - 1.0).abs() < 1e-10,
                "W[target_period, target_unit] should be 1.0, got {}",
                w[[target_period, target_unit]]);
        }
    }
}
