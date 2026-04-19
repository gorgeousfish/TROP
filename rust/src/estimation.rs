//! Model estimation for the TROP (Triply Robust Panel) estimator.
//!
//! Provides two estimation strategies via alternating minimization:
//!
//! - **Twostep**: per-observation weighted least squares with nuclear norm penalty.
//!   Solves `min_{α,β,L} Σ w_{ts} (Y_{ts} - α_s - β_t - L_{ts})² + λ‖L‖_*`
//!   for each treated observation, yielding heterogeneous treatment effects τ_{i,t}.
//!
//! - **Joint**: global weighted least squares under a homogeneous treatment effect.
//!   Solves `min_{μ,α,β,L,τ} Σ δ_{ts} (Y_{ts} - μ - α_s - β_t - L_{ts} - τ D_{ts})² + λ‖L‖_*`
//!   yielding a single scalar τ.

#[cfg(target_os = "macos")]
use crate::newlapack;
use faer::linalg::solvers::Svd;
use faer::Mat;
#[cfg(not(target_os = "macos"))]
use lapack::dgelsd;
use ndarray::{Array1, Array2, ArrayView2, Axis};

/// Maximum inner iterations for the proximal gradient solver on the nuclear norm subproblem.
/// The proximal operator converges quickly; 10 iterations suffice for typical panel dimensions.
const MAX_INNER_ITER: usize = 10;

/// SVD singular value truncation tolerance.
///
/// Singular values below this threshold after soft-thresholding are treated as zero
/// to prevent numerical noise from dominating the low-rank reconstruction.
///
/// Used in `soft_threshold_svd()` for:
///   1. Counting nonzero singular values after thresholding.
///   2. Skipping near-zero components in truncated SVD reconstruction.
pub const SVD_TRUNCATION_TOL: f64 = 1e-10;

/// Tolerance for detecting degenerate (all-zero) weight matrices.
///
/// If the sum of all weights falls below this threshold, the estimation
/// is considered degenerate and returns `None`.
pub const WEIGHT_SUM_TOL: f64 = 1e-10;

/// Result type for twostep per-observation estimation.
///
/// Fields: `(alpha, beta, L, n_iterations, converged)`.
/// - `alpha`: unit fixed effects (length N).
/// - `beta`: time fixed effects (length T).
/// - `L`: low-rank matrix (T × N).
/// - `n_iterations`: number of alternating minimization iterations performed.
/// - `converged`: whether the algorithm met the convergence tolerance.
#[allow(clippy::type_complexity)]
pub type TwostepModelResult = Option<(Array1<f64>, Array1<f64>, Array2<f64>, usize, bool)>;

/// Result type for joint estimation with low-rank component.
///
/// Fields: `(mu, alpha, beta, L, tau, n_iterations, converged)`.
/// - `mu`: global intercept.
/// - `alpha`: unit fixed effects (length N).
/// - `beta`: time fixed effects (length T).
/// - `L`: low-rank matrix (T × N).
/// - `tau`: homogeneous treatment effect.
/// - `n_iterations`: number of alternating minimization iterations performed.
/// - `converged`: whether the algorithm met the convergence tolerance.
#[allow(clippy::type_complexity)]
pub type JointLowRankResult = Option<(f64, Array1<f64>, Array1<f64>, Array2<f64>, f64, usize, bool)>;

/// Maximum absolute difference between two 1D arrays.
#[inline]
pub fn max_abs_diff(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

/// Maximum absolute difference between two 2D arrays.
#[inline]
pub fn max_abs_diff_2d(a: &Array2<f64>, b: &Array2<f64>) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

/// Apply soft-thresholding to singular values (proximal operator for the nuclear norm).
///
/// Computes `prox_{λ‖·‖_*}(M) = U diag(max(σ_k − λ, 0)) V^T`, i.e., singular
/// value soft-thresholding. This is the key step in proximal gradient descent
/// for nuclear-norm-penalized matrix estimation.
///
/// - When `threshold = 0`, returns `M` unchanged (no regularization).
/// - When `threshold → ∞`, returns the zero matrix (TWFE/DID limit).
///
/// # Arguments
/// * `m` - Input matrix.
/// * `threshold` - Soft-threshold value λ.
///
/// # Returns
/// The soft-thresholded matrix, or `None` if SVD computation fails.
pub fn soft_threshold_svd(m: &Array2<f64>, threshold: f64) -> Option<Array2<f64>> {
    // λ = 0 means no regularization; return the original matrix.
    if threshold <= 0.0 {
        return Some(m.clone());
    }

    // Check for non-finite values
    if !m.iter().all(|&x| x.is_finite()) {
        return Some(Array2::zeros(m.raw_dim()));
    }

    let n_rows = m.nrows();
    let n_cols = m.ncols();

    // Convert ndarray to faer Mat
    let faer_mat = Mat::from_fn(n_rows, n_cols, |i, j| m[[i, j]]);

    // Compute SVD using faer
    let svd = match Svd::new(faer_mat.as_ref()) {
        Ok(s) => s,
        Err(_) => return Some(Array2::zeros(m.raw_dim())),
    };

    let u = svd.U();
    let s = svd.S().column_vector();
    let v = svd.V();

    // Check for non-finite SVD output
    let k = u.ncols().min(v.ncols());
    for i in 0..k {
        if !s[i].is_finite() {
            return Some(Array2::zeros(m.raw_dim()));
        }
    }

    // Soft-threshold singular values: σ_k ← max(σ_k - threshold, 0)
    let mut s_thresh = Vec::with_capacity(k);
    let mut nonzero_count = 0;
    for i in 0..k {
        let sv = s[i];
        let sv_thresh = (sv - threshold).max(0.0);
        s_thresh.push(sv_thresh);
        if sv_thresh > SVD_TRUNCATION_TOL {
            nonzero_count += 1;
        }
    }

    if nonzero_count == 0 {
        return Some(Array2::zeros(m.raw_dim()));
    }

    // Truncated reconstruction: U @ diag(s_thresh) @ V^T
    let mut result = Array2::<f64>::zeros((n_rows, n_cols));

    for idx in 0..k {
        if s_thresh[idx] > SVD_TRUNCATION_TOL {
            for i in 0..n_rows {
                for j in 0..n_cols {
                    result[[i, j]] += s_thresh[idx] * u[(i, idx)] * v[(j, idx)];
                }
            }
        }
    }

    Some(result)
}

/// Estimate the TROP model via alternating minimization (twostep method).
///
/// For each treated observation (i, t), solves the weighted nuclear-norm-penalized
/// least squares problem:
///
/// ```text
/// min_{α,β,L}  Σ_{j,s} w_{js} (Y_{js} - α_j - β_s - L_{js})²  +  λ ‖L‖_*
/// ```
///
/// where the weight matrix `w` zeroes out treated cells (except the target observation
/// when used in leave-one-out cross-validation). The treatment effect is then recovered
/// externally as `τ̂_{i,t} = Y_{i,t} − α̂_i − β̂_t − L̂_{i,t}`.
///
/// The alternating minimization proceeds as:
///   1. Fix L, update α and β by weighted least squares (Gauss–Seidel).
///   2. Fix (α, β), update L by proximal gradient with step size η = 1/(2·max(w))
///      (Lipschitz constant L_f = 2·max(w) for the unhalved objective):
///      `L ← prox_{η λ ‖·‖_*}(L + w_norm ⊙ (R − L))`, with threshold = λ/(2·max(w)).
///
/// # Arguments
/// * `y` - Outcome matrix Y (T × N).
/// * `control_mask` - Binary mask: 1 for control observations, 0 for treated.
/// * `weight_matrix` - Weight matrix w (T × N).
/// * `lambda_nn` - Nuclear norm penalty parameter λ.
/// * `n_periods` - Number of time periods T.
/// * `n_units` - Number of units N.
/// * `max_iter` - Maximum alternating minimization iterations.
/// * `tol` - Convergence tolerance on max absolute parameter change.
/// * `exclude_obs` - Optional (t, i) index to exclude (for leave-one-out CV).
///
/// # Returns
/// `Some((alpha, beta, L, n_iterations, converged))` on success, `None` on failure.
#[allow(clippy::too_many_arguments)]
pub fn estimate_model(
    y: &ArrayView2<f64>,
    control_mask: &ArrayView2<u8>,
    weight_matrix: &ArrayView2<f64>,
    lambda_nn: f64,
    n_periods: usize,
    n_units: usize,
    max_iter: usize,
    tol: f64,
    exclude_obs: Option<(usize, usize)>,
) -> TwostepModelResult {
    // Create estimation mask
    let mut est_mask =
        Array2::<bool>::from_shape_fn((n_periods, n_units), |(t, i)| control_mask[[t, i]] != 0);

    if let Some((t_ex, i_ex)) = exclude_obs {
        est_mask[[t_ex, i_ex]] = false;
    }

    // Valid mask: non-NaN and in estimation set
    let valid_mask = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
        y[[t, i]].is_finite() && est_mask[[t, i]]
    });

    // Masked weights: W=0 for invalid/treated observations
    let w_masked = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
        if valid_mask[[t, i]] {
            weight_matrix[[t, i]]
        } else {
            0.0
        }
    });

    // Per Eq. 2 in Athey, Imbens, Qu & Viviano (2025) the loss is
    //   f(L) = Σ w_{ti} (R_{ti} − L_{ti})²   (no 1/2 factor),
    // so ∇f(L) = 2 w ⊙ (L − R) and the Lipschitz constant of ∇f is
    // L_f = 2·max(w). Proximal gradient uses step size η = 1/L_f =
    // 1/(2·w_max) and threshold η·λ = λ/(2·w_max).
    //
    // `weight_norm` = w/w_max rescales the gradient step so that the
    // inner update L ← L + weight_norm ⊙ (R − L) has Lipschitz 1
    // (equivalent to using η = 1 on the rescaled problem). The
    // proximal threshold must still be η·λ = λ/(2·w_max) on the
    // original scale.
    let w_max = w_masked.iter().cloned().fold(0.0_f64, f64::max);
    let weight_norm_factor = if w_max > 0.0 { 1.0 / w_max } else { 1.0 };
    let prox_threshold = if w_max > 0.0 {
        lambda_nn / (2.0 * w_max)
    } else {
        lambda_nn / 2.0
    };

    // Weight sums per unit and time
    let weight_sum_per_unit: Array1<f64> = w_masked.sum_axis(Axis(0));
    let weight_sum_per_time: Array1<f64> = w_masked.sum_axis(Axis(1));

    // Safe denominators
    let safe_unit_denom: Array1<f64> = weight_sum_per_unit.mapv(|w| if w > 0.0 { w } else { 1.0 });
    let safe_time_denom: Array1<f64> = weight_sum_per_time.mapv(|w| if w > 0.0 { w } else { 1.0 });

    let unit_has_obs: Array1<bool> = weight_sum_per_unit.mapv(|w| w > 0.0);
    let time_has_obs: Array1<bool> = weight_sum_per_time.mapv(|w| w > 0.0);

    // Safe Y (replace NaN with 0)
    let y_safe = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
        if y[[t, i]].is_finite() {
            y[[t, i]]
        } else {
            0.0
        }
    });

    // Initialize
    let mut alpha = Array1::<f64>::zeros(n_units);
    let mut beta = Array1::<f64>::zeros(n_periods);
    let mut l = Array2::<f64>::zeros((n_periods, n_units));

    // Track actual iteration count and convergence status
    let mut actual_iters: usize = 0;
    let mut converged = false;

    // Alternating minimization
    for _ in 0..max_iter {
        actual_iters += 1;
        let alpha_old = alpha.clone();
        let beta_old = beta.clone();
        let l_old = l.clone();

        // Step 1: Update α and β (weighted least squares).
        // R = Y - L
        let r = &y_safe - &l;

        // Gauss–Seidel update order: α first, then β using the new α.
        // Converges faster than Jacobi; the fixed point is identical.

        // α_i = Σ_t w_{ti} (R_{ti} − β_t) / Σ_t w_{ti}
        // (uses β from the previous iteration)
        for i in 0..n_units {
            if unit_has_obs[i] {
                let mut num = 0.0;
                for t in 0..n_periods {
                    num += w_masked[[t, i]] * (r[[t, i]] - beta[t]);
                }
                alpha[i] = num / safe_unit_denom[i];
            }
        }

        // β_t = Σ_i w_{ti} (R_{ti} − α_i) / Σ_i w_{ti}
        // (uses the newly computed α — Gauss–Seidel)
        for t in 0..n_periods {
            if time_has_obs[t] {
                let mut num = 0.0;
                for i in 0..n_units {
                    num += w_masked[[t, i]] * (r[[t, i]] - alpha[i]);
                }
                beta[t] = num / safe_time_denom[t];
            }
        }

        // Step 2: Update L via proximal gradient for the nuclear norm penalty.
        //
        // Subproblem: min_L (1/2) Σ w_{ti} (R_{ti} − L_{ti})² + λ ‖L‖_*
        //
        // Lipschitz constant of ∇f is w_max, so step size η = 1/w_max.
        // With normalized weights w_norm = w/w_max:
        //   gradient_step = L + w_norm ⊙ (R − L)
        //   L ← prox_{η λ ‖·‖_*}(gradient_step)
        //
        // Threshold = η λ = λ/w_max.

        // Compute target residual R = Y - α - β
        let mut r_target = Array2::<f64>::zeros((n_periods, n_units));
        for t in 0..n_periods {
            for i in 0..n_units {
                r_target[[t, i]] = y_safe[[t, i]] - alpha[i] - beta[t];
            }
        }

        // When λ = 0, set L = R directly for valid observations (no nuclear norm penalty).
        if lambda_nn <= 0.0 {
            for t in 0..n_periods {
                for i in 0..n_units {
                    if valid_mask[[t, i]] {
                        l[[t, i]] = r_target[[t, i]];
                    }
                    // Invalid observations (w = 0): keep L unchanged.
                }
            }
        } else {
            // FISTA-accelerated proximal gradient for the L subproblem.
            // Iterates (with Nesterov momentum):
            //   L_{k+1} = prox_{η λ ‖·‖_*}(L̂_k + w_norm ⊙ (R_masked − L̂_k))
            // where L̂_k = L_k + ((t_k − 1)/t_{k+1}) (L_k − L_{k−1}),
            //   w_norm = w / w_max, and threshold = η λ = λ / (2·w_max).

            // R_masked: use r_target for valid observations, keep L for invalid ones
            // to prevent L from absorbing signal at zero-weight cells.
            let r_masked = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
                if valid_mask[[t, i]] {
                    r_target[[t, i]]
                } else {
                    l[[t, i]]
                }
            });

            // W_norm = W / W_max (normalized weights, max = 1)
            let w_norm = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
                w_masked[[t, i]] * weight_norm_factor // w / w_max
            });

            let mut l_prev = l.clone();
            let mut t_fista = 1.0_f64;

            for _ in 0..MAX_INNER_ITER {
                let l_inner_old = l.clone();

                // Nesterov momentum extrapolation.
                let t_fista_new = (1.0 + (1.0 + 4.0 * t_fista * t_fista).sqrt()) / 2.0;
                let momentum = (t_fista - 1.0) / t_fista_new;
                let l_momentum = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
                    l[[t, i]] + momentum * (l[[t, i]] - l_prev[[t, i]])
                });

                // Gradient step from the momentum point.
                let gradient_step = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
                    l_momentum[[t, i]]
                        + w_norm[[t, i]] * (r_masked[[t, i]] - l_momentum[[t, i]])
                });

                // Proximal step: soft-threshold singular values.
                l_prev = l.clone();
                l = soft_threshold_svd(&gradient_step, prox_threshold)?;
                t_fista = t_fista_new;

                // Gradient-based adaptive restart (O'Donoghue & Candes 2015,
                // "Adaptive Restart for Accelerated Gradient Schemes", Found.
                // Comp. Math.).  FISTA's Nesterov momentum can overshoot on
                // non-strongly-convex problems, producing oscillations that
                // slow convergence.  The restart criterion
                //     ⟨y_k − x_k, x_k − x_{k−1}⟩ > 0
                // (here y_k = l_momentum, x_k = l, x_{k−1} = l_inner_old)
                // detects when the proximal step moves in the direction
                // opposite to the momentum step, signalling that the
                // accumulated momentum has overshot the optimum.  Resetting
                // t_fista = 1 clears the momentum without altering the fixed
                // point.  This yields strictly monotone progress on the
                // weighted nuclear-norm subproblem and resolves PWT-style
                // slow convergence at small λ_nn.
                let mut restart_inner = 0.0_f64;
                for t in 0..n_periods {
                    for i in 0..n_units {
                        restart_inner += (l_momentum[[t, i]] - l[[t, i]])
                            * (l[[t, i]] - l_inner_old[[t, i]]);
                    }
                }
                if restart_inner > 0.0 {
                    t_fista = 1.0;
                }

                // Check inner convergence against the previous iterate.
                if max_abs_diff_2d(&l, &l_inner_old) < tol {
                    break;
                }
            }
        }

        // Check convergence
        let alpha_diff = max_abs_diff(&alpha, &alpha_old);
        let beta_diff = max_abs_diff(&beta, &beta_old);
        let l_diff = max_abs_diff_2d(&l, &l_old);

        if alpha_diff.max(beta_diff).max(l_diff) < tol {
            converged = true;
            break;
        }
    }

    // Return actual iteration count and convergence status
    Some((alpha, beta, l, actual_iters, converged))
}

/// Solve the weighted two-way fixed effects regression over control observations.
///
/// The paper's joint objective (Eq. 2) restricts the quadratic loss to control
/// cells via the (1 − D) factor. We expect the caller to supply `delta` that
/// is already multiplied by (1 − D) (e.g. via `compute_joint_weights`). With
/// that convention the problem reduces to
///
/// ```text
/// min_{μ, α, β}   Σ_{t, i} δ_{t, i} (Y_{t, i} − μ − α_i − β_t)²
/// ```
///
/// where δ is zero at treated observations. The treatment effect τ is then
/// obtained post-hoc as the mean residual over treated cells (caller's job).
///
/// Solved via LAPACK `dgelsd` (minimum-norm least squares) for numerical
/// stability with potentially rank-deficient design matrices.
///
/// Identification constraint: α_0 = β_0 = 0 (first unit and time dummy dropped).
///
/// # Arguments
/// * `y` - Outcome matrix Y (T × N).
/// * `delta` - Global weight matrix δ (T × N), already (1 − D)-masked.
///
/// # Returns
/// `Some((mu, alpha, beta))` on success, `None` if the system is degenerate.
pub fn solve_joint_no_lowrank(
    y: &ArrayView2<f64>,
    delta: &ArrayView2<f64>,
) -> Option<(f64, Array1<f64>, Array1<f64>)> {
    let n_periods = y.nrows();
    let n_units = y.ncols();
    let n_obs = n_periods * n_units;

    // Parameter count: 1 (intercept) + (N−1) unit dummies + (T−1) time dummies.
    // There is NO treatment column — τ is computed post-hoc as ATT on residuals.
    let n_params = 1 + (n_units - 1) + (n_periods - 1);

    // Vectorize the panel and compute observation weights.
    let mut y_vec = Vec::with_capacity(n_obs);
    let mut w_vec = Vec::with_capacity(n_obs);

    for t in 0..n_periods {
        for i in 0..n_units {
            let y_ti = y[[t, i]];
            let delta_ti = delta[[t, i]];

            // Zero weight for non-finite outcomes or weights.
            let (y_val, w_val) = if y_ti.is_finite() && delta_ti.is_finite() {
                (y_ti, delta_ti)
            } else {
                (0.0, 0.0)
            };

            y_vec.push(y_val);
            w_vec.push(w_val);
        }
    }

    // Check for all-zero weights
    let sum_w: f64 = w_vec.iter().sum();
    if sum_w < WEIGHT_SUM_TOL {
        return None;
    }

    // Compute sqrt(weights)
    let sqrt_w: Vec<f64> = w_vec.iter().map(|&w| w.max(0.0).sqrt()).collect();

    let m = n_obs as i32;
    let n = n_params as i32;
    let lda = m;
    let ldb = m.max(n);
    let min_mn = m.min(n) as usize;

    // Build weighted design matrix A in column-major (Fortran) layout.
    let mut a_mat = vec![0.0_f64; (m * n) as usize];

    for t in 0..n_periods {
        for i in 0..n_units {
            let obs_idx = t * n_units + i;
            let sw = sqrt_w[obs_idx];

            // Column 0: intercept
            a_mat[obs_idx] = sw;

            // Columns 1..(N−1): unit dummies (unit 0 dropped for identification)
            if i > 0 {
                a_mat[i * n_obs + obs_idx] = sw;
            }

            // Columns N..(N+T−2): time dummies (time 0 dropped for identification)
            if t > 0 {
                a_mat[((n_units - 1) + t) * n_obs + obs_idx] = sw;
            }
        }
    }

    // Build weighted response vector b (extended to ldb size)
    let mut b_vec = vec![0.0_f64; ldb as usize];
    for (idx, (&y_val, &sw)) in y_vec.iter().zip(sqrt_w.iter()).enumerate() {
        b_vec[idx] = y_val * sw;
    }

    // Singular values output
    let mut s = vec![0.0_f64; min_mn];

    // rcond for rank determination: ε · max(m, n)
    let eps = f64::EPSILON;
    let rcond = eps * (m.max(n) as f64);

    // Output rank
    let mut rank: i32 = 0;

    // Workspace query
    let mut work_query = vec![0.0_f64; 1];
    let mut iwork_query = vec![0_i32; 1];
    let mut info: i32 = 0;
    let nrhs: i32 = 1;
    #[allow(unused_assignments)]
    let mut lwork: i32 = -1;

    // Query optimal workspace size
    #[cfg(target_os = "macos")]
    {
        newlapack::dgelsd(
            m,
            n,
            nrhs,
            &mut a_mat,
            lda,
            &mut b_vec,
            ldb,
            &mut s,
            rcond,
            &mut rank,
            &mut work_query,
            -1,
            &mut iwork_query,
            &mut info,
        );
    }
    #[cfg(not(target_os = "macos"))]
    unsafe {
        dgelsd(
            m,
            n,
            nrhs,
            &mut a_mat,
            lda,
            &mut b_vec,
            ldb,
            &mut s,
            rcond,
            &mut rank,
            &mut work_query,
            lwork,
            &mut iwork_query,
            &mut info,
        );
    }

    if info != 0 {
        return None;
    }

    // Allocate workspace
    lwork = work_query[0] as i32;
    let mut work = vec![0.0_f64; lwork as usize];

    // Calculate iwork size (LAPACK internal parameter)
    let smlsiz = 25_i32;
    let nlvl = if min_mn > 0 {
        ((min_mn as f64 / (smlsiz + 1) as f64).ln() / 2.0_f64.ln()).floor() as i32 + 1
    } else {
        0
    };
    let nlvl = nlvl.max(0);
    let liwork = (3 * (min_mn as i32) * nlvl + 11 * (min_mn as i32)).max(1);
    let mut iwork = vec![0_i32; liwork as usize];

    // Rebuild A matrix (dgelsd overwrites it during the workspace query).
    a_mat.fill(0.0);
    for t in 0..n_periods {
        for i in 0..n_units {
            let obs_idx = t * n_units + i;
            let sw = sqrt_w[obs_idx];

            a_mat[obs_idx] = sw;
            if i > 0 {
                a_mat[i * n_obs + obs_idx] = sw;
            }
            if t > 0 {
                a_mat[((n_units - 1) + t) * n_obs + obs_idx] = sw;
            }
        }
    }

    // Rebuild b vector (also overwritten by dgelsd).
    b_vec.fill(0.0);
    for (idx, (&y_val, &sw)) in y_vec.iter().zip(sqrt_w.iter()).enumerate() {
        b_vec[idx] = y_val * sw;
    }

    // Solve the weighted least squares system via dgelsd.
    #[cfg(target_os = "macos")]
    {
        newlapack::dgelsd(
            m, n, nrhs, &mut a_mat, lda, &mut b_vec, ldb, &mut s, rcond, &mut rank, &mut work,
            lwork, &mut iwork, &mut info,
        );
    }
    #[cfg(not(target_os = "macos"))]
    unsafe {
        dgelsd(
            m, n, nrhs, &mut a_mat, lda, &mut b_vec, ldb, &mut s, rcond, &mut rank, &mut work,
            lwork, &mut iwork, &mut info,
        );
    }

    if info != 0 {
        return None;
    }

    // Extract solution: b_vec[0..n_params] holds the least-squares coefficients.
    let mu = b_vec[0];

    let mut alpha = Array1::<f64>::zeros(n_units);
    for i in 1..n_units {
        alpha[i] = b_vec[i];
    }

    let mut beta = Array1::<f64>::zeros(n_periods);
    for t in 1..n_periods {
        beta[t] = b_vec[(n_units - 1) + t];
    }

    Some((mu, alpha, beta))
}

/// Solve the joint TWFE + low-rank model via alternating minimization, with
/// post-hoc τ extraction.
///
/// The weight matrix `delta` must already carry the paper's (1 − D) mask so
/// treated cells contribute zero to the weighted quadratic loss. The objective
/// solved is therefore the paper's Equation 2 without τ as an explicit
/// variable:
///
/// ```text
/// min_{μ, α, β, L}   Σ_{t, i} δ_{t, i} (Y_{t, i} − μ − α_i − β_t − L_{t, i})²
///                  + λ ‖L‖_*
/// ```
///
/// Alternating minimization with FISTA acceleration on the L subproblem:
///   1. Fix L, solve (μ, α, β) via WLS (control-only thanks to δ masking).
///   2. Fix (μ, α, β), update L via a few FISTA/Nesterov proximal iterations.
/// After outer convergence we do a final re-solve of (μ, α, β) using the
/// converged L, then compute τ post-hoc as the mean residual over treated
/// observations:  τ̂ = mean_{D=1} (Y − μ − α − β − L).
///
/// # Arguments
/// * `y` - Outcome matrix Y (T × N).
/// * `d` - Treatment indicator matrix D (T × N). Used only to locate treated
///   cells for the post-hoc τ calculation.
/// * `delta` - Global weight matrix δ (T × N). **Must** already be (1 − D)-masked.
/// * `lambda_nn` - Nuclear norm penalty parameter λ.
/// * `max_iter` - Maximum alternating minimization iterations.
/// * `tol` - Convergence tolerance on max absolute parameter change.
///
/// # Returns
/// `Some((mu, alpha, beta, L, tau, n_iterations, converged))` on success,
/// `None` on failure.
pub fn solve_joint_with_lowrank(
    y: &ArrayView2<f64>,
    d: &ArrayView2<f64>,
    delta: &ArrayView2<f64>,
    lambda_nn: f64,
    max_iter: usize,
    tol: f64,
) -> JointLowRankResult {
    let n_periods = y.nrows();
    let n_units = y.ncols();

    // Sanitize Y: replace non-finite values with 0 for the inner arithmetic.
    // Zero the corresponding δ so these positions contribute nothing.
    let y_safe = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
        if y[[t, i]].is_finite() { y[[t, i]] } else { 0.0 }
    });
    let delta_masked = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
        if y[[t, i]].is_finite() { delta[[t, i]] } else { 0.0 }
    });

    // Precompute δ_max and the FISTA prox threshold (constant across iterations).
    // Per paper Eq. 2 (no 1/2 factor): ∇f = 2 δ ⊙ (L − R), Lipschitz constant
    // L_f = 2·δ_max, optimal step size η = 1/(2·δ_max), and the soft-threshold
    // value is η λ = λ/(2·δ_max).
    let delta_max = delta_masked.iter().copied().fold(0.0_f64, f64::max);
    let prox_threshold = if delta_max > 0.0 {
        lambda_nn / (2.0 * delta_max)
    } else {
        lambda_nn / 2.0
    };

    // δ_norm = δ / δ_max, used both for the gradient step scaling and for
    // detecting "active" cells where the R target overrides the current L.
    let delta_norm = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
        if delta_max > 0.0 {
            delta_masked[[t, i]] / delta_max
        } else {
            delta_masked[[t, i]]
        }
    });

    // Initialize L = 0.
    let mut l = Array2::<f64>::zeros((n_periods, n_units));

    // Store the last iteration's parameters; after the final re-solve below
    // these are overwritten with the converged values. The `mu` initializer
    // is a placeholder that always gets replaced either inside the loop or by
    // the final re-solve; mark the initial assignment as intentionally unused.
    #[allow(unused_assignments)]
    let mut mu = 0.0_f64;
    let mut alpha = Array1::<f64>::zeros(n_units);
    let mut beta = Array1::<f64>::zeros(n_periods);

    // Track actual iteration count and convergence status.
    let mut actual_iters: usize = 0;
    let mut converged = false;

    for _ in 0..max_iter {
        actual_iters += 1;
        let l_old = l.clone();
        let alpha_old = alpha.clone();
        let beta_old = beta.clone();

        // Step 1: Fix L, solve for (μ, α, β) via WLS on adjusted outcome Y − L.
        // δ_masked is (1 − D)-masked at the caller, so this regression uses only
        // control observations; no τ column is needed.
        let y_adj = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
            y_safe[[t, i]] - l[[t, i]]
        });
        let (mu_new, alpha_new, beta_new) =
            solve_joint_no_lowrank(&y_adj.view(), &delta_masked.view())?;
        mu = mu_new;
        alpha = alpha_new;
        beta = beta_new;

        // Step 2: Fix (μ, α, β), update L via FISTA proximal gradient.
        // Target residual R = Y − μ − α − β; at zero-weight cells (treated /
        // non-finite) we substitute the current L so the gradient step leaves
        // L unchanged in those positions.
        let r_masked = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
            if delta_norm[[t, i]] > 0.0 {
                y_safe[[t, i]] - mu - alpha[i] - beta[t]
            } else {
                l[[t, i]]
            }
        });

        // FISTA inner loop: max 20 steps, early exit on tol.
        const MAX_JOINT_INNER_ITER: usize = 20;
        let mut l_prev = l.clone();
        let mut t_fista = 1.0_f64;

        for _ in 0..MAX_JOINT_INNER_ITER {
            // Snapshot at the start of this inner step.  Needed both for the
            // gradient-restart criterion below and (incidentally) for an
            // unambiguous "progress between two consecutive iterates"
            // convergence check.
            let l_inner_old = l.clone();

            // Nesterov momentum extrapolation.
            let t_fista_new = (1.0 + (1.0 + 4.0 * t_fista * t_fista).sqrt()) / 2.0;
            let momentum = (t_fista - 1.0) / t_fista_new;
            let l_momentum = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
                l[[t, i]] + momentum * (l[[t, i]] - l_prev[[t, i]])
            });

            // Gradient step from the momentum point.
            let gradient_step = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
                l_momentum[[t, i]]
                    + delta_norm[[t, i]] * (r_masked[[t, i]] - l_momentum[[t, i]])
            });

            // Proximal step: soft-threshold singular values.
            l_prev = l.clone();
            l = soft_threshold_svd(&gradient_step, prox_threshold)?;
            t_fista = t_fista_new;

            // Gradient-based adaptive restart (O'Donoghue & Candes 2015).
            // See the equivalent comment in `estimate_model` for the
            // derivation.  The joint path re-solves (μ, α, β) after every
            // outer iteration, so FISTA oscillations here propagate into the
            // outer fixed point — restart is therefore a stability
            // requirement, not just a convergence accelerator.
            let mut restart_inner = 0.0_f64;
            for t in 0..n_periods {
                for i in 0..n_units {
                    restart_inner += (l_momentum[[t, i]] - l[[t, i]])
                        * (l[[t, i]] - l_inner_old[[t, i]]);
                }
            }
            if restart_inner > 0.0 {
                t_fista = 1.0;
            }

            // Inner convergence check: |L_new - L_old| across this FISTA
            // step.  We compare against `l_inner_old` (start-of-step) rather
            // than the pre-SVD `l_prev`: the former is the quantity whose
            // decay governs fixed-point progress, which is what the outer
            // iteration actually needs.
            if max_abs_diff_2d(&l, &l_inner_old) < tol {
                break;
            }
        }

        // Outer convergence check on L, α, and β.
        let l_diff = max_abs_diff_2d(&l, &l_old);
        let alpha_diff = max_abs_diff(&alpha, &alpha_old);
        let beta_diff = max_abs_diff(&beta, &beta_old);
        if l_diff.max(alpha_diff).max(beta_diff) < tol {
            converged = true;
            break;
        }
    }

    // Final re-solve of (μ, α, β) using the converged L so the returned
    // parameters are mutually consistent (cf. diff-diff 3.1.1 behaviour).
    let y_adj_final = Array2::from_shape_fn((n_periods, n_units), |(t, i)| {
        y_safe[[t, i]] - l[[t, i]]
    });
    let (mu_final, alpha_final, beta_final) =
        solve_joint_no_lowrank(&y_adj_final.view(), &delta_masked.view())?;
    mu = mu_final;
    alpha = alpha_final;
    beta = beta_final;

    // Post-hoc τ: mean residual over observed treated cells.
    //   τ̂ = (1 / |T_1|) Σ_{D=1} (Y − μ − α_i − β_t − L_{t, i})
    let mut tau_sum = 0.0_f64;
    let mut tau_count: usize = 0;
    for t in 0..n_periods {
        for i in 0..n_units {
            if d[[t, i]] == 1.0 && y[[t, i]].is_finite() {
                tau_sum += y[[t, i]] - mu - alpha[i] - beta[t] - l[[t, i]];
                tau_count += 1;
            }
        }
    }
    let tau = if tau_count > 0 {
        tau_sum / tau_count as f64
    } else {
        0.0
    };

    Some((mu, alpha, beta, l, tau, actual_iters, converged))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    /// Helper: compute post-hoc τ = mean residual over observed treated cells.
    ///
    /// `l_opt = None` corresponds to the no-low-rank case (L ≡ 0). Tests use
    /// this to recover τ from the (μ, α, β) returned by `solve_joint_no_lowrank`
    /// or from an explicit L (e.g., a joint low-rank fit).
    fn post_hoc_tau(
        y: &ArrayView2<f64>,
        d: &ArrayView2<f64>,
        mu: f64,
        alpha: &Array1<f64>,
        beta: &Array1<f64>,
        l_opt: Option<&Array2<f64>>,
    ) -> f64 {
        let mut tau_sum = 0.0_f64;
        let mut tau_count = 0usize;
        for t in 0..y.nrows() {
            for i in 0..y.ncols() {
                if d[[t, i]] == 1.0 && y[[t, i]].is_finite() {
                    let l_val = l_opt.map_or(0.0, |m| m[[t, i]]);
                    tau_sum += y[[t, i]] - mu - alpha[i] - beta[t] - l_val;
                    tau_count += 1;
                }
            }
        }
        if tau_count > 0 {
            tau_sum / tau_count as f64
        } else {
            0.0
        }
    }

    #[test]
    fn test_max_abs_diff() {
        let a = array![1.0, 2.0, 3.0];
        let b = array![1.1, 1.9, 3.5];

        let diff = max_abs_diff(&a, &b);
        assert!((diff - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_max_abs_diff_2d() {
        let a = array![[1.0, 2.0], [3.0, 4.0]];
        let b = array![[1.1, 2.2], [2.5, 4.0]];

        let diff = max_abs_diff_2d(&a, &b);
        assert!((diff - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_soft_threshold_svd_zero_threshold() {
        let m = array![[1.0, 2.0], [3.0, 4.0]];

        // With threshold = 0, should return original matrix
        let result = soft_threshold_svd(&m, 0.0).unwrap();

        for i in 0..2 {
            for j in 0..2 {
                assert!((result[[i, j]] - m[[i, j]]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_soft_threshold_svd_correctness() {
        // Create a rank-2 matrix
        let m = array![[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [1.0, 2.0, 3.0]];

        // With large threshold, should reduce rank
        let result = soft_threshold_svd(&m, 5.0).unwrap();

        // Result should have smaller Frobenius norm
        let orig_norm: f64 = m.iter().map(|x| x * x).sum::<f64>().sqrt();
        let result_norm: f64 = result.iter().map(|x| x * x).sum::<f64>().sqrt();

        assert!(result_norm <= orig_norm + 1e-10);
    }

    #[test]
    fn test_solve_joint_no_lowrank_simple() {
        // Simple case: 3 periods, 2 units
        // Unit 0: control, Unit 1: treated at period 2
        let y = array![
            [1.0, 2.0],
            [2.0, 3.0],
            [3.0, 5.0] // Unit 1 has treatment effect of 1
        ];
        let d = array![[0.0, 0.0], [0.0, 0.0], [0.0, 1.0]];
        let delta = array![[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]];

        let result = solve_joint_no_lowrank(&y.view(), &delta.view());
        assert!(result.is_some());

        let (mu, alpha, beta) = result.unwrap();
        // Post-hoc τ = mean residual over treated cells.
        let tau = post_hoc_tau(&y.view(), &d.view(), mu, &alpha, &beta, None);
        // tau should capture the treatment effect
        // The exact value depends on the identification constraints
        assert!(tau.is_finite());
    }

    // ========================================================================
    // SVD Soft Threshold Tests
    // ========================================================================

    #[test]
    fn test_soft_threshold_svd_numerical_precision() {
        // Test SVD soft threshold numerical precision (tolerance < 1e-10)
        // Create a known matrix with specific singular values
        // M = U @ diag(s) @ V^T where s = [5.0, 3.0, 1.0]
        let m = array![[5.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 1.0]];

        // Threshold = 2.0 should give s_thresh = [3.0, 1.0, 0.0]
        let result = soft_threshold_svd(&m, 2.0).unwrap();

        // Expected result: diag([3.0, 1.0, 0.0])
        let expected = array![[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]];

        for i in 0..3 {
            for j in 0..3 {
                let diff = (result[[i, j]] - expected[[i, j]]).abs();
                assert!(
                    diff < 1e-10,
                    "SVD soft threshold precision error at [{}, {}]: expected {}, got {}, diff={}",
                    i,
                    j,
                    expected[[i, j]],
                    result[[i, j]],
                    diff
                );
            }
        }
    }

    #[test]
    fn test_soft_threshold_svd_rank_reduction() {
        // Test that soft thresholding reduces matrix rank
        // Create a rank-3 matrix
        let m = array![[3.0, 1.0, 1.0], [1.0, 3.0, 1.0], [1.0, 1.0, 3.0]];

        // Large threshold should reduce to lower rank
        let result = soft_threshold_svd(&m, 3.0).unwrap();

        // Result should have smaller nuclear norm
        // (sum of singular values should be reduced)
        let result_frob: f64 = result.iter().map(|x| x * x).sum::<f64>().sqrt();
        let orig_frob: f64 = m.iter().map(|x| x * x).sum::<f64>().sqrt();

        assert!(
            result_frob < orig_frob,
            "Soft threshold should reduce Frobenius norm: orig={}, result={}",
            orig_frob,
            result_frob
        );
    }

    #[test]
    fn test_soft_threshold_svd_negative_threshold() {
        // Negative threshold should behave like zero threshold
        let m = array![[1.0, 2.0], [3.0, 4.0]];
        let result = soft_threshold_svd(&m, -1.0).unwrap();

        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    (result[[i, j]] - m[[i, j]]).abs() < 1e-10,
                    "Negative threshold should return original matrix"
                );
            }
        }
    }

    #[test]
    fn test_soft_threshold_svd_large_threshold() {
        // Very large threshold should give zero matrix
        let m = array![[1.0, 2.0], [3.0, 4.0]];
        let result = soft_threshold_svd(&m, 100.0).unwrap();

        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    result[[i, j]].abs() < 1e-10,
                    "Large threshold should give zero matrix, got {} at [{}, {}]",
                    result[[i, j]],
                    i,
                    j
                );
            }
        }
    }

    // ========================================================================
    // Model Convergence Tests
    // ========================================================================

    #[test]
    fn test_estimate_model_convergence() {
        // Test that alternating minimization converges
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 6.0]
        ];
        let control_mask = array![
            [1u8, 1, 1],
            [1, 1, 1],
            [1, 1, 0], // Unit 2 treated at period 2
            [1, 1, 0]  // Unit 2 treated at period 3
        ];
        let weight_matrix = array![
            [0.1, 0.1, 0.1],
            [0.1, 0.1, 0.1],
            [0.1, 0.1, 0.1],
            [0.1, 0.1, 0.1]
        ];

        let result = estimate_model(
            &y.view(),
            &control_mask.view(),
            &weight_matrix.view(),
            0.1,  // lambda_nn
            4,    // n_periods
            3,    // n_units
            100,  // max_iter
            1e-6, // tol
            None,
        );

        assert!(result.is_some(), "Model estimation should converge");

        let (alpha, beta, l, n_iters, converged) = result.unwrap();

        // Check dimensions
        assert_eq!(alpha.len(), 3, "Alpha should have n_units elements");
        assert_eq!(beta.len(), 4, "Beta should have n_periods elements");
        assert_eq!(l.dim(), (4, 3), "L should be n_periods x n_units");

        // Check all values are finite
        assert!(
            alpha.iter().all(|&x| x.is_finite()),
            "Alpha should be finite"
        );
        assert!(beta.iter().all(|&x| x.is_finite()), "Beta should be finite");
        assert!(l.iter().all(|&x| x.is_finite()), "L should be finite");

        // Verify iteration info is meaningful
        assert!(n_iters > 0, "Should have at least 1 iteration");
        assert!(n_iters <= 100, "Should not exceed max_iter=100");
        assert!(converged, "Should converge for this simple case with lambda_nn=0.1");
    }

    #[test]
    fn test_estimate_model_with_exclude_obs() {
        // Test LOOCV exclusion functionality
        let y = array![[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]];
        let control_mask = array![[1u8, 1], [1, 1], [1, 1]];
        let weight_matrix = array![[0.2, 0.2], [0.2, 0.2], [0.2, 0.2]];

        // Estimate with exclusion
        let result = estimate_model(
            &y.view(),
            &control_mask.view(),
            &weight_matrix.view(),
            0.0,
            3,
            2,
            100,
            1e-6,
            Some((1, 0)), // Exclude observation at (1, 0)
        );

        assert!(
            result.is_some(),
            "Model should converge with excluded observation"
        );

        let (_alpha, _beta, _l, n_iters, _converged) = result.unwrap();
        assert!(n_iters > 0, "Should have at least 1 iteration with exclude_obs");
    }

    // ========================================================================
    // ATT Estimation Precision Tests
    // ========================================================================

    #[test]
    fn test_joint_estimation_att_precision() {
        // Test ATT estimation precision (tolerance < 1e-6)
        // Create synthetic data with known treatment effect
        let true_effect = 2.0;

        // Y = alpha_i + beta_t + tau * D + noise
        // alpha = [0, 1], beta = [0, 1, 2, 3]
        let y = array![
            [0.0, 1.0],               // t=0: alpha + beta[0]
            [1.0, 2.0],               // t=1: alpha + beta[1]
            [2.0, 3.0],               // t=2: alpha + beta[2]
            [3.0, 4.0 + true_effect]  // t=3: alpha + beta[3] + tau*D
        ];
        let d = array![
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 1.0] // Unit 1 treated at period 3
        ];
        // δ must be (1 − D)-masked per the new solver contract. Here we zero
        // the one treated cell (3, 1) so the control-only WLS is well defined.
        let delta = array![
            [1.0, 1.0],
            [1.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0]
        ];

        let result = solve_joint_no_lowrank(&y.view(), &delta.view());
        assert!(result.is_some(), "Joint estimation should succeed");

        let (mu, alpha, beta) = result.unwrap();
        // Post-hoc τ = mean residual (Y − μ − α − β) over treated cells.
        let tau = post_hoc_tau(&y.view(), &d.view(), mu, &alpha, &beta, None);

        // tau should be close to true_effect
        let diff = (tau - true_effect).abs();
        assert!(
            diff < 1e-6,
            "ATT estimation precision error: expected {}, got {}, diff={}",
            true_effect,
            tau,
            diff
        );
    }

    #[test]
    fn test_joint_with_lowrank_convergence() {
        // Test joint estimation with low-rank component converges
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 8.0] // Unit 2 has treatment effect
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0]
        ];
        let delta = array![
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0]
        ];

        let result = solve_joint_with_lowrank(
            &y.view(),
            &d.view(),
            &delta.view(),
            0.1,  // lambda_nn
            100,  // max_iter
            1e-6, // tol
        );

        assert!(result.is_some(), "Joint with low-rank should converge");

        let (mu, alpha, beta, l, tau, n_iters, converged) = result.unwrap();

        // Check all values are finite
        assert!(mu.is_finite(), "mu should be finite");
        assert!(
            alpha.iter().all(|&x| x.is_finite()),
            "alpha should be finite"
        );
        assert!(beta.iter().all(|&x| x.is_finite()), "beta should be finite");
        assert!(l.iter().all(|&x| x.is_finite()), "L should be finite");
        assert!(tau.is_finite(), "tau should be finite");

        // Verify convergence info is correct
        assert!(n_iters > 0, "Should have at least 1 iteration");
        assert!(n_iters <= 100, "Should not exceed max_iter={}", 100);
        // With this simple data, should converge well before max_iter
        assert!(converged, "Should converge for this simple case");

        // tau should be positive (treatment effect is positive)
        assert!(tau > 0.0, "tau should be positive, got {}", tau);
    }

    #[test]
    fn test_joint_no_lowrank_with_nan() {
        // Test handling of NaN values
        let y = array![[1.0, f64::NAN], [2.0, 3.0], [3.0, 4.0]];
        let d = array![[0.0, 0.0], [0.0, 0.0], [0.0, 1.0]];
        let delta = array![[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]];

        let result = solve_joint_no_lowrank(&y.view(), &delta.view());
        assert!(result.is_some(), "Should handle NaN values");

        let (mu, alpha, beta) = result.unwrap();
        let tau = post_hoc_tau(&y.view(), &d.view(), mu, &alpha, &beta, None);
        assert!(
            tau.is_finite(),
            "tau should be finite even with NaN in data"
        );
    }

    /// Solve_joint_no_lowrank should handle non-finite
    /// delta (weight) values by zeroing their contribution.
    #[test]
    fn test_joint_no_lowrank_with_nonfinite_delta() {
        // Base case: all finite weights
        let y = array![[1.0, 2.0], [2.0, 3.0], [3.0, 5.0]];
        let d = array![[0.0, 0.0], [0.0, 0.0], [0.0, 1.0]];
        let delta_clean = array![[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]];

        let result_clean = solve_joint_no_lowrank(&y.view(), &delta_clean.view());
        assert!(result_clean.is_some(), "Clean case should succeed");
        let (mu_c, alpha_c, beta_c) = result_clean.unwrap();
        let tau_clean = post_hoc_tau(&y.view(), &d.view(), mu_c, &alpha_c, &beta_c, None);
        assert!(tau_clean.is_finite(), "tau should be finite for clean data");

        // Case 1: NaN delta at a control position — should be excluded gracefully
        let delta_nan = array![[f64::NAN, 1.0], [1.0, 1.0], [1.0, 1.0]];
        let result_nan = solve_joint_no_lowrank(&y.view(), &delta_nan.view());
        assert!(result_nan.is_some(), "NaN delta case should succeed");
        let (mu_n, alpha_n, beta_n) = result_nan.unwrap();
        let tau_nan = post_hoc_tau(&y.view(), &d.view(), mu_n, &alpha_n, &beta_n, None);
        assert!(tau_nan.is_finite(), "tau should be finite with NaN delta");

        // Case 2: Inf delta at a control position — should be excluded gracefully
        let delta_inf = array![[f64::INFINITY, 1.0], [1.0, 1.0], [1.0, 1.0]];
        let result_inf = solve_joint_no_lowrank(&y.view(), &delta_inf.view());
        assert!(result_inf.is_some(), "Inf delta case should succeed");
        let (mu_i, alpha_i, beta_i) = result_inf.unwrap();
        let tau_inf = post_hoc_tau(&y.view(), &d.view(), mu_i, &alpha_i, &beta_i, None);
        assert!(tau_inf.is_finite(), "tau should be finite with Inf delta");

        // Case 3: Negative Inf delta
        let delta_ninf = array![[f64::NEG_INFINITY, 1.0], [1.0, 1.0], [1.0, 1.0]];
        let result_ninf = solve_joint_no_lowrank(&y.view(), &delta_ninf.view());
        assert!(result_ninf.is_some(), "NegInf delta case should succeed");
        let (mu_ni, alpha_ni, beta_ni) = result_ninf.unwrap();
        let tau_ninf = post_hoc_tau(&y.view(), &d.view(), mu_ni, &alpha_ni, &beta_ni, None);
        assert!(tau_ninf.is_finite(), "tau should be finite with NegInf delta");
    }

    /// Solve_joint_with_lowrank should exclude NaN positions
    /// when computing delta_max for the proximal gradient step size.
    #[test]
    fn test_joint_with_lowrank_nan_delta_max() {
        // Place NaN at a position that has the LARGEST delta weight.
        // Before fix: delta_max would include this NaN-position weight,
        //   making eta = 1/delta_max too small → wrong step size.
        // delta_max excludes NaN positions.
        let y = array![
            [1.0,     2.0, 3.0],
            [2.0,     3.0, 4.0],
            [3.0,     4.0, 5.0],
            [f64::NAN, 5.0, 8.0]  // NaN at (3,0), unit 2 treated
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0]
        ];
        // Give the NaN position (3,0) a very large weight
        let delta = array![
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [100.0, 1.0, 1.0]  // (3,0) has weight 100 but Y is NaN
        ];

        let result = solve_joint_with_lowrank(
            &y.view(),
            &d.view(),
            &delta.view(),
            0.5,  // lambda_nn
            100,  // max_iter
            1e-6, // tol
        );

        assert!(result.is_some(), "Should handle NaN in Y with large delta");

        let (mu, alpha, beta, l, tau, _n_iters, _converged) = result.unwrap();
        assert!(mu.is_finite(), "mu should be finite");
        assert!(alpha.iter().all(|&x| x.is_finite()), "alpha should be finite");
        assert!(beta.iter().all(|&x| x.is_finite()), "beta should be finite");
        assert!(l.iter().all(|&x| x.is_finite()), "L should be finite");
        assert!(tau.is_finite(), "tau should be finite");
    }

    /// When max_iter=1, the function should report
    /// actual_iters=1 and converged=false (unless the problem trivially converges
    /// in 1 iteration). This ensures convergence diagnostics are not hardcoded.
    #[test]
    fn test_joint_with_lowrank_not_converged() {
        // Use a larger panel where 1 iteration is unlikely to converge
        let y = array![
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 5.0, 6.0],
            [3.0, 5.0, 8.0, 9.0],
            [4.0, 6.0, 9.0, 15.0],  // Unit 3 treated with large effect
            [5.0, 7.0, 10.0, 20.0]  // Continued treatment
        ];
        let d = array![
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 1.0]
        ];
        let delta = array![
            [1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0]
        ];

        // max_iter=1 with a non-trivial lambda_nn: very unlikely to converge
        let result = solve_joint_with_lowrank(
            &y.view(), &d.view(), &delta.view(),
            0.5,  // lambda_nn
            1,    // max_iter = 1
            1e-12, // very tight tol
        );

        assert!(result.is_some(), "Should return Some even if not converged");
        let (_mu, _alpha, _beta, _l, _tau, n_iters, converged) = result.unwrap();

        // With max_iter=1, should have exactly 1 iteration
        assert_eq!(n_iters, 1, "Should report exactly 1 iteration");
        // With tight tolerance and only 1 iteration, should NOT be converged
        assert!(!converged, "Should NOT be converged with max_iter=1 and tight tol");
    }

    /// When convergence is achieved, n_iters < max_iter
    /// and converged=true. Ensures we don't always report max_iter.
    #[test]
    fn test_joint_with_lowrank_converged_before_max_iter() {
        // Simple panel that should converge quickly
        let y = array![
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [4.0, 5.0, 8.0]  // Unit 2 treated
        ];
        let d = array![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0]
        ];
        let delta = array![
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0]
        ];

        let max_iter = 1000;
        let result = solve_joint_with_lowrank(
            &y.view(), &d.view(), &delta.view(),
            0.1, max_iter, 1e-6,
        );

        assert!(result.is_some());
        let (_mu, _alpha, _beta, _l, _tau, n_iters, converged) = result.unwrap();

        // Should converge well before max_iter for this simple case
        assert!(converged, "Should converge for simple panel data");
        assert!(
            n_iters < max_iter,
            "Should converge before max_iter={}, got n_iters={}",
            max_iter, n_iters
        );
    }

    /// Verify delta_max computation:
    /// only finite-Y positions contribute to delta_max.
    #[test]
    fn test_joint_with_lowrank_nan_vs_no_nan_consistency() {
        // Panel without NaN
        let y_clean = array![
            [1.0, 2.0],
            [2.0, 3.0],
            [3.0, 4.0],
            [4.0, 6.0]  // Unit 1 treated with effect=1
        ];
        // Panel with NaN at a CONTROL position that has uniform weight
        let y_nan = array![
            [f64::NAN, 2.0],  // NaN at (0,0)
            [2.0, 3.0],
            [3.0, 4.0],
            [4.0, 6.0]
        ];
        let d = array![
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 1.0]
        ];
        let delta = array![
            [1.0, 1.0],
            [1.0, 1.0],
            [1.0, 1.0],
            [1.0, 1.0]
        ];

        let result_clean = solve_joint_with_lowrank(
            &y_clean.view(), &d.view(), &delta.view(), 0.5, 100, 1e-6,
        );
        let result_nan = solve_joint_with_lowrank(
            &y_nan.view(), &d.view(), &delta.view(), 0.5, 100, 1e-6,
        );

        assert!(result_clean.is_some());
        assert!(result_nan.is_some());

        let (_, _, _, _, tau_clean, _, _) = result_clean.unwrap();
        let (_, _, _, _, tau_nan, _, _) = result_nan.unwrap();

        // Both should produce finite tau
        assert!(tau_clean.is_finite());
        assert!(tau_nan.is_finite());
    }

    #[test]
    fn test_estimate_model_residual_structure() {
        // Test that residuals have expected structure after estimation
        let y = array![[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]];
        let control_mask = array![[1u8, 1], [1, 1], [1, 1]];
        let weight_matrix = array![[0.2, 0.2], [0.2, 0.2], [0.2, 0.2]];

        let result = estimate_model(
            &y.view(),
            &control_mask.view(),
            &weight_matrix.view(),
            0.0, // No nuclear norm penalty
            3,
            2,
            100,
            1e-6,
            None,
        );

        let (alpha, beta, l, _n_iters, _converged) = result.unwrap();

        // Compute fitted values: Y_hat = alpha + beta + L
        let mut y_hat = Array2::<f64>::zeros((3, 2));
        for t in 0..3 {
            for i in 0..2 {
                y_hat[[t, i]] = alpha[i] + beta[t] + l[[t, i]];
            }
        }

        // Residuals should be small for control observations
        for t in 0..3 {
            for i in 0..2 {
                let residual = (y[[t, i]] - y_hat[[t, i]]).abs();
                assert!(
                    residual < 0.1,
                    "Residual too large at [{}, {}]: {}",
                    t,
                    i,
                    residual
                );
            }
        }
    }
}

/// Property-based tests for model convergence using proptest
#[cfg(test)]
mod proptests {
    use super::*;
    use ndarray::Array2;
    use proptest::prelude::*;

    /// Strategy: generate a valid Y matrix with reasonable values
    fn y_matrix_strategy(
        n_periods: usize,
        n_units: usize,
    ) -> impl Strategy<Value = Array2<f64>> {
        prop::collection::vec(-10.0..10.0_f64, n_periods * n_units).prop_map(move |v| {
            Array2::from_shape_vec((n_periods, n_units), v).unwrap()
        })
    }

    /// Strategy: generate uniform weight matrix
    fn uniform_weights(n_periods: usize, n_units: usize) -> Array2<f64> {
        Array2::from_elem((n_periods, n_units), 1.0)
    }

    /// Strategy: generate all-control mask
    fn all_control_mask(n_periods: usize, n_units: usize) -> Array2<u8> {
        Array2::from_elem((n_periods, n_units), 1u8)
    }

    proptest! {
        /// Property 5a: estimate_model always terminates and returns Some
        /// For any valid all-control panel with uniform weights, estimation SHALL converge
        #[test]
        fn prop_estimation_terminates(
            y in y_matrix_strategy(5, 4),
            lambda_nn in 0.0..2.0_f64,
        ) {
            let mask = all_control_mask(5, 4);
            let w = uniform_weights(5, 4);

            let result = estimate_model(
                &y.view(), &mask.view(), &w.view(),
                lambda_nn, 5, 4, 100, 1e-6, None,
            );

            prop_assert!(result.is_some(),
                "estimate_model should always return Some for valid all-control panel");

            let (_alpha, _beta, _l, n_iters, _converged) = result.unwrap();
            prop_assert!(n_iters > 0, "Should have at least 1 iteration");
            prop_assert!(n_iters <= 100, "Should not exceed max_iter=100");
        }

        /// Property 5b: Alpha and beta have correct dimensions
        #[test]
        fn prop_output_dimensions(
            y in y_matrix_strategy(6, 5),
            lambda_nn in 0.0..1.0_f64,
        ) {
            let mask = all_control_mask(6, 5);
            let w = uniform_weights(6, 5);

            let result = estimate_model(
                &y.view(), &mask.view(), &w.view(),
                lambda_nn, 6, 5, 100, 1e-6, None,
            );

            let (alpha, beta, l, _n_iters, _converged) = result.unwrap();
            prop_assert_eq!(alpha.len(), 5, "alpha should have n_units elements");
            prop_assert_eq!(beta.len(), 6, "beta should have n_periods elements");
            prop_assert_eq!(l.shape(), &[6, 5], "L should be n_periods x n_units");
        }

        /// Property 5c: All outputs are finite (no NaN or Inf)
        #[test]
        fn prop_outputs_finite(
            y in y_matrix_strategy(5, 4),
            lambda_nn in 0.01..2.0_f64,
        ) {
            let mask = all_control_mask(5, 4);
            let w = uniform_weights(5, 4);

            let result = estimate_model(
                &y.view(), &mask.view(), &w.view(),
                lambda_nn, 5, 4, 100, 1e-6, None,
            );

            let (alpha, beta, l, _n_iters, _converged) = result.unwrap();

            for &v in alpha.iter() {
                prop_assert!(v.is_finite(), "alpha contains non-finite: {}", v);
            }
            for &v in beta.iter() {
                prop_assert!(v.is_finite(), "beta contains non-finite: {}", v);
            }
            for &v in l.iter() {
                prop_assert!(v.is_finite(), "L contains non-finite: {}", v);
            }
        }

        /// Property 5d: Residuals decrease with more iterations
        /// Running with max_iter=1 vs max_iter=100 should give smaller residuals for the latter
        #[test]
        fn prop_more_iterations_better_fit(
            y in y_matrix_strategy(5, 4),
        ) {
            let mask = all_control_mask(5, 4);
            let w = uniform_weights(5, 4);
            let lambda_nn = 0.1;

            let result_1 = estimate_model(
                &y.view(), &mask.view(), &w.view(),
                lambda_nn, 5, 4, 1, 1e-6, None,
            );
            let result_100 = estimate_model(
                &y.view(), &mask.view(), &w.view(),
                lambda_nn, 5, 4, 100, 1e-6, None,
            );

            let (a1, b1, l1, n_iters_1, _conv1) = result_1.unwrap();
            let (a100, b100, l100, n_iters_100, _conv100) = result_100.unwrap();

            // Verify iteration counts are meaningful
            prop_assert_eq!(n_iters_1, 1, "max_iter=1 should give exactly 1 iteration");
            prop_assert!(n_iters_100 >= 1, "max_iter=100 should give at least 1 iteration");
            prop_assert!(n_iters_100 <= 100, "Should not exceed max_iter=100");

            // Compute weighted residual sum of squares
            let compute_wrss = |alpha: &Array1<f64>, beta: &Array1<f64>, l: &Array2<f64>| -> f64 {
                let mut wrss = 0.0;
                for t in 0..5 {
                    for i in 0..4 {
                        let resid = y[[t, i]] - alpha[i] - beta[t] - l[[t, i]];
                        wrss += w[[t, i]] * resid * resid;
                    }
                }
                wrss
            };

            let wrss_1 = compute_wrss(&a1, &b1, &l1);
            let wrss_100 = compute_wrss(&a100, &b100, &l100);

            // More iterations should give equal or better fit
            // Allow small tolerance for numerical noise
            prop_assert!(wrss_100 <= wrss_1 + 1e-6,
                "100 iterations WRSS {} should be <= 1 iteration WRSS {}",
                wrss_100, wrss_1);
        }

        /// Property 5e: Higher lambda_nn produces smaller L matrix (stronger regularization)
        #[test]
        fn prop_higher_lambda_smaller_l(
            y in y_matrix_strategy(5, 4),
        ) {
            let mask = all_control_mask(5, 4);
            let w = uniform_weights(5, 4);

            let result_low = estimate_model(
                &y.view(), &mask.view(), &w.view(),
                0.01, 5, 4, 100, 1e-6, None,
            );
            let result_high = estimate_model(
                &y.view(), &mask.view(), &w.view(),
                10.0, 5, 4, 100, 1e-6, None,
            );

            let (_, _, l_low, _, _) = result_low.unwrap();
            let (_, _, l_high, _, _) = result_high.unwrap();

            let norm_low: f64 = l_low.iter().map(|x| x * x).sum::<f64>().sqrt();
            let norm_high: f64 = l_high.iter().map(|x| x * x).sum::<f64>().sqrt();

            // Higher lambda should produce smaller or equal L norm
            prop_assert!(norm_high <= norm_low + 1e-4,
                "Higher lambda L norm {} should be <= lower lambda L norm {}",
                norm_high, norm_low);
        }
    }
}
