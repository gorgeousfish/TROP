/*==============================================================================
  trop_estimator_core.mata

  Simplified Mata implementation of the TROP per-observation estimator.
  Provides a self-contained reference path for the alternating minimization
  algorithm when the compiled plugin back-end is unavailable.

  Structures defined here:
    LambdaSearchResult  — return type for lambda search
    fe_result           — unit and time fixed-effect estimates
    ADMMResult          — nuclear-norm regularised low-rank estimate

  Exported entry point:
    trop_estimate_simple()  — per-observation model fitting and ATT estimation
==============================================================================*/

version 17
mata:
mata set matastrict on

/*==============================================================================
  Part 1: Auxiliary structures
==============================================================================*/

// Return value from the cyclic coordinate descent lambda search.
struct LambdaSearchResult {
    real rowvector lambda_opt   // optimal (lambda_time, lambda_unit, lambda_nn)
    real scalar cv_loss_opt     // cross-validation loss at optimum
    real scalar converged       // convergence indicator (0/1)
    real scalar iterations      // number of iterations performed
}

// Fixed-effect estimates from the alternating minimisation step.
struct fe_result {
    real colvector alpha        // unit fixed effects (N x 1)
    real colvector beta         // time fixed effects (T x 1)
}

// Result of nuclear-norm regularised matrix estimation.
struct ADMMResult {
    real matrix mu              // low-rank component L-hat (N x T)
    real scalar converged       // convergence indicator (0/1)
    real scalar iterations      // number of iterations performed
}

/*==============================================================================
  Part 2: Helper functions
==============================================================================*/

/*------------------------------------------------------------------------------
  _tec_lambda_search()

  Stub for the three-stage cyclic coordinate descent lambda search.
  Returns default values without performing actual LOOCV optimisation.
  The full search is handled by the compiled plugin back-end.
------------------------------------------------------------------------------*/

struct LambdaSearchResult scalar _tec_lambda_search(
    real matrix Y,
    real matrix W,
    string scalar strategy,
    string scalar cv_mode,
    real scalar maxiter
)
{
    struct LambdaSearchResult scalar result

    printf("{err}Warning: _tec_lambda_search() is a stub; returning defaults.\n")
    printf("{err}Use the compiled plugin back-end for production LOOCV.\n")

    result.lambda_opt = (1.0, 1.0, 0.1)
    result.cv_loss_opt = .
    result.converged = 0
    result.iterations = 0

    return(result)
}

/*------------------------------------------------------------------------------
  _tec_compute_time_weights()

  Compute time weights theta_s for a given target period t_target.

    theta_s = exp(-lambda_time * |s - t_target|)

  Weights are assigned only at control periods indicated by `m_control`.
  The returned vector is unnormalised.
------------------------------------------------------------------------------*/

real colvector _tec_compute_time_weights(
    real scalar lambda_time,
    real scalar T,
    real scalar t_target,
    real colvector m_control
)
{
    real colvector theta
    real scalar s

    theta = J(T, 1, 0)
    for (s = 1; s <= T; s++) {
        if (m_control[s] == 1) {
            theta[s] = exp(-lambda_time * abs(s - t_target))
        }
    }

    return(theta)
}

/*------------------------------------------------------------------------------
  _tec_compute_unit_distance()

  Compute the cross-unit distance between units i_target and j,
  excluding the target period t_target:

    dist_{-t}(j, i) = sqrt( sum_{u != t} (1-W_{iu})(1-W_{ju})(Y_{iu}-Y_{ju})^2
                            / sum_{u != t} (1-W_{iu})(1-W_{ju}) )

  Returns missing (.) when no jointly-untreated periods exist.
------------------------------------------------------------------------------*/

real scalar _tec_compute_unit_distance(
    real matrix Y,
    real matrix W,
    real scalar i_target,
    real scalar j,
    real scalar t_target
)
{
    real scalar T, s, n_valid, sum_sq

    T = cols(Y)
    n_valid = 0
    sum_sq = 0

    for (s = 1; s <= T; s++) {
        if (s == t_target) continue
        if (W[i_target, s] == 0 & W[j, s] == 0) {
            sum_sq = sum_sq + (Y[i_target, s] - Y[j, s])^2
            n_valid++
        }
    }

    if (n_valid > 0) {
        return(sqrt(sum_sq / n_valid))
    }
    else {
        return(.)
    }
}

/*------------------------------------------------------------------------------
  _tec_compute_unit_weights()

  Compute unit weights omega_j using an exponential kernel:

    omega_j = exp(-lambda_unit * dist_j)

  When lambda_unit == 0 the kernel reduces to uniform weights.
  Only units flagged in `control_mask` with finite distances receive
  positive weights.  The returned vector is unnormalised.
------------------------------------------------------------------------------*/

real colvector _tec_compute_unit_weights(
    real scalar lambda_unit,
    real colvector dist_row,
    real colvector control_mask
)
{
    real colvector l
    real scalar N, j

    N = rows(dist_row)
    l = J(N, 1, 0)

    for (j = 1; j <= N; j++) {
        if (control_mask[j] == 1 & dist_row[j] < .) {
            if (lambda_unit == 0) {
                l[j] = 1
            }
            else {
                l[j] = exp(-lambda_unit * dist_row[j])
            }
        }
    }

    return(l)
}

/*------------------------------------------------------------------------------
  _tec_estimate_fe()

  Estimate unit and time fixed effects by weighted least squares.

    alpha_i = sum_s w_{is} (Y_{is} - L_{is}) / sum_s w_{is}
    beta_s  = sum_j w_{js} (Y_{js} - alpha_j - L_{js}) / sum_j w_{js}
------------------------------------------------------------------------------*/

struct fe_result scalar _tec_estimate_fe(
    real matrix Y,
    real matrix W_weight,
    real matrix mu
)
{
    struct fe_result scalar res
    real scalar N, T, i, s
    real scalar numer, denom
    real matrix R

    N = rows(Y)
    T = cols(Y)

    R = Y - mu

    // Unit fixed effects: weighted row means of R.
    res.alpha = J(N, 1, 0)
    for (i = 1; i <= N; i++) {
        numer = 0
        denom = 0
        for (s = 1; s <= T; s++) {
            if (W_weight[i, s] > 0) {
                numer = numer + W_weight[i, s] * R[i, s]
                denom = denom + W_weight[i, s]
            }
        }
        if (denom > 0) {
            res.alpha[i] = numer / denom
        }
    }

    // Time fixed effects: weighted column means of (R - alpha).
    res.beta = J(T, 1, 0)
    for (s = 1; s <= T; s++) {
        numer = 0
        denom = 0
        for (i = 1; i <= N; i++) {
            if (W_weight[i, s] > 0) {
                numer = numer + W_weight[i, s] * (R[i, s] - res.alpha[i])
                denom = denom + W_weight[i, s]
            }
        }
        if (denom > 0) {
            res.beta[s] = numer / denom
        }
    }

    return(res)
}

/*------------------------------------------------------------------------------
  _tec_nuclear_norm()

  Minimise the nuclear-norm penalised weighted least squares objective:

    min_L  sum_{j,s} w_{js} (R_{js} - L_{js})^2  +  lambda_nn * ||L||_*

  Uses a single-step SVD soft-thresholding approximation (proximal operator).

    sigma_k_hat = max(0, sigma_k - lambda_nn)
    L_hat       = U * diag(sigma_hat) * V'
------------------------------------------------------------------------------*/

struct ADMMResult scalar _tec_nuclear_norm(
    real matrix R,
    real matrix W_weight,
    real scalar lambda_nn
)
{
    struct ADMMResult scalar res
    real matrix U, Vt, S_mat
    real colvector s_vals, s_thresh
    real scalar N, T, rank, k

    N = rows(R)
    T = cols(R)

    // Mask the residual matrix to observed (weighted) entries.
    real matrix R_weighted
    R_weighted = R :* (W_weight :> 0)

    // SVD soft-thresholding.
    _svd(R_weighted, s_vals=., U=., Vt=.)

    rank = min((N, T))
    if (length(s_vals) < rank) {
        rank = length(s_vals)
    }

    s_thresh = J(rank, 1, 0)
    for (k = 1; k <= rank; k++) {
        s_thresh[k] = max((0, s_vals[k] - lambda_nn))
    }

    // Reconstruct L_hat = U * diag(sigma_hat) * V'.
    S_mat = diag(s_thresh)
    if (rows(U) > 0 & cols(Vt) > 0) {
        res.mu = U[., 1..rank] * S_mat * Vt[1..rank, .]
    }
    else {
        res.mu = J(N, T, 0)
    }

    res.converged = 1
    res.iterations = 1

    return(res)
}

/*==============================================================================
  Part 3: Main estimation function

  trop_estimate_simple()

  Per-observation model fitting (Algorithm 2).  For each treated
  observation (i, t) an independent set of kernel weights and a separate
  additive model  Y = alpha + beta + L  are estimated.  Individual
  treatment effects  tau_{it} = Y_{it} - alpha_i - beta_t - L_{it}
  are averaged to produce the ATT.
==============================================================================*/

void trop_estimate_simple(
    string scalar depvar,
    string scalar treatvar,
    string scalar panel_idx_var,
    string scalar time_idx_var,
    string scalar touse_var,
    real scalar lambda_time_user,
    real scalar lambda_unit_user,
    real scalar lambda_nn_user,
    real scalar need_cv,
    string scalar cv_mode,
    string scalar strategy,
    real scalar maxiter,
    real scalar tau_true
)
{
    real matrix Y, W, W_weight, mu, mu_old, R, Y_hat_0
    real matrix treated_pairs, alpha_sum, beta_sum, mu_sum
    real colvector Y_vec, W_vec, panel_idx, time_idx
    real colvector alpha, beta, theta, l, dist_row, tau_it_vec
    real vector m_control, control_mask_t
    real rowvector lambda_opt
    real scalar N, T, i, j, idx
    real scalar tau, se, rmse_control, bias
    real scalar n_treat, n_control, sum_sq_err
    real scalar var_tau, t_stat, pvalue, ci_lower, ci_upper
    real scalar alpha_level, df_pvalue
    real scalar t_target, i_target
    real scalar iter, converged, max_outer_iter, delta_mu
    real scalar cv_loss_final, converged_final, iterations_final
    real scalar n_treated_pairs, pair_idx, n_success, total_converged
    struct LambdaSearchResult scalar search_result
    struct fe_result scalar fe_res
    struct ADMMResult scalar admm_res

    // Initialise weight vectors to empty.
    theta = J(0, 1, .)
    l = J(0, 1, .)

    // ---- Load data from Stata -----------------------------------------------
    Y_vec = st_data(., depvar, touse_var)
    W_vec = st_data(., treatvar, touse_var)
    panel_idx = st_data(., panel_idx_var, touse_var)
    time_idx = st_data(., time_idx_var, touse_var)

    N = max(panel_idx)
    T = max(time_idx)

    // Reshape long-format vectors into N x T panel matrices.
    Y = J(N, T, .)
    W = J(N, T, .)
    for (idx = 1; idx <= rows(Y_vec); idx++) {
        i = panel_idx[idx]
        j = time_idx[idx]
        Y[i, j] = Y_vec[idx]
        W[i, j] = W_vec[idx]
    }

    n_treat = sum(W)
    n_control = N * T - n_treat

    // ---- Lambda selection ---------------------------------------------------
    if (need_cv) {
        printf("\n========== Lambda search (cyclic coordinate descent) ==========\n")

        search_result = _tec_lambda_search(Y, W, strategy, cv_mode, maxiter)

        lambda_opt = search_result.lambda_opt
        cv_loss_final = search_result.cv_loss_opt
        converged_final = search_result.converged
        iterations_final = search_result.iterations
    }
    else {
        lambda_opt = (lambda_time_user, lambda_unit_user, lambda_nn_user)
        cv_loss_final = .
        converged_final = 1
        iterations_final = 0
    }

    // ---- Per-observation estimation -----------------------------------------
    printf("\n========== Estimation with selected lambda ==========\n")
    printf("  lambda = (%g, %g, %g)\n", lambda_opt[1], lambda_opt[2], lambda_opt[3])

    // Enumerate all treated observations (i, t).
    treated_pairs = J(0, 2, .)
    for (i=1; i<=N; i++) {
        for (j=1; j<=T; j++) {
            if (W[i,j] == 1) {
                treated_pairs = treated_pairs \ (i, j)
            }
        }
    }
    n_treated_pairs = rows(treated_pairs)

    // Default initialisations.
    alpha = J(N, 1, 0)
    beta = J(T, 1, 0)
    mu = J(N, T, 0)
    t_target = .
    i_target = .
    n_success = 0
    total_converged = 1

    if (n_treated_pairs == 0) {
        printf("  Error: no treated observations found\n")
        tau = 0
        se = 0
        tau_it_vec = J(0, 1, .)
    }
    else {
        printf("  Treated observations: %g\n", n_treated_pairs)

        tau_it_vec = J(n_treated_pairs, 1, .)
        alpha_sum = J(N, 1, 0)
        beta_sum = J(T, 1, 0)
        mu_sum = J(N, T, 0)
        max_outer_iter = 10

        for (pair_idx=1; pair_idx<=n_treated_pairs; pair_idx++) {
            i_target = treated_pairs[pair_idx, 1]
            t_target = treated_pairs[pair_idx, 2]

            // Time weights for (i_target, t_target).
            m_control = J(T, 1, 0)
            for (j = 1; j <= T; j++) {
                if (W[i_target, j] == 0) m_control[j] = 1
            }
            theta = _tec_compute_time_weights(lambda_opt[1], T, t_target, m_control)

            // Unit distances and weights for (i_target, t_target).
            control_mask_t = J(N, 1, 0)
            for (i = 1; i <= N; i++) {
                if (W[i, t_target] == 0) control_mask_t[i] = 1
            }
            dist_row = J(N, 1, .)
            for (j = 1; j <= N; j++) {
                if (j != i_target) {
                    dist_row[j] = _tec_compute_unit_distance(Y, W, i_target, j, t_target)
                }
                else {
                    dist_row[j] = 0
                }
            }
            l = _tec_compute_unit_weights(lambda_opt[2], dist_row, control_mask_t)

            // Combined weight matrix:  w_{js} = omega_j * theta_s * (1 - W_{js}).
            W_weight = (l * theta') :* (1 :- W)

            // Alternating minimisation for this (i, t).
            alpha = J(N, 1, 0)
            beta = J(T, 1, 0)
            mu = J(N, T, 0)
            converged = 0

            for (iter=1; iter<=max_outer_iter; iter++) {
                mu_old = mu

                fe_res = _tec_estimate_fe(Y, W_weight, mu)
                alpha = fe_res.alpha
                beta = fe_res.beta

                if (lambda_opt[3] >= 1e9) {
                    mu = J(N, T, 0)
                }
                else {
                    R = Y :- alpha * J(1,T,1) :- J(N,1,1) * beta'
                    admm_res = _tec_nuclear_norm(R, W_weight, lambda_opt[3])
                    mu = admm_res.mu
                }

                delta_mu = sqrt(sum((mu :- mu_old) :^ 2)) / (sqrt(sum(mu :^ 2)) + 1e-10)
                if (delta_mu < 1e-4) {
                    converged = 1
                    break
                }
            }

            if (!converged) total_converged = 0

            // Individual treatment effect: tau_{it} = Y_{it} - alpha_i - beta_t - L_{it}.
            tau_it_vec[pair_idx] = Y[i_target, t_target] - alpha[i_target] - beta[t_target] - mu[i_target, t_target]

            // Accumulate parameter estimates for subsequent averaging.
            alpha_sum = alpha_sum + alpha
            beta_sum = beta_sum + beta
            mu_sum = mu_sum + mu
            n_success = n_success + 1
        }

        printf("  Per-observation fitting complete: %g/%g succeeded", n_success, n_treated_pairs)
        if (total_converged) {
            printf(" (all converged)\n")
        } else {
            printf(" (some did not converge)\n")
        }

        // Average parameter estimates across treated observations.
        if (n_success > 0) {
            alpha = alpha_sum / n_success
            beta = beta_sum / n_success
            mu = mu_sum / n_success
        }

        // ATT and standard error.
        if (n_success > 0) {
            tau = mean(tau_it_vec)
            if (n_success > 1) {
                var_tau = variance(tau_it_vec)
                se = sqrt(var_tau / n_success)
            }
            else {
                se = 0
            }
        }
        else {
            tau = 0
            se = 0
        }
    }

    // ---- Counterfactual predictions and RMSE --------------------------------
    Y_hat_0 = J(N, T, 0)
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= T; j++) {
            Y_hat_0[i, j] = alpha[i] + beta[j] + mu[i, j]
        }
    }

    sum_sq_err = 0
    n_control = 0
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= T; j++) {
            if (W[i, j] == 0) {
                sum_sq_err = sum_sq_err + (Y[i, j] - Y_hat_0[i, j])^2
                n_control++
            }
        }
    }
    if (n_control > 0) {
        rmse_control = sqrt(sum_sq_err / n_control)
    } else {
        rmse_control = .
    }

    // Bias relative to a known true effect (when available).
    if (tau_true != -999) {
        bias = tau - tau_true
    } else {
        bias = .
    }

    // ---- Inference -----------------------------------------------------------
    // Confidence level is read from c(level); default to 95% if invalid.
    alpha_level = 1 - st_numscalar("c(level)") / 100
    if (alpha_level <= 0 | alpha_level >= 1) alpha_level = 0.05

    if (se > 0 & se < .) {
        t_stat = tau / se
        if (n_treated_pairs >= 1) {
            df_pvalue = max((1, n_treated_pairs - 1))
            pvalue = 2 * ttail(df_pvalue, abs(t_stat))
            ci_lower = tau - invttail(df_pvalue, alpha_level / 2) * se
            ci_upper = tau + invttail(df_pvalue, alpha_level / 2) * se
        }
        else {
            pvalue = 2 * normal(-abs(t_stat))
            ci_lower = tau - invnormal(1 - alpha_level / 2) * se
            ci_upper = tau + invnormal(1 - alpha_level / 2) * se
        }
    }
    else {
        t_stat = .
        pvalue = .
        ci_lower = .
        ci_upper = .
    }

    // ---- Store results in e() -----------------------------------------------
    st_numscalar("e(N)", N)
    st_numscalar("e(T)", T)
    st_numscalar("e(N_obs)", N * T)
    st_numscalar("e(N_treat)", n_treat)
    st_numscalar("e(N_control)", n_control)

    st_numscalar("e(lambda_time)", lambda_opt[1])
    st_numscalar("e(lambda_unit)", lambda_opt[2])
    st_numscalar("e(lambda_nn)", lambda_opt[3])

    st_numscalar("e(tau)", tau)
    st_numscalar("e(se)", se)
    st_numscalar("e(t)", t_stat)
    st_numscalar("e(pvalue)", pvalue)
    st_numscalar("e(ci_lower)", ci_lower)
    st_numscalar("e(ci_upper)", ci_upper)

    st_numscalar("e(alpha_level)", alpha_level)
    st_numscalar("e(N_treated_pairs)", n_treated_pairs)
    st_numscalar("e(rmse_control)", rmse_control)
    st_numscalar("e(bias)", bias)

    if (need_cv) {
        st_numscalar("e(cv_loss)", cv_loss_final)
        st_numscalar("e(cv_converged)", converged_final)
        st_numscalar("e(cv_iterations)", iterations_final)
    } else {
        st_numscalar("e(cv_loss)", .)
        st_numscalar("e(cv_converged)", 1)
        st_numscalar("e(cv_iterations)", 0)
    }

    st_matrix("e(alpha)", alpha)
    st_matrix("e(beta)", beta)
    st_matrix("e(mu)", mu)
    if (t_target < .) {
        // Weights from the last treated observation (representative snapshot).
        st_matrix("e(theta)", theta)
        st_matrix("e(l)", l)
    }

    st_numscalar("e(converged)", total_converged)
    st_numscalar("e(iterations)", n_success)

    st_global("e(cmd)", "trop_estimator")
    st_global("e(depvar)", depvar)
    st_global("e(treatvar)", treatvar)
    st_global("e(cv_mode)", cv_mode)
    st_global("e(strategy)", strategy)
}

end
