/*──────────────────────────────────────────────────────────────────────────────
  trop_ereturn_store.mata

  Post-estimation result storage for the TROP estimator.

  Reads intermediate scalars and matrices written by the computational
  backend (temporary Stata scalars/matrices with the __trop_ prefix) and
  posts them into Stata's e() return structure.

  The TROP model decomposes the outcome as

      Y_{it} = alpha_i + beta_t + mu_{it} + tau_{it} * W_{it}

  where alpha_i are unit fixed effects, beta_t are time fixed effects,
  mu_{it} is a low-rank (nuclear-norm penalised) component, and tau_{it}
  is the treatment effect.

  e(mu) semantics differ by estimation method:
    - Twostep:  e(mu) = .  (no global intercept; model is Y = alpha + beta + L)
    - Joint:    e(mu) = mu (scalar intercept with identification constraint
                alpha_1 = beta_1 = 0)

  Contents
    trop_store_results()                main e() storage entry point
    _trop_safe_read_scalar()            safe scalar reader
    _trop_display_warnings()            estimation diagnostic messages
    _trop_compute_effective_rank()      SVD-based effective rank
    _trop_display_bootstrap_warnings()  bootstrap diagnostic messages
    _trop_store_lambda_grids()          LOOCV grid and CV curve storage
──────────────────────────────────────────────────────────────────────────────*/

version 17
mata:
mata set matastrict on

/*──────────────────────────────────────────────────────────────────────────────
  trop_store_results()

  Main entry point for e() result storage.

  Arguments
    method   "twostep" or "joint"

  Stored e() returns
    Scalars   att, se, t, ci_lower, ci_upper, pvalue, df_r, mu,
              lambda_time, lambda_unit, lambda_nn,
              loocv_first_failed_t, loocv_first_failed_i,
              n_iterations, converged, n_obs_estimated, n_obs_failed,
              n_bootstrap_valid, level, ci_lower_percentile,
              ci_upper_percentile, N_treated, N_treated_obs,
              N_treated_units, effective_rank, alpha_level
    Matrices  alpha, beta, tau, factor_matrix, bootstrap_estimates,
              theta/omega (twostep), delta_time/delta_unit (joint),
              lambda_grid, cv_curve
    Globals   title, predict
──────────────────────────────────────────────────────────────────────────────*/

void trop_store_results(string scalar method)
{
    real scalar att, se, ci_lower, ci_upper, pvalue
    real scalar lambda_time, lambda_unit, lambda_nn
    real scalar n_iterations, converged, bootstrap_reps, alpha_level
    real scalar n_units, n_periods, n_treated, effective_rank
    real scalar loocv_first_failed_t, loocv_first_failed_i
    real scalar loocv_score, loocv_n_valid, loocv_n_attempted
    real matrix alpha, beta, factor_matrix, tau, bootstrap_estimates
    real scalar mu
    real scalar n_bootstrap_valid, level
    real scalar tstat
    real scalar n_treated_obs, n_treated_units, df_pvalue
    real scalar ci_lower_pct, ci_upper_pct
    real scalar n_obs_estimated, n_obs_failed
    real matrix temp_mat, temp_scalar

    /* ── verify that the core ATT scalar exists ──────────────────────── */
    temp_scalar = st_numscalar("__trop_att")
    if (rows(temp_scalar) == 0) {
        errprintf("Error: TROP estimation did not produce results.\n")
        errprintf("       __trop_att is missing. Check if plugin executed correctly.\n")
        _error(3200)
    }
    att = temp_scalar

    /* ── core estimation results ─────────────────────────────────────── */

    se = _trop_safe_read_scalar("__trop_se")

    /* Resolve significance level: prefer bootstrap-specific alpha,
       fall back to estimation-level alpha, then default 0.05. */
    alpha_level = _trop_safe_read_scalar("__trop_bs_alpha")
    if (alpha_level >= . | alpha_level <= 0 | alpha_level >= 1) {
        alpha_level = _trop_safe_read_scalar("__trop_alpha_level")
    }
    if (alpha_level >= . | alpha_level <= 0 | alpha_level >= 1) {
        alpha_level = 0.05
    }

    /* Number of treated observations (needed for degrees of freedom).
       Twostep: length of the tau vector (one per treated cell).
       Joint:   count of W_{it}=1 cells in the panel. */
    n_treated_obs = _trop_safe_read_scalar("__trop_n_treated")
    if (n_treated_obs >= .) {
        temp_mat = st_matrix("__trop_tau")
        if (rows(temp_mat) > 0) {
            n_treated_obs = rows(temp_mat)
        }
        else {
            n_treated_obs = .
        }
    }

    /* CI and p-value use the same reference distribution:
         - t-distribution with df = max(1, n_treated - 1) when df is known
         - standard normal as asymptotic fallback */
    if (se > 0 && se < .) {
        tstat = att / se
        if (n_treated_obs < . && n_treated_obs >= 1) {
            df_pvalue = max((1, n_treated_obs - 1))
            pvalue = 2 * ttail(df_pvalue, abs(tstat))
            ci_lower = att - invttail(df_pvalue, alpha_level / 2) * se
            ci_upper = att + invttail(df_pvalue, alpha_level / 2) * se
        }
        else {
            pvalue = 2 * normal(-abs(tstat))
            ci_lower = att - invnormal(1 - alpha_level / 2) * se
            ci_upper = att + invnormal(1 - alpha_level / 2) * se
        }
    }
    else {
        tstat = .
        pvalue = .
        ci_lower = .
        ci_upper = .
    }

    /* ── regularization hyperparameters ──────────────────────────────── */
    lambda_time = _trop_safe_read_scalar("__trop_lambda_time")
    lambda_unit = _trop_safe_read_scalar("__trop_lambda_unit")
    lambda_nn = _trop_safe_read_scalar("__trop_lambda_nn")

    /* ── LOOCV diagnostics ───────────────────────────────────────────── */
    loocv_score = _trop_safe_read_scalar("__trop_loocv_score")
    loocv_n_valid = _trop_safe_read_scalar("__trop_loocv_n_valid")
    loocv_n_attempted = _trop_safe_read_scalar("__trop_loocv_n_attempted")
    loocv_first_failed_t = _trop_safe_read_scalar("__trop_loocv_first_failed_t")
    loocv_first_failed_i = _trop_safe_read_scalar("__trop_loocv_first_failed_i")

    /* ── convergence information ─────────────────────────────────────── */
    n_iterations = _trop_safe_read_scalar("__trop_n_iterations")
    converged = _trop_safe_read_scalar("__trop_converged")
    n_obs_estimated = _trop_safe_read_scalar("__trop_n_obs_estimated")
    n_obs_failed = _trop_safe_read_scalar("__trop_n_obs_failed")

    /* ── bootstrap configuration ─────────────────────────────────────── */
    bootstrap_reps = _trop_safe_read_scalar("__trop_n_bootstrap")

    /* ── bootstrap diagnostics ───────────────────────────────────────── */
    n_bootstrap_valid = _trop_safe_read_scalar("__trop_n_bootstrap_valid")
    level = _trop_safe_read_scalar("__trop_level")

    /* Percentile CI from the bootstrap distribution (diagnostic only).
       The authoritative CI uses the t-distribution (or normal fallback).
       Large discrepancy between percentile and parametric CI may indicate
       asymmetry in the bootstrap distribution. */
    ci_lower_pct = _trop_safe_read_scalar("__trop_ci_lower_percentile")
    ci_upper_pct = _trop_safe_read_scalar("__trop_ci_upper_percentile")

    /* ── sample information ──────────────────────────────────────────── */
    n_units = _trop_safe_read_scalar("__trop_n_units")
    n_periods = _trop_safe_read_scalar("__trop_n_periods")
    effective_rank = .

    /* ── read matrices from temporary storage ────────────────────────── */

    temp_mat = st_matrix("__trop_alpha")
    if (rows(temp_mat) > 0) alpha = temp_mat
    else alpha = J(0, 1, .)

    temp_mat = st_matrix("__trop_beta")
    if (rows(temp_mat) > 0) beta = temp_mat
    else beta = J(0, 1, .)

    temp_mat = st_matrix("__trop_factor_matrix")
    if (rows(temp_mat) > 0) factor_matrix = temp_mat
    else factor_matrix = J(0, 0, .)

    if (bootstrap_reps > 0) {
        temp_mat = st_matrix("__trop_bootstrap_estimates")
        if (rows(temp_mat) > 0) {
            /* Filter trailing missing values: the matrix is pre-allocated
               at full size; only the first n_valid rows are populated. */
            bootstrap_estimates = select(temp_mat, temp_mat :< .)
        }
        else {
            bootstrap_estimates = J(0, 1, .)
        }
    }
    else {
        bootstrap_estimates = J(0, 1, .)
    }

    /* ═══════════════════ store into e() ══════════════════════════════ */

    /* ── core estimation results ─────────────────────────────────────── */
    st_numscalar("e(att)", att)
    st_numscalar("e(se)", se)
    st_numscalar("e(t)", tstat)
    st_numscalar("e(ci_lower)", ci_lower)
    st_numscalar("e(ci_upper)", ci_upper)
    st_numscalar("e(pvalue)", pvalue)
    if (n_treated_obs < . && n_treated_obs >= 1) {
        st_numscalar("e(df_r)", max((1, n_treated_obs - 1)))
    }
    else {
        st_numscalar("e(df_r)", .)
    }

    /* ── regularization hyperparameters ──────────────────────────────── */
    st_numscalar("e(lambda_time)", lambda_time)
    st_numscalar("e(lambda_unit)", lambda_unit)
    st_numscalar("e(lambda_nn)", lambda_nn)

    /* ── LOOCV diagnostics (first-failure indices) ───────────────────── */
    st_numscalar("e(loocv_first_failed_t)", loocv_first_failed_t)
    st_numscalar("e(loocv_first_failed_i)", loocv_first_failed_i)

    /* ── convergence information ─────────────────────────────────────── */
    st_numscalar("e(n_iterations)", n_iterations)
    st_numscalar("e(converged)", converged)
    if (n_obs_estimated < .) {
        st_numscalar("e(n_obs_estimated)", n_obs_estimated)
    }
    if (n_obs_failed < . && n_obs_failed > 0) {
        st_numscalar("e(n_obs_failed)", n_obs_failed)
    }

    /* ── bootstrap diagnostics ───────────────────────────────────────── */
    st_numscalar("e(n_bootstrap_valid)", n_bootstrap_valid)
    st_numscalar("e(level)", level)
    if (ci_lower_pct < . & ci_upper_pct < .) {
        st_numscalar("e(ci_lower_percentile)", ci_lower_pct)
        st_numscalar("e(ci_upper_percentile)", ci_upper_pct)
    }

    /* ── parameter matrices ──────────────────────────────────────────── */
    st_matrix("e(alpha)", alpha)
    st_matrix("e(beta)", beta)
    st_matrix("e(factor_matrix)", factor_matrix)

    /* Effective rank via continuous SVD criterion: sum(s) / s[1] */
    effective_rank = _trop_compute_effective_rank(factor_matrix)
    if (effective_rank < .) {
        st_numscalar("e(effective_rank)", effective_rank)
    }

    /* ── bootstrap replicate distribution ────────────────────────────── */
    if (bootstrap_reps > 0 && rows(bootstrap_estimates) > 0) {
        st_matrix("e(bootstrap_estimates)", bootstrap_estimates)
    }

    /* ── lambda grids and CV curve ───────────────────────────────────── */
    _trop_store_lambda_grids(lambda_time, lambda_unit, lambda_nn, loocv_score)

    /* ── method-specific storage ─────────────────────────────────────── */
    n_treated_units = _trop_safe_read_scalar("__trop_n_treated_units")
    if (n_treated_units < .) {
        st_numscalar("e(N_treated_units)", n_treated_units)
    }

    if (method == "twostep") {
        /* Twostep: no global intercept; model is Y = alpha + beta + L */
        st_numscalar("e(mu)", .)

        /* Individual treatment effect vector tau_{i,t} */
        temp_mat = st_matrix("__trop_tau")
        if (rows(temp_mat) > 0) {
            tau = temp_mat
            st_matrix("e(tau)", tau)
            n_treated = rows(tau)
            st_numscalar("e(N_treated)", n_treated)
            st_numscalar("e(N_treated_obs)", n_treated)
        }

        /* Time weights theta (T x 1) and unit weights omega (N x 1) */
        temp_mat = st_matrix("__trop_theta")
        if (rows(temp_mat) > 0 && any(temp_mat)) {
            st_matrix("e(theta)", temp_mat)
        }
        temp_mat = st_matrix("__trop_omega")
        if (rows(temp_mat) > 0 && any(temp_mat)) {
            st_matrix("e(omega)", temp_mat)
        }
    }
    else if (method == "joint") {
        /* Joint: global intercept mu with alpha_1 = beta_1 = 0 */
        mu = _trop_safe_read_scalar("__trop_mu")
        st_numscalar("e(mu)", mu)

        n_treated = _trop_safe_read_scalar("__trop_n_treated")
        st_numscalar("e(N_treated)", n_treated)
        st_numscalar("e(N_treated_obs)", n_treated)

        /* Time weights delta_time (T x 1) and unit weights delta_unit (N x 1) */
        temp_mat = st_matrix("__trop_delta_time")
        if (rows(temp_mat) > 0 && any(temp_mat)) {
            st_matrix("e(delta_time)", temp_mat)
        }
        temp_mat = st_matrix("__trop_delta_unit")
        if (rows(temp_mat) > 0 && any(temp_mat)) {
            st_matrix("e(delta_unit)", temp_mat)
        }
    }

    /* ── diagnostic warnings ─────────────────────────────────────────── */
    _trop_display_warnings(loocv_n_valid, loocv_n_attempted,
        converged, n_iterations, method, n_obs_estimated, n_obs_failed)

    if (bootstrap_reps > 0 && n_bootstrap_valid < bootstrap_reps) {
        _trop_display_bootstrap_warnings(n_bootstrap_valid, bootstrap_reps)
    }

    /* ── command metadata ────────────────────────────────────────────── */
    st_global("e(title)", "TROP Estimator")
    st_global("e(predict)", "trop_p")
}

/*──────────────────────────────────────────────────────────────────────────────
  _trop_safe_read_scalar()

  Reads a Stata scalar by name.  Returns missing (.) when the scalar
  does not exist (st_numscalar() yields a 0x0 matrix in that case).

  Arguments
    name   Stata scalar name

  Returns
    real scalar   the scalar value, or missing
──────────────────────────────────────────────────────────────────────────────*/

real scalar _trop_safe_read_scalar(string scalar name)
{
    real matrix temp
    temp = st_numscalar(name)
    if (rows(temp) == 0) {
        return(.)
    }
    return(temp)
}

/*──────────────────────────────────────────────────────────────────────────────
  _trop_display_warnings()

  Emits diagnostic messages for estimation issues:
    - Per-observation estimation failures (twostep only)
    - Non-convergence warnings

  LOOCV failure-rate warnings are handled separately by
  check_loocv_fail_rate().

  Arguments
    loocv_n_valid      successful LOOCV fits
    loocv_n_attempted  attempted LOOCV fits
    converged          convergence flag (1 = converged)
    n_iterations       iteration count
    method             "twostep" or "joint"
    n_obs_estimated    successfully estimated observations
    n_obs_failed       failed observations
──────────────────────────────────────────────────────────────────────────────*/

void _trop_display_warnings(
    real scalar loocv_n_valid,
    real scalar loocv_n_attempted,
    real scalar converged,
    real scalar n_iterations,
    string scalar method,
    real scalar n_obs_estimated,
    real scalar n_obs_failed
)
{
    if (method == "twostep" && n_obs_failed < . && n_obs_failed > 0) {
        printf("{res}Warning: %g of %g treated observations failed to estimate.{txt}\n",
               n_obs_failed, n_obs_estimated + n_obs_failed)
        printf("{res}         ATT is based on %g successfully estimated observations.{txt}\n",
               n_obs_estimated)
    }

    if (converged == 0) {
        if (method == "twostep") {
            printf("{res}Warning: Not all treated observations converged (max iterations=%g){txt}\n", n_iterations)
        }
        else {
            printf("{res}Warning: Estimation did not converge (iterations=%g){txt}\n", n_iterations)
        }
        printf("{res}         Results may be unreliable.{txt}\n")
    }
}

/*──────────────────────────────────────────────────────────────────────────────
  _trop_compute_effective_rank()

  Computes the continuous effective rank of a matrix via SVD:

      effective_rank = sum(s) / s[1]

  where s is the vector of singular values.  This measures the degree
  to which the low-rank component mu concentrates on a small number of
  factors.

  Arguments
    L   real matrix (T x N factor matrix)

  Returns
    real scalar   effective rank, 0 for empty/zero matrices, . on failure
──────────────────────────────────────────────────────────────────────────────*/

real scalar _trop_compute_effective_rank(real matrix L)
{
    real matrix Lcopy, Vt
    real colvector s
    real scalar s_sum, min_dim, rc, r, c

    if (rows(L) == 0 || cols(L) == 0) return(0)

    min_dim = min((rows(L), cols(L)))
    if (min_dim < 1) return(0)

    Lcopy = L

    /* Replace missing values with zero (SVD requires finite entries) */
    for (r = 1; r <= rows(Lcopy); r++) {
        for (c = 1; c <= cols(Lcopy); c++) {
            if (Lcopy[r, c] >= .) {
                Lcopy[r, c] = 0
            }
        }
    }

    if (max(abs(Lcopy)) == 0) return(0)

    /* _svd() requires rows >= cols; transpose if necessary */
    if (rows(Lcopy) < cols(Lcopy)) {
        Lcopy = Lcopy'
    }
    _svd(Lcopy, s, Vt)
    if (length(s) == 0 || hasmissing(s)) return(.)

    if (length(s) == 0) return(0)
    if (s[1] <= 0) return(0)

    s_sum = sum(s)
    return(s_sum / s[1])
}

/*──────────────────────────────────────────────────────────────────────────────
  _trop_display_bootstrap_warnings()

  Emits warnings or errors based on the bootstrap failure rate:
    > 50%   error message and deferred fatal exit flag
    > 10%   warning that standard errors may be less reliable

  Arguments
    n_valid   number of successful bootstrap replications
    n_total   total number of bootstrap replications attempted
──────────────────────────────────────────────────────────────────────────────*/

void _trop_display_bootstrap_warnings(
    real scalar n_valid,
    real scalar n_total
)
{
    real scalar failure_rate

    if (n_total <= 0) return

    failure_rate = (n_total - n_valid) / n_total

    if (failure_rate > 0.50) {
        errprintf("Error: Bootstrap failure rate exceeds 50%% (%.1f%%)\n", failure_rate * 100)
        errprintf("       %g of %g bootstrap iterations failed.\n",
                  n_total - n_valid, n_total)
        errprintf("       Check data quality or reduce bootstrap replications.\n")
        st_numscalar("__trop_fatal_rc", 504)
    }
    else if (failure_rate > 0.10) {
        printf("{res}Warning: Bootstrap failure rate is %.1f%% (%g/%g successful){txt}\n",
               failure_rate * 100, n_valid, n_total)
        printf("{res}         Standard errors may be less reliable.{txt}\n")
    }
}

/*──────────────────────────────────────────────────────────────────────────────
  _trop_store_lambda_grids()

  Constructs e(lambda_grid) and e(cv_curve) from the individual
  per-dimension lambda grids produced during LOOCV.

  TROP uses cycling (coordinate descent) LOOCV, which does not evaluate
  the full Cartesian product of the grid.  Consequently, e(cv_curve)
  records the LOOCV score Q(lambda_hat) only at the row closest to the
  selected optimal lambdas; all other rows have missing CV loss.

  Arguments
    opt_time    selected lambda_time
    opt_unit    selected lambda_unit
    opt_nn      selected lambda_nn
    opt_score   LOOCV score at the optimum
──────────────────────────────────────────────────────────────────────────────*/

void _trop_store_lambda_grids(
    real scalar opt_time,
    real scalar opt_unit,
    real scalar opt_nn,
    real scalar opt_score
)
{
    real matrix lt_grid, lu_grid, ln_grid
    real scalar n_lt, n_lu, n_ln, n_combo
    real matrix lambda_grid, cv_curve
    real scalar idx, ilt, ilu, iln
    real scalar best_dist, dist, best_row

    lt_grid = st_matrix("__trop_lambda_time_grid")
    lu_grid = st_matrix("__trop_lambda_unit_grid")
    ln_grid = st_matrix("__trop_lambda_nn_grid")

    if (rows(lt_grid) == 0 || rows(lu_grid) == 0 || rows(ln_grid) == 0) {
        return
    }

    n_lt = cols(lt_grid)
    n_lu = cols(lu_grid)
    n_ln = cols(ln_grid)

    /* Cartesian product: e(lambda_grid) is N_combo x 3 */
    n_combo = n_lt * n_lu * n_ln
    lambda_grid = J(n_combo, 3, .)
    cv_curve = J(n_combo, 4, .)

    idx = 0
    for (ilt = 1; ilt <= n_lt; ilt++) {
        for (ilu = 1; ilu <= n_lu; ilu++) {
            for (iln = 1; iln <= n_ln; iln++) {
                idx++
                lambda_grid[idx, 1] = lt_grid[1, ilt]
                lambda_grid[idx, 2] = lu_grid[1, ilu]
                lambda_grid[idx, 3] = ln_grid[1, iln]
                cv_curve[idx, 1] = lt_grid[1, ilt]
                cv_curve[idx, 2] = lu_grid[1, ilu]
                cv_curve[idx, 3] = ln_grid[1, iln]
            }
        }
    }

    /* Locate the grid row closest to the optimal lambdas */
    best_dist = 1e100
    best_row = 1
    for (idx = 1; idx <= n_combo; idx++) {
        dist = (lambda_grid[idx, 1] - opt_time)^2 +
               (lambda_grid[idx, 2] - opt_unit)^2 +
               (lambda_grid[idx, 3] - opt_nn)^2
        if (dist < best_dist) {
            best_dist = dist
            best_row = idx
        }
    }
    cv_curve[best_row, 4] = opt_score

    st_matrix("e(lambda_grid)", lambda_grid)
    st_matrix("e(cv_curve)", cv_curve)
}

end
