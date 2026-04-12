*! Counterfactual outcome prediction for trop
/*
    trop_predict_y0

    Generates predicted counterfactual outcome Y(0) for each observation.

    The additive model is:

        Y(0) = mu + alpha_i + beta_t + L_{t,i}

    where alpha_i are unit fixed effects, beta_t are time fixed effects,
    and L_{t,i} is the low-rank interaction component.

    Two-step method:
        Treated units:  Y(0) = Y_obs - tau_hat  (exact decomposition)
        Control units:  Y(0) = alpha_i + beta_t + L_{t,i}

    Joint method:
        All units:      Y(0) = mu + alpha_i + beta_t + L_{t,i}

    Required stored estimates:
        e(alpha)          unit fixed effects (N x 1)
        e(beta)           time fixed effects (T x 1)
        e(factor_matrix)  low-rank interaction (T x N, column-major)
        e(mu)             intercept (joint method)
        e(tau)            observation-level treatment effects (two-step method)
*/


program define trop_predict_y0
    version 17
    syntax newvarname [if] [in]
    
    marksample touse, novarlist
    
    
    // Initialize prediction variable.
    qui gen double `varlist' = .
    
    // Construct panel and time indices consistent with estimation sample.
    tempvar panel_idx time_idx
    qui egen `panel_idx' = group(`e(panelvar)') if e(sample)
    qui egen `time_idx' = group(`e(timevar)') if e(sample)
    
    // Retrieve estimation method; default is two-step.
    local method "`e(method)'"
    if "`method'" == "" {
        local method "twostep"
    }
    
    // Dispatch to Mata.
    local treatvar "`e(treatvar)'"
    local depvar "`e(depvar)'"
    mata: _trop_predict_y0("`varlist'", "`panel_idx'", "`time_idx'", ///
        "`method'", "`touse'", "`treatvar'", "`depvar'")
    
    label variable `varlist' "Counterfactual prediction Y(0)"
end

// ---------------------------------------------------------------------------
// Mata: counterfactual prediction
// ---------------------------------------------------------------------------
version 17
mata:
mata set matastrict on

void _trop_predict_y0(
    string scalar varname,
    string scalar panel_var,
    string scalar time_var,
    string scalar method,
    string scalar touse_var,
    string scalar treatvar,
    string scalar depvar
)
{
    real colvector panel_idx, time_idx
    real colvector alpha, beta
    real matrix L
    real colvector pred
    real scalar mu, n, i, i_idx, t_idx
    
    panel_idx = st_data(., panel_var, touse_var)
    time_idx = st_data(., time_var, touse_var)
    
    alpha = st_matrix("e(alpha)")
    beta = st_matrix("e(beta)")
    
    // T x N low-rank interaction matrix, column-major layout.
    L = st_matrix("e(factor_matrix)")
    
    n = rows(panel_idx)
    pred = J(n, 1, .)
    
    if (method == "twostep") {
        _trop_predict_y0_twostep(varname, panel_var, time_var, 
            touse_var, treatvar, depvar, alpha, beta, L, n)
    }
    else {
        // Joint method: reconstruct Y(0) with global intercept.
        mu = st_numscalar("e(mu)")
        if (mu == .) {
            mu = 0
            printf("{txt}(note: e(mu) not found; defaulting to 0)\n")
        }
        
        for (i = 1; i <= n; i++) {
            i_idx = panel_idx[i]
            t_idx = time_idx[i]
            
            if (i_idx < 1 | i_idx > rows(alpha)) continue
            if (t_idx < 1 | t_idx > rows(beta)) continue
            if (t_idx > rows(L) | i_idx > cols(L)) continue
            
            pred[i] = mu + alpha[i_idx] + beta[t_idx] + L[t_idx, i_idx]
        }
        
        st_store(., varname, touse_var, pred)
    }
}

void _trop_predict_y0_twostep(
    string scalar varname,
    string scalar panel_var,
    string scalar time_var,
    string scalar touse_var,
    string scalar treatvar,
    string scalar depvar,
    real colvector alpha,
    real colvector beta,
    real matrix L,
    real scalar n
)
{
/*
    Two-step counterfactual prediction.

    Treated observations:
        Y(0) = Y_obs - tau_hat
        This preserves the exact identity Y_obs - Y(0) = tau_hat.

    Control observations:
        Y(0) = alpha[i] + beta[t] + L[t,i]

    e(tau) is stored in time-major order (t = 1..T within i = 1..N for
    treated cells). Treated observations are sorted accordingly before
    tau assignment.
*/
    real colvector panel_idx, time_idx, treat, y_data, pred
    real matrix tau_vec
    real scalar i, i_idx, t_idx, k, n_tau, n_treated, obs_i
    real matrix treated_info
    
    panel_idx = st_data(., panel_var, touse_var)
    time_idx = st_data(., time_var, touse_var)
    treat = st_data(., treatvar, touse_var)
    y_data = st_data(., depvar, touse_var)
    
    pred = J(n, 1, .)
    
    tau_vec = st_matrix("e(tau)")
    n_tau = rows(tau_vec)
    
    // Control observations: reconstruct from additive model.
    for (i = 1; i <= n; i++) {
        i_idx = panel_idx[i]
        t_idx = time_idx[i]
        
        if (i_idx < . && t_idx < . && treat[i] != 1) {
            if (i_idx >= 1 & i_idx <= rows(alpha) & 
                t_idx >= 1 & t_idx <= rows(beta) &
                t_idx <= rows(L) & i_idx <= cols(L)) {
                pred[i] = alpha[i_idx] + beta[t_idx] + L[t_idx, i_idx]
            }
        }
    }
    
    // Treated observations: subtract observation-level treatment effect.
    if (n_tau > 0) {
        n_treated = 0
        for (i = 1; i <= n; i++) {
            if (treat[i] == 1 && panel_idx[i] < . && time_idx[i] < .) {
                n_treated++
            }
        }
        
        if (n_treated == n_tau) {
            // Collect treated observations with row indices.
            treated_info = J(n_treated, 3, .)
            k = 0
            for (i = 1; i <= n; i++) {
                if (treat[i] == 1 && panel_idx[i] < . && time_idx[i] < .) {
                    k++
                    treated_info[k, 1] = time_idx[i]
                    treated_info[k, 2] = panel_idx[i]
                    treated_info[k, 3] = i
                }
            }
            
            // Sort by (time, panel) to align with time-major order of e(tau).
            treated_info = sort(treated_info, (1, 2))
            
            // Assign Y(0) = Y_obs - tau_hat.
            for (k = 1; k <= n_treated; k++) {
                obs_i = treated_info[k, 3]
                if (y_data[obs_i] < . && tau_vec[k, 1] < .) {
                    pred[obs_i] = y_data[obs_i] - tau_vec[k, 1]
                }
            }
        }
        else {
            // Count mismatch: fall back to additive model for treated units.
            printf("{txt}(note: treated count %g differs from e(tau) length %g; " +
                   "using additive model)\n", n_treated, n_tau)
            for (i = 1; i <= n; i++) {
                if (treat[i] == 1) {
                    i_idx = panel_idx[i]
                    t_idx = time_idx[i]
                    if (i_idx >= 1 & i_idx <= rows(alpha) & 
                        t_idx >= 1 & t_idx <= rows(beta) &
                        t_idx <= rows(L) & i_idx <= cols(L)) {
                        pred[i] = alpha[i_idx] + beta[t_idx] + L[t_idx, i_idx]
                    }
                }
            }
        }
    }
    else {
        // No observation-level effects available; use additive model.
        for (i = 1; i <= n; i++) {
            if (treat[i] == 1) {
                i_idx = panel_idx[i]
                t_idx = time_idx[i]
                if (i_idx >= 1 & i_idx <= rows(alpha) & 
                    t_idx >= 1 & t_idx <= rows(beta) &
                    t_idx <= rows(L) & i_idx <= cols(L)) {
                    pred[i] = alpha[i_idx] + beta[t_idx] + L[t_idx, i_idx]
                }
            }
        }
    }
    
    st_store(., varname, touse_var, pred)
}

end
