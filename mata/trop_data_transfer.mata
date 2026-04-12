/*──────────────────────────────────────────────────────────────────────────────
  trop_data_transfer.mata

  Data transfer layer between Mata and the native-code plugin.

  Stores variable indices, regularization grids, estimation options, and
  pre-allocated output matrices into the Stata namespace (scalars, matrices,
  global macros) so that the plugin can retrieve them via the Stata Plugin
  Interface (SPI).  Observation-level data is read by the plugin directly
  through SF_vdata(); this module does not construct large Mata matrices.

  Contents
    trop_prepare_data()              variable-index and panel-dimension setup
    trop_prepare_lambda_grids()      regularization-parameter grids
    trop_prepare_options()           estimation and inference options
    trop_prepare_output_matrices()   pre-allocate plugin output matrices
    trop_prepare_bootstrap()         bootstrap parameters and output matrix
    trop_cleanup_temp_vars()         drop all __trop_* temporaries
──────────────────────────────────────────────────────────────────────────────*/

version 17
mata:
mata set matastrict on

/*──────────────────────────────────────────────────────────────────────────────
  trop_prepare_data()

  Maps variable names to positional indices within the plugin-call varlist
  and stores them as Stata scalars.  Records panel dimensions (N, T) and
  constructs the T x T time-distance matrix used by the exponential
  time-decay kernel  theta_s(t) = exp(-lambda_time * |t - s|).

  Arguments
    y_varname       outcome variable name          (Y)
    d_varname       binary treatment indicator name (W)
    panel_varname   panel (unit) identifier name
    time_varname    time period identifier name
    n_units         number of cross-sectional units (N)
    n_periods       number of time periods          (T)

  Stored scalars
    __trop_y_varindex, __trop_d_varindex, __trop_ctrl_varindex,
    __trop_panel_varindex, __trop_time_varindex,
    __trop_n_units, __trop_n_periods

  Stored globals
    __trop_varlist, __trop_y_varname, __trop_d_varname,
    __trop_panel_varname, __trop_time_varname

  Stored matrix
    __trop_time_dist   T x T absolute period distances
──────────────────────────────────────────────────────────────────────────────*/

void trop_prepare_data(
    string scalar y_varname,
    string scalar d_varname,
    string scalar panel_varname,
    string scalar time_varname,
    real scalar n_units,
    real scalar n_periods
)
{
    /* Varlist order for `plugin call`: Y, W, panel_id, time_id.
       Variable indices are therefore fixed at 1..4. */
    st_global("__trop_varlist", y_varname + " " + d_varname + " " + panel_varname + " " + time_varname)

    st_numscalar("__trop_y_varindex", 1)
    st_numscalar("__trop_d_varindex", 2)

    /* Control indicator reuses the treatment column; the plugin
       identifies control observations as those with W_{it} = 0. */
    st_numscalar("__trop_ctrl_varindex", 2)

    st_numscalar("__trop_panel_varindex", 3)
    st_numscalar("__trop_time_varindex", 4)

    st_numscalar("__trop_n_units", n_units)
    st_numscalar("__trop_n_periods", n_periods)

    st_global("__trop_y_varname", y_varname)
    st_global("__trop_d_varname", d_varname)
    st_global("__trop_panel_varname", panel_varname)
    st_global("__trop_time_varname", time_varname)

    _trop_prepare_time_dist_matrix(n_periods)
}

/*──────────────────────────────────────────────────────────────────────────────
  trop_prepare_lambda_grids()

  Stores candidate grids for the three regularization parameters
  (lambda_time, lambda_unit, lambda_nn) as Stata matrices, together with
  the corresponding grid lengths as scalars.

  Arguments
    lambda_time_grid   row vector of lambda_time candidates
    lambda_unit_grid   row vector of lambda_unit candidates
    lambda_nn_grid     row vector of lambda_nn   candidates
──────────────────────────────────────────────────────────────────────────────*/

void trop_prepare_lambda_grids(
    real rowvector lambda_time_grid,
    real rowvector lambda_unit_grid,
    real rowvector lambda_nn_grid
)
{
    st_matrix("__trop_lambda_time_grid", lambda_time_grid)
    st_matrix("__trop_lambda_unit_grid", lambda_unit_grid)
    st_matrix("__trop_lambda_nn_grid", lambda_nn_grid)

    st_numscalar("__trop_n_lambda_time", cols(lambda_time_grid))
    st_numscalar("__trop_n_lambda_unit", cols(lambda_unit_grid))
    st_numscalar("__trop_n_lambda_nn", cols(lambda_nn_grid))
}

/*──────────────────────────────────────────────────────────────────────────────
  _trop_prepare_time_dist_matrix()   [private]

  Constructs a T x T symmetric matrix of absolute period distances:

      time_dist[t, s] = |t - s|

  These distances enter the exponential time-decay kernel

      theta_s^{i,t}(lambda) = exp(-lambda_time * |t - s|)

  which down-weights periods far from the target period t.

  Argument
    n_periods   number of time periods (T)
──────────────────────────────────────────────────────────────────────────────*/

void _trop_prepare_time_dist_matrix(real scalar n_periods)
{
    real matrix time_dist
    real scalar t, s

    time_dist = J(n_periods, n_periods, 0)

    for (t = 1; t <= n_periods; t++) {
        for (s = 1; s <= n_periods; s++) {
            time_dist[t, s] = abs(t - s)
        }
    }

    st_matrix("__trop_time_dist", time_dist)
}

/*──────────────────────────────────────────────────────────────────────────────
  trop_prepare_output_matrices()

  Pre-allocates Stata matrices that the plugin fills via SF_mat_store().
  Must be called before `plugin call`.

  Arguments
    n_units     number of units   (N)
    n_periods   number of periods (T)
    n_treated   number of treated observations (optional; default 1)

  Pre-allocated matrices
    __trop_alpha           N x 1          unit fixed effects   (alpha_i)
    __trop_beta            T x 1          time fixed effects   (beta_t)
    __trop_tau             n_treated x 1  treatment effects    (tau_{i,t})
    __trop_factor_matrix   T x N          low-rank component   (mu)
    __trop_theta           T x 1          time weights   (two-step)
    __trop_omega           N x 1          unit weights   (two-step)
    __trop_delta_time      T x 1          time weights   (joint)
    __trop_delta_unit      N x 1          unit weights   (joint)
──────────────────────────────────────────────────────────────────────────────*/

void trop_prepare_output_matrices(
    real scalar n_units,
    real scalar n_periods,
    | real scalar n_treated
)
{
    if (args() < 3 || n_treated < 1) {
        n_treated = 1
    }

    /* Parameter estimates */
    st_matrix("__trop_alpha", J(n_units, 1, 0))
    st_matrix("__trop_beta", J(n_periods, 1, 0))
    st_matrix("__trop_tau", J(n_treated, 1, 0))
    st_matrix("__trop_factor_matrix", J(n_periods, n_units, 0))

    /* Weight vectors for both estimator variants.
       SF_mat_store() requires the target matrix to exist before the
       plugin call. */
    st_matrix("__trop_theta", J(n_periods, 1, 0))
    st_matrix("__trop_omega", J(n_units, 1, 0))
    st_matrix("__trop_delta_time", J(n_periods, 1, 0))
    st_matrix("__trop_delta_unit", J(n_units, 1, 0))
}

/*──────────────────────────────────────────────────────────────────────────────
  trop_prepare_options()

  Stores estimation options as Stata scalars for the plugin.

  Arguments
    max_loocv_samples   maximum LOOCV evaluation points
    max_iter            maximum iterations for alternating minimization
    tol                 convergence tolerance
    seed                PRNG seed for bootstrap
    bootstrap_reps      number of bootstrap replications
    alpha_level         significance level (e.g. 0.05 for 95% CI)
    verbose             verbosity level (0 = silent)
──────────────────────────────────────────────────────────────────────────────*/

void trop_prepare_options(
    real scalar max_loocv_samples,
    real scalar max_iter,
    real scalar tol,
    real scalar seed,
    real scalar bootstrap_reps,
    real scalar alpha_level,
    real scalar verbose
)
{
    st_numscalar("__trop_max_loocv_samples", max_loocv_samples)
    st_numscalar("__trop_max_iter", max_iter)
    st_numscalar("__trop_tol", tol)
    st_numscalar("__trop_seed", seed)
    st_numscalar("__trop_n_bootstrap", bootstrap_reps)
    st_numscalar("__trop_alpha_level", alpha_level)
    st_numscalar("__trop_verbose", verbose)
}

/*──────────────────────────────────────────────────────────────────────────────
  trop_prepare_bootstrap()

  Stores bootstrap parameters and pre-allocates the output matrix that
  collects replicate estimates.

  Arguments
    n_bootstrap   number of bootstrap replications
    alpha         significance level for confidence intervals
    seed          PRNG seed
    lambda_time   selected lambda_time
    lambda_unit   selected lambda_unit
    lambda_nn     selected lambda_nn
    max_iter      maximum iterations per replicate
    tol           convergence tolerance per replicate

  Pre-allocated matrix
    __trop_bootstrap_estimates   n_bootstrap x 1, initialised to missing
──────────────────────────────────────────────────────────────────────────────*/

void trop_prepare_bootstrap(
    real scalar n_bootstrap,
    real scalar alpha,
    real scalar seed,
    real scalar lambda_time,
    real scalar lambda_unit,
    real scalar lambda_nn,
    real scalar max_iter,
    real scalar tol
)
{
    st_numscalar("__trop_n_bootstrap", n_bootstrap)
    /* Named __trop_bs_alpha to avoid collision with the unit-fixed-effects
       matrix __trop_alpha. */
    st_numscalar("__trop_bs_alpha", alpha)
    st_numscalar("__trop_seed", seed)

    st_numscalar("__trop_lambda_time", lambda_time)
    st_numscalar("__trop_lambda_unit", lambda_unit)
    st_numscalar("__trop_lambda_nn", lambda_nn)

    st_numscalar("__trop_max_iter", max_iter)
    st_numscalar("__trop_tol", tol)

    st_matrix("__trop_bootstrap_estimates", J(n_bootstrap, 1, .))
}

/*──────────────────────────────────────────────────────────────────────────────
  trop_cleanup_temp_vars()

  Drops every __trop_* scalar, matrix, and global macro from the Stata
  namespace.  Uses st_dir() to enumerate names programmatically so that
  newly introduced temporaries are covered without code changes.
──────────────────────────────────────────────────────────────────────────────*/

void trop_cleanup_temp_vars()
{
    string colvector names
    real scalar i

    /* Scalars */
    names = st_dir("global", "numscalar", "__trop_*")
    for (i = 1; i <= length(names); i++) {
        stata("capture scalar drop " + names[i])
    }

    /* Matrices */
    names = st_dir("global", "matrix", "__trop_*")
    for (i = 1; i <= length(names); i++) {
        stata("capture matrix drop " + names[i])
    }

    /* Global macros */
    names = st_dir("global", "macro", "__trop_*")
    for (i = 1; i <= length(names); i++) {
        st_global(names[i], "")
    }
    st_global("__trop_touse_var", "")
    st_global("__trop_varlist", "")
    st_global("__trop_y_varname", "")
    st_global("__trop_d_varname", "")
    st_global("__trop_panel_varname", "")
    st_global("__trop_time_varname", "")

    /* Clean up the shared random seed scalar (not __trop_ prefixed) */
    stata("capture scalar drop TROP_GLOBAL_SEED")
}

end
