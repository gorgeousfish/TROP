/*──────────────────────────────────────────────────────────────────────────────
  trop_rust_interface.mata

  Plugin interface layer for the TROP estimator.

  Mata wrappers for the compiled native plugin.  This module handles
  platform detection, plugin discovery and loading, and thin call-through
  functions for every plugin entry point.  No estimation logic resides
  here; all numerical computation is delegated to the plugin binary.

  Contents
  ────────
  1. Platform detection and plugin loading
       _trop_get_plugin_name()        platform-specific plugin filename
       _trop_get_plugin_path()        locate plugin binary on disk
       _trop_load_plugin()            load plugin into Stata
       _trop_call_plugin()            issue a plugin call with varlist
       _trop_call_plugin_simple()     issue a plugin call without varlist
       trop_rust_available()          query plugin availability

  2. LOOCV interface
       trop_loocv_twostep()             two-step LOOCV grid search
       trop_loocv_joint()               joint LOOCV coordinate-descent search
       trop_loocv_joint_exhaustive()    joint LOOCV Cartesian grid search

  3. Estimation interface
       trop_estimate_twostep()        two-step point estimation
       trop_estimate_joint()          joint point estimation

  4. Bootstrap interface
       trop_bootstrap_twostep()       two-step bootstrap variance
       trop_bootstrap_joint()         joint bootstrap variance

  5. Distance matrix
       trop_compute_unit_distance()   inter-unit RMSE distance matrix

  6. LOOCV diagnostics
       store_loocv_diagnostics()      post LOOCV diagnostics to e()
       check_loocv_fail_rate()        abort or warn on high failure rate
       display_loocv_verbose()        print LOOCV summary to console

  7. Orchestration
       _trop_main()                   LOOCV -> estimation -> bootstrap
──────────────────────────────────────────────────────────────────────────────*/

version 17
mata:
mata set matastrict on

/* ═══════════════════════════════════════════════════════════════════════════
   1. Platform detection and plugin loading
   ═══════════════════════════════════════════════════════════════════════════ */

/*──────────────────────────────────────────────────────────────────────────────
  _trop_get_plugin_name()

  Maps the current OS and CPU architecture to the corresponding plugin
  filename.

  Supported targets
    macOS  Apple Silicon   trop_macos_arm64.plugin
    macOS  Intel x86-64    trop_macos_x64.plugin
    Linux  x86-64          trop_linux_x64.plugin
    Windows x86-64         trop_windows_x64.plugin

  Returns "" when the platform cannot be identified.
──────────────────────────────────────────────────────────────────────────────*/

string scalar _trop_get_plugin_name()
{
    string scalar os, machine

    os      = c("os")
    machine = c("machine_type")

    /* macOS: c(os) may report "MacOSX" or "Unix" across Stata versions */
    if (os == "MacOSX" | (os == "Unix" & strpos(machine, "Mac") > 0)) {
        if (strpos(machine, "Apple Silicon") > 0 |
            strpos(machine, "arm64") > 0) {
            return("trop_macos_arm64.plugin")
        }
        return("trop_macos_x64.plugin")
    }

    /* Linux: Unix that is not macOS */
    if (os == "Unix" & strpos(machine, "Mac") == 0) {
        return("trop_linux_x64.plugin")
    }

    /* Windows */
    if (os == "Windows") {
        return("trop_windows_x64.plugin")
    }

    return("")
}

/*──────────────────────────────────────────────────────────────────────────────
  _trop_get_plugin_path()

  Searches for the plugin binary in the following order:
    1. Relative to trop.ado via adopath (CWD-independent)
    2. <cwd>/trop_stata/plugin/
    3. <cwd>/plugin/
    4. <cwd>/../plugin/
    5. PLUS   sysdir  c(sysdir_plus)t/
    6. PERSONAL sysdir  c(sysdir_personal)
    7. Current working directory

  Returns the full path if found, "" otherwise.
──────────────────────────────────────────────────────────────────────────────*/

string scalar _trop_get_plugin_path()
{
    string scalar plugin_name, plugin_path
    string scalar sysdir_plus, sysdir_personal, pwd
    string scalar ado_path, pkg_dir

    plugin_name = _trop_get_plugin_name()
    if (plugin_name == "") return("")

    /* adopath-based lookup */
    (void) _stata("capture findfile trop.ado", 1)
    ado_path = st_global("r(fn)")
    if (ado_path != "") {
        pkg_dir = subinstr(ado_path, "/ado/trop.ado", "", 1)
        if (pkg_dir != ado_path) {
            plugin_path = pkg_dir + "/plugin/" + plugin_name
            if (fileexists(plugin_path)) return(plugin_path)
        }
    }

    /* CWD-relative paths */
    pwd = st_global("c(pwd)")

    plugin_path = pwd + "/trop_stata/plugin/" + plugin_name
    if (fileexists(plugin_path)) return(plugin_path)

    plugin_path = pwd + "/plugin/" + plugin_name
    if (fileexists(plugin_path)) return(plugin_path)

    plugin_path = pwd + "/../plugin/" + plugin_name
    if (fileexists(plugin_path)) return(plugin_path)

    /* System directories */
    sysdir_plus = st_global("c(sysdir_plus)")
    plugin_path = sysdir_plus + "t/" + plugin_name
    if (fileexists(plugin_path)) return(plugin_path)

    sysdir_personal = st_global("c(sysdir_personal)")
    plugin_path = sysdir_personal + plugin_name
    if (fileexists(plugin_path)) return(plugin_path)

    plugin_path = pwd + "/" + plugin_name
    if (fileexists(plugin_path)) return(plugin_path)

    return("")
}

/*──────────────────────────────────────────────────────────────────────────────
  _trop_load_plugin()

  Loads the plugin into Stata via `program ... plugin using()`.
  If the plugin is already loaded (rc 110), the error is silently
  absorbed.  The plugin need only be loaded once per session.

  Returns 0 on success, 601 if the binary cannot be found.
──────────────────────────────────────────────────────────────────────────────*/

real scalar _trop_load_plugin()
{
    string scalar plugin_path, cmd

    plugin_path = _trop_get_plugin_path()
    if (plugin_path == "") {
        errprintf("plugin not found; check installation\n")
        return(601)
    }

    cmd = `"capture program _trop_plugin, plugin using(""' ///
        + plugin_path + `"")"'
    (void) _stata(cmd, 1)

    return(0)
}

/*──────────────────────────────────────────────────────────────────────────────
  _trop_call_plugin()

  Dispatches a command to the plugin with an optional variable list.

  Arguments
    command   name forwarded to the plugin dispatcher
    varlist   space-separated variable names (optional; defaults to
              the global macro __trop_varlist)

  The plugin receives variables indexed relative to the supplied varlist,
  not by their absolute dataset position.  When the global
  __trop_touse_var is set, an `if` qualifier is appended so that
  SF_ifobs() filters observations accordingly.

  Returns 0 on success.
──────────────────────────────────────────────────────────────────────────────*/

real scalar _trop_call_plugin(string scalar command, | string scalar varlist)
{
    string scalar cmd, vars, touse_var
    real scalar rc

    rc = _trop_load_plugin()
    if (rc != 0) return(rc)

    if (args() < 2 | varlist == "") {
        vars = st_global("__trop_varlist")
        if (vars == "") {
            errprintf("no variable list specified\n")
            return(198)
        }
    }
    else {
        vars = varlist
    }

    touse_var = st_global("__trop_touse_var")
    if (touse_var != "") {
        cmd = "plugin call _trop_plugin " + vars ///
            + " if " + touse_var + " , " + command
    }
    else {
        cmd = "plugin call _trop_plugin " + vars + " , " + command
    }
    rc = _stata(cmd, 1)

    return(rc)
}

/*──────────────────────────────────────────────────────────────────────────────
  _trop_call_plugin_simple()

  Dispatches a command to the plugin without a variable list.

  Arguments
    command   name forwarded to the plugin dispatcher

  Returns 0 on success.
──────────────────────────────────────────────────────────────────────────────*/

real scalar _trop_call_plugin_simple(string scalar command)
{
    string scalar cmd
    real scalar rc

    rc = _trop_load_plugin()
    if (rc != 0) return(rc)

    cmd = "plugin call _trop_plugin , " + command
    rc = _stata(cmd, 1)

    return(rc)
}

/*──────────────────────────────────────────────────────────────────────────────
  trop_rust_available()

  Tests whether the plugin binary exists and can be loaded.

  Returns 1 if available, 0 otherwise.
──────────────────────────────────────────────────────────────────────────────*/

real scalar trop_rust_available()
{
    if (_trop_get_plugin_path() == "") return(0)
    return(_trop_load_plugin() == 0)
}

/* ═══════════════════════════════════════════════════════════════════════════
   2. LOOCV interface

   Leave-one-out cross-validation selects the regularization triplet
   (lambda_time, lambda_unit, lambda_nn) by minimising the criterion

       Q(lambda) = sum_{i,t: W_{it}=0} (tau_hat_{it}(lambda))^2

   over a grid of candidate values.
   ═══════════════════════════════════════════════════════════════════════════ */

/*──────────────────────────────────────────────────────────────────────────────
  trop_loocv_twostep()

  Two-step LOOCV grid search.  For each control observation (i,t) and
  each candidate lambda, solves the leave-one-out penalised regression
  to obtain tau_hat_{it}(lambda), then evaluates Q(lambda).

  Returns 0 on success.
──────────────────────────────────────────────────────────────────────────────*/

real scalar trop_loocv_twostep()
{
    return(_trop_call_plugin("loocv_twostep"))
}

/*──────────────────────────────────────────────────────────────────────────────
  trop_loocv_joint()

  Joint LOOCV coordinate-descent search (default).

  Two-stage strategy adapted from Footnote 2 of Athey et al. (2025):
    Stage 1 — univariate sweeps with extreme fixed values;
    Stage 2 — cyclic coordinate descent until convergence.
  Complexity O(|grid| * max_cycles); favoured for large grids.

  Returns 0 on success.
──────────────────────────────────────────────────────────────────────────────*/

real scalar trop_loocv_joint()
{
    return(_trop_call_plugin("loocv_joint"))
}

/*──────────────────────────────────────────────────────────────────────────────
  trop_loocv_joint_exhaustive()

  Joint LOOCV exhaustive (Cartesian) grid search.

  Evaluates every (lambda_time, lambda_unit, lambda_nn) combination in
  parallel; complexity is O(|Lambda_time| * |Lambda_unit| * |Lambda_nn|).
  Matches the Python reference (diff_diff.trop_global, v3.1.1) exactly.

  Returns 0 on success.
──────────────────────────────────────────────────────────────────────────────*/

real scalar trop_loocv_joint_exhaustive()
{
    return(_trop_call_plugin("loocv_joint_exhaustive"))
}

/* ═══════════════════════════════════════════════════════════════════════════
   3. Estimation interface

   Given the selected lambda_hat, estimate the treatment effect(s).
   The two-step method yields observation-level tau_{i,t}; the joint
   method yields a single scalar tau.
   ═══════════════════════════════════════════════════════════════════════════ */

/*──────────────────────────────────────────────────────────────────────────────
  trop_estimate_twostep()

  Two-step point estimation of observation-level treatment effects.

  Returns 0 on success.
──────────────────────────────────────────────────────────────────────────────*/

real scalar trop_estimate_twostep()
{
    return(_trop_call_plugin("estimate_twostep"))
}

/*──────────────────────────────────────────────────────────────────────────────
  trop_estimate_joint()

  Joint point estimation of the average treatment effect on the treated.

  Returns 0 on success.
──────────────────────────────────────────────────────────────────────────────*/

real scalar trop_estimate_joint()
{
    return(_trop_call_plugin("estimate_joint"))
}

/* ═══════════════════════════════════════════════════════════════════════════
   4. Bootstrap interface

   Variance estimation via unit-level resampling.  Each replicate
   re-draws N units with replacement and re-estimates tau, yielding
   an empirical distribution from which standard errors and percentile
   confidence intervals are computed.
   ═══════════════════════════════════════════════════════════════════════════ */

/*──────────────────────────────────────────────────────────────────────────────
  trop_bootstrap_twostep()

  Two-step bootstrap variance estimation.

  Returns 0 on success.
──────────────────────────────────────────────────────────────────────────────*/

real scalar trop_bootstrap_twostep()
{
    return(_trop_call_plugin("bootstrap_twostep"))
}

/*──────────────────────────────────────────────────────────────────────────────
  trop_bootstrap_joint()

  Joint bootstrap variance estimation.

  Returns 0 on success.
──────────────────────────────────────────────────────────────────────────────*/

real scalar trop_bootstrap_joint()
{
    return(_trop_call_plugin("bootstrap_joint"))
}

/* ═══════════════════════════════════════════════════════════════════════════
   5. Distance matrix

   The unit distance metric is

       dist_{-t}(j,i) = sqrt( sum_{u!=t} (1-W_{iu})(1-W_{ju})
                               (Y_{iu}-Y_{ju})^2
                             / sum_{u!=t} (1-W_{iu})(1-W_{ju}) )

   which measures the root-mean-square outcome difference over jointly
   observed control periods, excluding the target period t.
   ═══════════════════════════════════════════════════════════════════════════ */

/*──────────────────────────────────────────────────────────────────────────────
  trop_compute_unit_distance()

  Computes the N x N inter-unit distance matrix and stores it in the
  Stata matrix __trop_unit_dist.

  Returns 0 on success.
──────────────────────────────────────────────────────────────────────────────*/

real scalar trop_compute_unit_distance()
{
    return(_trop_call_plugin("distance_matrix"))
}

/* ═══════════════════════════════════════════════════════════════════════════
   6. LOOCV diagnostics
   ═══════════════════════════════════════════════════════════════════════════ */

/*──────────────────────────────────────────────────────────────────────────────
  store_loocv_diagnostics()

  Collects LOOCV diagnostic scalars written by the plugin, derives the
  failure rate, and posts everything to e():

    e(loocv_score)        optimal cross-validation score  Q(lambda_hat)
    e(loocv_n_valid)      successful leave-one-out fits
    e(loocv_n_attempted)  attempted leave-one-out fits (=total D=0 cells)
    e(loocv_fail_rate)    fraction of failed fits
    e(seed)               RNG seed (used for bootstrap)

  Arguments
    method  "twostep" or "joint"
    seed    RNG seed (recorded only; LOOCV itself is deterministic)
──────────────────────────────────────────────────────────────────────────────*/

void store_loocv_diagnostics(
    string scalar method,
    real scalar seed)
{
    real scalar loocv_score, n_valid, n_attempted, fail_rate
    real matrix temp

    temp = st_numscalar("__trop_loocv_score")
    if (rows(temp) > 0) loocv_score = temp
    else loocv_score = .

    temp = st_numscalar("__trop_loocv_n_valid")
    if (rows(temp) > 0) n_valid = temp
    else n_valid = .

    temp = st_numscalar("__trop_loocv_n_attempted")
    if (rows(temp) > 0) n_attempted = temp
    else n_attempted = .

    if (n_attempted == . | n_attempted == 0) {
        fail_rate = .
    }
    else {
        fail_rate = (n_attempted - n_valid) / n_attempted
    }

    st_numscalar("e(loocv_score)",       loocv_score)
    st_numscalar("e(loocv_n_valid)",     n_valid)
    st_numscalar("e(loocv_n_attempted)", n_attempted)
    st_numscalar("e(loocv_fail_rate)",   fail_rate)
    st_numscalar("e(seed)",              seed)
}

/*──────────────────────────────────────────────────────────────────────────────
  check_loocv_fail_rate()

  Inspects e(loocv_fail_rate) and acts accordingly:
    > 50%    abort with rc 498 (results unreliable)
    > 10%    issue a warning, continue
    <= 10%   silent

  Returns 0 to continue, 498 to abort.
──────────────────────────────────────────────────────────────────────────────*/

real scalar check_loocv_fail_rate()
{
    real scalar fail_rate, first_t, first_i
    real matrix temp

    temp = st_numscalar("e(loocv_fail_rate)")
    if (rows(temp) == 0) return(0)
    fail_rate = temp

    first_t = .
    first_i = .
    temp = st_numscalar("e(loocv_first_failed_t)")
    if (rows(temp) > 0) first_t = temp
    temp = st_numscalar("e(loocv_first_failed_i)")
    if (rows(temp) > 0) first_i = temp

    if (fail_rate > 0.50) {
        errprintf("LOOCV failure rate too high: %5.1f%%\n", fail_rate * 100)
        errprintf("results would be unreliable; check data quality\n")
        if (first_t != . & first_i != . & first_t >= 0 & first_i >= 0) {
            /* Rust indices are 0-based; display 1-based to match Stata. */
            errprintf(
                "first failing LOO fit: period %g, unit %g (1-based);\n",
                first_t + 1, first_i + 1)
            errprintf(
                "  cross-reference e(panelvar), e(timevar), and the estimation sample\n")
        }
        return(498)
    }
    if (fail_rate > 0.10) {
        printf("{txt}warning: LOOCV failure rate is %5.1f%%\n", fail_rate * 100)
        if (first_t != . & first_i != . & first_t >= 0 & first_i >= 0) {
            printf(
                "{txt}  first failing LOO fit at period %g, unit %g (1-based);\n",
                first_t + 1, first_i + 1)
            printf(
                "{txt}  see e(loocv_first_failed_t), e(loocv_first_failed_i)\n")
        }
    }

    return(0)
}

/*──────────────────────────────────────────────────────────────────────────────
  display_loocv_verbose()

  Prints a one-line summary of LOOCV activity to the Results window:
  number of control observations evaluated and the success/failure
  breakdown.  LOOCV always sums over every D=0 cell (paper Eq. 5).
──────────────────────────────────────────────────────────────────────────────*/

void display_loocv_verbose()
{
    real scalar n_attempted, n_valid, fail_pct
    real matrix temp

    temp = st_numscalar("e(loocv_n_attempted)")
    if (rows(temp) > 0) n_attempted = temp
    else n_attempted = .

    temp = st_numscalar("e(loocv_n_valid)")
    if (rows(temp) > 0) n_valid = temp
    else n_valid = .

    if (n_attempted == . | n_attempted == 0) return

    printf("LOOCV: %g control obs evaluated\n", n_attempted)

    if (n_valid == .) n_valid = n_attempted
    fail_pct = (n_attempted - n_valid) / n_attempted * 100
    printf("LOOCV: %g/%g fits successful (%5.1f%% failure rate)\n",
           n_valid, n_attempted, fail_pct)
}

/* ═══════════════════════════════════════════════════════════════════════════
   7. Orchestration
   ═══════════════════════════════════════════════════════════════════════════ */

/*──────────────────────────────────────────────────────────────────────────────
  _trop_main()

  Executes the full estimation pipeline:
    1. Verify plugin availability
    2. LOOCV hyperparameter search   (if requested)
    3. Point estimation
    4. Bootstrap variance estimation (if requested)

  The e() result storage is handled by the calling ado-file after
  ereturn post.

  Arguments
    method           "twostep" or "joint"
    do_loocv         1 to run LOOCV, 0 to skip
    do_bootstrap     1 to run bootstrap, 0 to skip
    bootstrap_reps   number of bootstrap replications

  Returns 0 on success.
──────────────────────────────────────────────────────────────────────────────*/

real scalar _trop_main(
    string scalar method,
    real scalar do_loocv,
    real scalar do_bootstrap,
    real scalar bootstrap_reps)
{
    real scalar rc
    real scalar lambda_time, lambda_unit, lambda_nn
    real scalar max_iter, tol, seed, alpha_level
    string scalar joint_mode

    /* ── plugin check ────────────────────────────────────────────────── */
    if (!trop_rust_available()) {
        errprintf("plugin not available; reinstall the package\n")
        return(198)
    }

    /* ── LOOCV ───────────────────────────────────────────────────────── */
    if (do_loocv) {
        if (method == "twostep") {
            rc = trop_loocv_twostep()
        }
        else {
            /* joint: dispatch on __trop_joint_loocv_mode.
               "exhaustive" -> full Cartesian search;
               anything else (including missing) -> coordinate descent. */
            joint_mode = st_global("__trop_joint_loocv_mode")
            if (joint_mode == "exhaustive") rc = trop_loocv_joint_exhaustive()
            else rc = trop_loocv_joint()
        }
        if (rc != 0) return(rc)
    }

    /* ── point estimation ────────────────────────────────────────────── */
    if (method == "twostep") rc = trop_estimate_twostep()
    else rc = trop_estimate_joint()
    if (rc != 0) return(rc)

    /* ── bootstrap variance estimation ───────────────────────────────── */
    if (do_bootstrap & bootstrap_reps > 0) {
        lambda_time = st_numscalar("__trop_lambda_time")
        lambda_unit = st_numscalar("__trop_lambda_unit")
        lambda_nn   = st_numscalar("__trop_lambda_nn")
        max_iter    = st_numscalar("__trop_max_iter")
        tol         = st_numscalar("__trop_tol")
        seed        = st_numscalar("__trop_seed")
        alpha_level = st_numscalar("__trop_alpha_level")

        trop_prepare_bootstrap(
            bootstrap_reps, alpha_level, seed,
            lambda_time, lambda_unit, lambda_nn,
            max_iter, tol)

        if (method == "twostep") rc = trop_bootstrap_twostep()
        else rc = trop_bootstrap_joint()
        if (rc != 0) return(rc)
    }

    return(0)
}

end
