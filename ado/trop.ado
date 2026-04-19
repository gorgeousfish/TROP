*! trop -- Triply Robust Panel (TROP) estimator

/*
    trop -- Estimate average treatment effects on the treated (ATT) in panel
    data using the TROP framework: unit weights, time weights, and a
    nuclear-norm-penalized low-rank regression adjustment.

    Hyperparameters (lambda_time, lambda_unit, lambda_nn) are selected via
    leave-one-out cross-validation (LOOCV) unless fixedlambda() is specified.
    Inference is obtained through a cluster bootstrap that resamples units
    within the treated and control groups separately.

    Syntax
    ------
    trop depvar treatvar [if] [in], panelvar(varname) timevar(varname)
        [method(twostep|joint|local|global) grid_style(default|extended)
         joint_loocv(cycling|exhaustive)
         lambda_time_grid(numlist) lambda_unit_grid(numlist)
         lambda_nn_grid(numlist) fixedlambda(numlist)
         tol(real) maxiter(integer)
         bootstrap(integer) bsalpha(real) seed(integer)
         verbose level(cilevel)]
*/


program define trop, eclass
    version 17.0

    // --- Version subcommand -----------------------------------------------
    // Intercept `trop, version` before the main syntax parser, which
    // requires exactly two variables.  This allows users to verify the
    // installed version without supplying a varlist.
    if `"`0'"' == ", version" | `"`0'"' == ",version" | `"`0'"' == "version" {
        di as txt "trop version 1.0.0"
        di as txt "Triply Robust Panel Estimator"
        di as txt "Athey, Imbens, Qu & Viviano (2025)"
        di as txt ""
        di as txt "Stata implementation by Xuanyu Cai & Wenli Xu"
        di as txt "License: AGPL-3.0"
        exit
    }

    // --- Syntax parsing ---------------------------------------------------
    syntax varlist(min=2 max=2) [if] [in], ///
        Panelvar(varname)                   /// panel (unit) identifier
        Timevar(varname)                    /// time period identifier
        [                                   ///
        Method(string)                      /// "twostep" (per-obs tau) or "joint" (scalar tau)
        Grid_style(string)                  /// grid layout for LOOCV: "default" or "extended"
        JOint_loocv(string)                 /// joint LOOCV strategy: "cycling" (default) or "exhaustive"
        LAMbda_time_grid(numlist)                   /// user-supplied lambda_time grid
        LAMbda_unit_grid(numlist)                   /// user-supplied lambda_unit grid
        LAMbda_nn_grid(numlist missingokay)         /// user-supplied lambda_nn grid (`.` = inf)
        FIXEDlambda(numlist missingokay min=3 max=3) /// fixed (l_time l_unit l_nn); `.` = inf for l_nn
        TOL(real 1e-6)                      /// convergence tolerance for iterative estimation
        MAXiter(integer 100)                /// maximum iterations for iterative estimation
        BOOTstrap(integer 200)              /// number of bootstrap replications (paper Alg 3 default; 0 = skip)
        BSalpha(real -1)                    /// deprecated; retained for backward compatibility
        SEED(integer 42)                    /// RNG seed for bootstrap
        VERbose                             /// display progress and diagnostic messages
        Level(cilevel)                      /// confidence level for bootstrap CI
        ]

    // --- Variable extraction ----------------------------------------------
    gettoken depvar treatvar : varlist
    // gettoken preserves the separating whitespace; strip it so downstream
    // consumers (e(treatvar), Mata _st_varindex, etc.) see a clean name.
    local treatvar = trim("`treatvar'")

    // --- Default values ---------------------------------------------------
    if "`method'" == "" {
        local method "twostep"
    }
    // Accept paper terminology as aliases:
    //   method("local")  == method("twostep")  (per-observation weights)
    //   method("global") == method("joint")    (single scalar tau)
    if "`method'" == "local" {
        local method "twostep"
    }
    else if "`method'" == "global" {
        local method "joint"
    }
    if !inlist("`method'", "twostep", "joint") {
        di as error "method() must be {bf:twostep} (alias {bf:local}) or " ///
            "{bf:joint} (alias {bf:global})"
        exit 198
    }
    if "`grid_style'" == "" {
        local grid_style "default"
    }
    if "`joint_loocv'" == "" {
        // Default to exhaustive Cartesian search to match diff-diff==3.1.1
        // `_fit_global`, the numerical ground truth.  Users who need speed can
        // opt into cycling coordinate descent via joint_loocv(cycling).
        local joint_loocv "exhaustive"
    }
    if "`level'" == "" {
        local level = c(level)
    }

    // --- Validate joint_loocv() -------------------------------------------
    // Accept only "exhaustive" (Cartesian product, O(|grid|^3), matches
    // Python diff-diff 3.1.1 `_fit_global`, default) or "cycling" (coordinate
    // descent, O(|grid|*max_cycles), faster but may select a different lambda
    // because the Q(lambda) surface is non-convex).  The option is parsed
    // for any method() but only takes effect for method("joint") with LOOCV
    // enabled.
    if !inlist("`joint_loocv'", "cycling", "exhaustive") {
        di as error "joint_loocv() must be either {bf:exhaustive} or {bf:cycling}"
        di as error "  exhaustive (default): Cartesian product, O(|grid|^3),"
        di as error "                        matches Python diff-diff 3.1.1"
        di as error "  cycling             : coordinate descent, O(|grid|*cycles),"
        di as error "                        faster but may select a different lambda"
        exit 198
    }

    // bsalpha() is superseded by level().  When the sentinel value (-1)
    // indicates that bsalpha() was not specified, alpha is derived from
    // level().  Otherwise the user-supplied bsalpha() is honored with a
    // deprecation notice.
    if `bsalpha' != -1 {
        di as txt "{it:Note: bsalpha() is deprecated. Use level() instead.}"
        di as txt "{it:  bsalpha(`bsalpha') is equivalent to level(" %3.0f (1-`bsalpha')*100 ")}"
    }
    else {
        local bsalpha = 1 - `level'/100
    }

    // --- Estimation sample -------------------------------------------------
    marksample touse
    markout `touse' `panelvar' `timevar'

    // --- Header -----------------------------------------------------------
    di as txt _n "{hline 78}"
    di as txt "Triply Robust Panel Estimator (TROP)"
    di as txt "{hline 78}"

    // --- Load auxiliary ado-files -----------------------------------------
    // Each helper is loaded on demand: first from the PLUS sysdir, then the
    // current working directory, and finally via findfile along the adopath.
    foreach prog in _trop_load_plugin _trop_set_grid ///
                    _trop_validate_params ///
                    trop_handle_error trop_validate {
        capture program list `prog'
        if _rc {
            local _fl = substr("`prog'", 1, 1)
            capture quietly run "`c(sysdir_plus)'`_fl'/`prog'.ado"
            if _rc {
                capture quietly run "`c(pwd)'/`prog'.ado"
                if _rc {
                    capture quietly findfile `prog'.ado
                    if !_rc {
                        quietly run "`r(fn)'"
                    }
                }
            }
        }
    }

    // --- Plugin loading ---------------------------------------------------
    if "`verbose'" != "" {
        di as txt _n "Loading plugin..."
    }

    capture _trop_load_plugin
    if _rc {
        di as error "Error: TROP plugin not found. The compiled plugin is required."
        di as error ""
        di as error "Platform: `c(os)' `c(machine_type)'"
        di as error "Please install the TROP plugin for your platform."
        di as error "See {help trop##installation:trop installation} for details."
        exit 601
    }
    else {
        if "`verbose'" != "" {
            di as txt "  Platform: `_platform_desc'"
            di as txt "  Plugin: `_plugin_name'"
        }
    }
    
    // --- Lambda grid construction ------------------------------------------
    // Build the three-dimensional grid (lambda_time, lambda_unit, lambda_nn)
    // over which LOOCV minimizes Q(lambda).
    if "`verbose'" != "" {
        di as txt _n "Setting lambda grids..."
    }

    _trop_set_grid, grid_style(`grid_style') ///
        lambda_time_grid(`lambda_time_grid') ///
        lambda_unit_grid(`lambda_unit_grid') ///
        lambda_nn_grid(`lambda_nn_grid')

    local lambda_time_grid "`_lambda_time_grid'"
    local lambda_unit_grid "`_lambda_unit_grid'"
    local lambda_nn_grid "`_lambda_nn_grid'"
    local grid_style "`_grid_style'"
    local n_combinations = `_n_combinations'
    local n_per_cycle = `_n_per_cycle'

    if "`verbose'" != "" {
        di as txt "  Grid style: `_grid_style'"
        di as txt "  lambda_time: `_n_time' values"
        di as txt "  lambda_unit: `_n_unit' values"
        di as txt "  lambda_nn: `_n_nn' values"
        di as txt "  Grid points (Cartesian): `n_combinations'"
        di as txt "  Evaluations per cycle (coordinate descent): `n_per_cycle'"
    }

    // --- Parameter validation ---------------------------------------------
    if "`verbose'" != "" {
        di as txt _n "Validating parameters..."
    }

    _trop_validate_params, depvar(`depvar') treatvar(`treatvar') ///
        panelvar(`panelvar') timevar(`timevar') ///
        method(`method') grid_style(`grid_style') ///
        lambda_time_grid(`lambda_time_grid') ///
        lambda_unit_grid(`lambda_unit_grid') ///
        lambda_nn_grid(`lambda_nn_grid') ///
        bootstrap(`bootstrap') bsalpha(`bsalpha') ///
        tol(`tol') maxiter(`maxiter') ///
        seed(`seed') ///
        touse(`touse')

    if "`verbose'" != "" {
        di as txt "  All parameters valid"
    }

    // --- Mata function loading --------------------------------------------
    // Mata routines must be available before data validation, which calls
    // into Mata for panel structure checks.
    _trop_load_mata
    
    // --- Data validation ---------------------------------------------------
    if "`verbose'" != "" {
        di as txt _n "Validating data..."
    }

    // Drop temporary variables that may persist from a prior estimation run
    _trop_cleanup_vars

    if "`verbose'" != "" {
        capture noisily trop_validate `depvar' `treatvar' if `touse', ///
            panelvar(`panelvar') timevar(`timevar')
    }
    else {
        capture trop_validate `depvar' `treatvar' if `touse', ///
            panelvar(`panelvar') timevar(`timevar')
    }

    if _rc != 0 {
        if "`verbose'" == "" {
            // Re-run with output so the user can see what failed
            noisily trop_validate `depvar' `treatvar' if `touse', ///
                panelvar(`panelvar') timevar(`timevar')
        }
        di as error _n "Data validation failed."
        exit _rc
    }

    if e(data_validated) != 1 {
        di as error _n "Data validation did not pass (e(data_validated)=0)"
        exit 459
    }

    local N = e(N)
    local T_val = e(T)
    local N_obs = e(N_obs)
    local N_treat = e(N_treat)
    local N_control = `N_obs' - `N_treat'

    // Cache validation diagnostics in locals; ereturn post (below) clears
    // all e() scalars and macros.
    local balanced = e(balanced)
    local miss_rate = e(miss_rate)
    local N_treated_units = e(N_treated_units)
    local T_treat_periods = e(T_treat_periods)
    local treatment_pattern "`e(treatment_pattern)'"
    local N_control_units = e(N_control_units)
    local min_pre_treated = e(min_pre_treated)
    local min_valid_pairs = e(min_valid_pairs)
    local has_switching = e(has_switching)
    local max_switches = e(max_switches)
    local data_validated = e(data_validated)
    local time_min = e(time_min)
    local time_max = e(time_max)
    local time_range = e(time_range)
    local n_pre_periods = e(n_pre_periods)
    local n_post_periods = e(n_post_periods)

    if "`verbose'" != "" {
        di as txt "  Data validation passed"
        di as txt "  Units (N): `N'"
        di as txt "  Periods (T): `T_val'"
        di as txt "  Observations: `N_obs'"
        di as txt "  Treated: `N_treat' (" %4.1f 100*`N_treat'/`N_obs' "%)"
        di as txt "  Control: `N_control' (" %4.1f 100*`N_control'/`N_obs' "%)"
    }

    // --- Method-pattern compatibility check -------------------------------
    // The joint method solves a single weighted least-squares problem for a
    // scalar tau, which requires all treated units to share the same adoption
    // period.  Staggered adoption violates this assumption.
    if "`method'" == "joint" & "`treatment_pattern'" == "staggered_adoption" {
        di as error _n "{bf:ERROR: method('joint') requires simultaneous treatment adoption}"
        di as error "Your data shows staggered adoption (units first treated at different periods)."
        di as error ""
        di as error "The joint method estimates a single scalar tau with shared weights,"
        di as error "which assumes all treated units begin treatment at the same period."
        di as error "Staggered adoption violates this assumption and would produce"
        di as error "incorrect results."
        di as error ""
        di as error "Solution: Use method('twostep') which estimates individual treatment"
        di as error "effects per (unit, period) pair and properly handles staggered designs."
        di as error ""
        di as error "  trop `depvar' `treatvar', panelvar(`panelvar') timevar(`timevar') method(twostep)"
        exit 459
    }
    
    // --- Panel and time index variables ------------------------------------
    // Create consecutive integer indices for the plugin, which expects
    // 1-based contiguous identifiers.
    tempvar panel_idx time_idx
    qui egen `panel_idx' = group(`panelvar') if `touse'
    qui egen `time_idx' = group(`timevar') if `touse'
    sort `panel_idx' `time_idx'

    // --- Transfer lambda grids to Stata matrices -------------------------
    // The plugin reads grids from named matrices __trop_lambda_*_grid.
    // Boundary encoding (inf -> 0 for time/unit, inf -> 1e10 for nn) is
    // performed by the C bridge layer, not here.
    local n_time : word count `lambda_time_grid'
    tempname mat_time_tmp mat_unit_tmp mat_nn_tmp
    matrix `mat_time_tmp' = J(1, `n_time', .)
    local i = 1
    foreach val of local lambda_time_grid {
        matrix `mat_time_tmp'[1, `i'] = `val'
        local ++i
    }
    matrix __trop_lambda_time_grid = `mat_time_tmp'

    local n_unit : word count `lambda_unit_grid'
    matrix `mat_unit_tmp' = J(1, `n_unit', .)
    local i = 1
    foreach val of local lambda_unit_grid {
        matrix `mat_unit_tmp'[1, `i'] = `val'
        local ++i
    }
    matrix __trop_lambda_unit_grid = `mat_unit_tmp'

    local n_nn : word count `lambda_nn_grid'
    matrix `mat_nn_tmp' = J(1, `n_nn', .)
    local i = 1
    foreach val of local lambda_nn_grid {
        matrix `mat_nn_tmp'[1, `i'] = `val'
        local ++i
    }
    matrix __trop_lambda_nn_grid = `mat_nn_tmp'

    // --- Fixed lambda option ----------------------------------------------
    // fixedlambda(lambda_time lambda_unit lambda_nn) bypasses LOOCV and
    // uses the supplied triplet directly for estimation.
    local run_cv = 1
    local lambda_time_fixed = 0
    local lambda_unit_fixed = 0
    local lambda_nn_fixed = 0
    local lambda_time_val = .
    local lambda_unit_val = .
    local lambda_nn_val = .

    if "`fixedlambda'" != "" {
        local n_fixed : word count `fixedlambda'
        if `n_fixed' != 3 {
            di as error "fixedlambda() requires exactly 3 values: lambda_time lambda_unit lambda_nn"
            di as error "  Example: fixedlambda(0.5 1.0 0.1)"
            exit 198
        }
        local lambda_time_val : word 1 of `fixedlambda'
        local lambda_unit_val : word 2 of `fixedlambda'
        local lambda_nn_val : word 3 of `fixedlambda'

        // λ_time and λ_unit must be non-negative finite numbers.  Per paper
        // Eq 3 the weights θ,ω are exp(−λ·dist); λ_time/λ_unit = ∞ would
        // collapse all weight to the target period/unit only, which is
        // degenerate and not a supported estimator configuration.
        foreach lname in lambda_time_val lambda_unit_val {
            // Detect the literal "." (Stata missing) before `confirm number`,
            // which rejects "." even though numlist(missingokay) accepts it.
            if "``lname''" == "." {
                di as error "fixedlambda() does not accept `.' (missing / inf) for lambda_time or lambda_unit."
                di as error "  Paper Eq 3: θ_s = exp(−lambda_time·|t−s|), ω_j = exp(−lambda_unit·dist)."
                di as error "  lambda_time = lambda_unit = 0 recovers uniform weights; use finite values ≥ 0."
                exit 198
            }
            capture confirm number ``lname''
            if _rc {
                di as error "fixedlambda() values must be numeric: found '``lname''"
                exit 198
            }
            if ``lname'' < 0 {
                di as error "fixedlambda() values must be non-negative"
                exit 198
            }
        }

        // λ_nn = missing / inf is the DID/TWFE special case (L ≡ 0).
        // Map to the large-finite sentinel 1e10 recognized by the C bridge.
        if "`lambda_nn_val'" == "." {
            local lambda_nn_val = 1e10
            if "`verbose'" != "" {
                di as txt "  lambda_nn = . (inf) → mapped to 1e10 (DID/TWFE special case: L ≡ 0)"
            }
        }
        else {
            capture confirm number `lambda_nn_val'
            if _rc {
                di as error "fixedlambda() values must be numeric: found '`lambda_nn_val''"
                exit 198
            }
            if `lambda_nn_val' < 0 {
                di as error "fixedlambda() lambda_nn must be non-negative (use . for +inf)"
                exit 198
            }
        }

        local run_cv = 0
        local lambda_time_fixed = 1
        local lambda_unit_fixed = 1
        local lambda_nn_fixed = 1

        if "`verbose'" != "" {
            di as txt "  Fixed lambda (LOOCV skipped):"
            di as txt "    lambda_time = `lambda_time_val'"
            di as txt "    lambda_unit = `lambda_unit_val'"
            di as txt "    lambda_nn   = `lambda_nn_val'"
        }
    }
    
    // --- Estimation --------------------------------------------------------
    if "`verbose'" != "" {
        di as txt _n "Running TROP estimation..."
        di as txt "  Method: `method'"
        if `run_cv' {
            di as txt "  LOOCV: enabled"
        }
        else {
            di as txt "  LOOCV: disabled (using fixed lambda)"
        }
        if `bootstrap' > 0 {
            di as txt "  Bootstrap: `bootstrap' replications (paper Alg 3)"
        }
        else {
            di as txt "  Bootstrap: skipped (user requested bootstrap(0))"
        }
    }

    local verbose_flag = ("`verbose'" != "")
    local cv_mode = "exact"

    // The seed is forwarded to the plugin's internal RNG; Stata's global
    // RNG state is left unchanged.

    // Transfer joint_loocv() mode to the Mata layer via a global macro.
    // The Mata dispatcher _trop_main reads __trop_joint_loocv_mode and
    // routes to either trop_loocv_joint() or trop_loocv_joint_exhaustive()
    // when method == "joint" and LOOCV is enabled.
    //
    // The global is set via Mata's st_global() because Stata's `global`
    // command rejects names with a leading underscore (r(198)), whereas
    // Mata's st_global() accepts any identifier.
    mata: st_global("__trop_joint_loocv_mode", "`joint_loocv'")

    // Dispatch to the compiled plugin through the Mata interface layer.
    // trop_main() returns 0 on success or a Stata return code on failure.
    mata: st_local("_trop_rc", strofreal(trop_main("`depvar'", "`treatvar'", "`panel_idx'", ///
        "`time_idx'", "`touse'", "`method'", ///
        `lambda_time_val', `lambda_unit_val', `lambda_nn_val', ///
        `run_cv', `bootstrap', `seed', `verbose_flag', ///
        `maxiter', `tol', `bsalpha')))

    // Clear the global once the Mata call returns so that stale state does
    // not leak into the next estimation.  Setting it to "" via st_global is
    // equivalent to dropping it from the Mata perspective.
    mata: st_global("__trop_joint_loocv_mode", "")

    // --- Post e(b) and e(V) ----------------------------------------------
    // ereturn post clears all prior e() contents, so it must precede the
    // storage of any other e() scalars or macros.

    if `_trop_rc' != 0 {
        di as error "Estimation failed (error code `_trop_rc')"
        exit `_trop_rc'
    }

    // Retrieve point estimate and standard error from plugin temporaries
    local _att_val = .
    capture local _att_val = scalar(__trop_att)
    local _se_val = .
    capture local _se_val = scalar(__trop_se)

    if missing(`_att_val') {
        di as error "Estimation failed: ATT is missing"
        exit 498
    }

    // Build b (1x1) and, when available, V (1x1).  No closed-form
    // asymptotic variance exists for the TROP estimator; V is populated
    // only when bootstrap inference has been performed.
    tempname _b
    matrix `_b' = (`_att_val')
    matrix colnames `_b' = att
    // `ereturn post … esample(var)` consumes `var` (it is dropped from the
    // data and preserved only as e(sample)).  Downstream Mata code still
    // needs to read the touse variable, so clone it into a fresh tempvar
    // before posting.  The clone is dropped automatically at program exit.
    tempvar _touse_esample
    qui gen byte `_touse_esample' = `touse'
    if `_se_val' > 0 & !missing(`_se_val') {
        tempname _V
        matrix `_V' = (`_se_val' ^ 2)
        matrix colnames `_V' = att
        matrix rownames `_V' = att
        ereturn post `_b' `_V', esample(`_touse_esample')
    }
    else {
        ereturn post `_b', esample(`_touse_esample')
    }

    // Expose the tempvar names so trop_store_results can build the
    // (time x unit) tau matrix by re-reading treatment coordinates.
    // Global macros are cleared automatically by trop_cleanup_temp_vars.
    mata: st_global("__trop_treatvar", "`treatvar'")
    mata: st_global("__trop_panel_idx_var", "`panel_idx'")
    mata: st_global("__trop_time_idx_var", "`time_idx'")

    // Transfer remaining estimation results from plugin temporaries to e()
    mata: trop_store_results("`method'")
    
    // --- Deferred error flag ------------------------------------------------
    local _fatal_rc = 0

    // --- LOOCV diagnostics ------------------------------------------------
    capture confirm scalar __trop_loocv_n_attempted
    if !_rc {
        mata: store_loocv_diagnostics("`method'", `seed')

        if "`verbose'" != "" {
            mata: display_loocv_verbose()
        }

        // A high LOOCV failure rate signals numerical problems; defer the
        // error exit until all e() results have been stored.
        mata: st_local("loocv_rc", strofreal(check_loocv_fail_rate()))
        if `loocv_rc' != 0 {
            local _fatal_rc = `loocv_rc'
        }
    }

    // --- Capture deferred fatal errors ------------------------------------
    // Bootstrap or LOOCV may flag a fatal condition via __trop_fatal_rc;
    // read it before temporary scalars are dropped.
    if `_fatal_rc' == 0 {
        capture confirm scalar __trop_fatal_rc
        if !_rc {
            local _fatal_rc = scalar(__trop_fatal_rc)
        }
    }

    // --- Store lambda grids in e() ----------------------------------------
    capture ereturn matrix lambda_time_grid = __trop_lambda_time_grid
    capture ereturn matrix lambda_unit_grid = __trop_lambda_unit_grid
    capture ereturn matrix lambda_nn_grid = __trop_lambda_nn_grid

    // --- Drop plugin temporaries ------------------------------------------
    // All code that reads __trop_* scalars or matrices must appear above.
    capture scalar drop __trop_att __trop_se
    capture mata: trop_cleanup_temp_vars()

    // --- Store estimation metadata ----------------------------------------
    ereturn local method "`method'"
    ereturn local grid_style "`grid_style'"
    ereturn local joint_loocv "`joint_loocv'"
    ereturn local cmd "trop"
    ereturn local estat_cmd "trop_estat"
    ereturn local cmdline "trop `0'"
    ereturn local depvar "`depvar'"
    ereturn local treatvar "`treatvar'"
    ereturn local panelvar "`panelvar'"
    ereturn local timevar "`timevar'"

    ereturn scalar N_units = `N'
    ereturn scalar N_periods = `T_val'
    ereturn scalar bootstrap_reps = `bootstrap'
    ereturn scalar alpha_level = `bsalpha'
    ereturn scalar loocv_used = `run_cv'
    if `bootstrap' > 0 {
        ereturn local vcetype "Bootstrap"
    }
    else {
        ereturn local vcetype ""
    }

    // LOOCV configuration scalars
    ereturn scalar seed = `seed'

    // Lambda grid dimensions
    ereturn scalar n_lambda_time = `_n_time'
    ereturn scalar n_lambda_unit = `_n_unit'
    ereturn scalar n_lambda_nn = `_n_nn'
    ereturn scalar n_grid_combinations = `n_combinations'
    // Coordinate-descent LOOCV evaluates n_time + n_unit + n_nn grid points
    // per cycle (one univariate sweep per regularization parameter).
    ereturn scalar n_grid_per_cycle = `n_per_cycle'

    // --- Restore validation diagnostics -----------------------------------
    // These were cached in locals before ereturn post cleared e().
    ereturn scalar N_obs = `N_obs'
    ereturn scalar balanced = `balanced'
    ereturn scalar miss_rate = `miss_rate'
    ereturn scalar N_treat = `N_treat'
    ereturn scalar N_control = `N_control'
    ereturn scalar N_treated_units = `N_treated_units'
    ereturn scalar T_treat_periods = `T_treat_periods'
    ereturn local treatment_pattern "`treatment_pattern'"
    ereturn scalar N_control_units = `N_control_units'
    ereturn scalar min_pre_treated = `min_pre_treated'
    ereturn scalar min_valid_pairs = `min_valid_pairs'
    ereturn scalar has_switching = `has_switching'
    ereturn scalar max_switches = `max_switches'
    if !missing(`n_pre_periods') {
        ereturn scalar n_pre_periods = `n_pre_periods'
    }
    if !missing(`n_post_periods') {
        ereturn scalar n_post_periods = `n_post_periods'
    }
    if !missing(`data_validated') {
        ereturn scalar data_validated = `data_validated'
    }
    if !missing(`time_min') {
        ereturn scalar time_min = `time_min'
    }
    if !missing(`time_max') {
        ereturn scalar time_max = `time_max'
    }
    if !missing(`time_range') {
        ereturn scalar time_range = `time_range'
    }

    // --- Cleanup temporary dataset variables ------------------------------
    _trop_cleanup_vars

    // --- Data integrity for predict consistency ----------------------------
    // The _depvar_checksum (below) provides a lightweight mechanism for
    // detecting whether the data has changed since estimation.  A full
    // datasignature is not stored because tempvars still exist in-scope
    // when this code runs, producing a false-positive mismatch after the
    // program returns and Stata drops those tempvars.

    // Lightweight checksum on the dependent variable
    capture {
        qui sum `depvar' if e(sample), meanonly
        ereturn scalar _depvar_checksum = r(sum)
    }

    // --- Deferred fatal error exit ----------------------------------------
    // Issued after all e() storage so that captured runs retain partial
    // results for diagnostic inspection.
    if `_fatal_rc' != 0 {
        exit `_fatal_rc'
    }
    
    // --- Display results ----------------------------------------------------
    di as txt _n "{hline 78}"
    di as txt "TROP Estimation Results"
    di as txt "{hline 78}"

    di as txt "Method:" _col(20) as res "`method'"
    di as txt "Grid style:" _col(20) as res "`grid_style'" as txt " (`n_per_cycle' grid points/cycle, coordinate descent)"
    if "`method'" == "joint" & `run_cv' {
        if "`joint_loocv'" == "exhaustive" {
            di as txt "Joint LOOCV:" _col(20) as res "exhaustive" as txt " (Cartesian product; default, matches Python 3.1.1)"
        }
        else {
            di as txt "Joint LOOCV:" _col(20) as res "cycling" as txt " (coordinate descent; faster but may differ from Python)"
        }
    }
    di as txt ""
    di as txt "Panel dimensions:" _col(20) as res "N = `N', T = `T_val'"
    di as txt "Observations:" _col(20) as res "`N_obs'"
    di as txt "Treated:" _col(20) as res "`N_treat'" as txt " (" %4.1f 100*`N_treat'/`N_obs' "%)"

    // Selected or fixed hyperparameters
    if !missing(e(lambda_time)) {
        if `run_cv' {
            di as txt _n "Selected hyperparameters (via LOOCV):"
        }
        else {
            di as txt _n "Fixed hyperparameters (LOOCV skipped):"
        }
        di as txt "  lambda_time = " as res %8.4f e(lambda_time)
        di as txt "  lambda_unit = " as res %8.4f e(lambda_unit)
        di as txt "  lambda_nn   = " as res %8.4f e(lambda_nn)

        if !missing(e(loocv_score)) {
            di as txt "  Q(lambda_hat) = " as res %10.6f e(loocv_score)
        }
    }

    // Point estimate and inference
    di as txt _n "Treatment Effect (ATT):"

    local att = e(att)
    local se = e(se)
    local ci_lower = e(ci_lower)
    local ci_upper = e(ci_upper)
    local pvalue = e(pvalue)
    local tstat = e(t)

    // Percentile CI from the bootstrap distribution (paper Algorithm 3).
    // Displayed alongside the t-based CI so users can compare parametric and
    // distribution-free intervals and detect bootstrap asymmetry.  The
    // t-based CI is retained as the primary one (e(ci_lower), e(ci_upper))
    // to preserve backward compatibility with downstream commands such as
    // `test` and `lincom`.
    local ci_lower_pct = e(ci_lower_percentile)
    local ci_upper_pct = e(ci_upper_percentile)

    di as txt "  tau     = " as res %12.6f `att'
    if `bootstrap' > 0 {
        di as txt "  SE     = " as res %12.6f `se'
        di as txt "  t      = " as res %12.4f `tstat'
        di as txt "  p-value= " as res %12.4f `pvalue'
        di as txt "  `level'% CI (t-based)   " _col(27) "= [" ///
            as res %12.6f `ci_lower' as txt ", " ///
            as res %12.6f `ci_upper' as txt "]"
        if !missing(`ci_lower_pct') & !missing(`ci_upper_pct') {
            di as txt "  `level'% CI (percentile)" _col(27) "= [" ///
                as res %12.6f `ci_lower_pct' as txt ", " ///
                as res %12.6f `ci_upper_pct' as txt "]" ///
                as txt "  {it:(paper Alg 3)}"
        }
    }
    else {
        di as txt "  SE     = " as res "(not computed)" ///
            as txt "  {it:use {bf:bootstrap()} option}"
        di as txt "  t      = " as res "(not computed)"
        di as txt "  p-value= " as res "(not computed)"
        di as txt "  `level'% CI = " as res "(not computed)"
    }

    // Global intercept (joint method only)
    if "`method'" == "joint" & !missing(e(mu)) {
        di as txt _n "Global intercept:"
        di as txt "  mu     = " as res %12.6f e(mu)
    }

    // Convergence diagnostics
    if !missing(e(converged)) {
        di as txt _n "Convergence:"
        local conv_status = cond(e(converged)==1, "Yes", "No")
        di as txt "  Converged: " as res "`conv_status'"
        if !missing(e(n_iterations)) {
            // twostep: maximum across per-observation fits;
            // joint: total model iterations.
            if "`method'" == "twostep" {
                di as txt "  Iterations (max per-obs): " as res e(n_iterations)
            }
            else {
                di as txt "  Iterations: " as res e(n_iterations)
            }
        }
        if "`method'" == "twostep" {
            capture confirm scalar e(n_obs_estimated)
            if !_rc & !missing(e(n_obs_estimated)) {
                di as txt "  Obs estimated: " as res e(n_obs_estimated)
                capture confirm scalar e(n_obs_failed)
                if !_rc & !missing(e(n_obs_failed)) & e(n_obs_failed) > 0 {
                    di as txt "  Obs failed:    " as res e(n_obs_failed)
                }
            }
        }
    }

    // Bootstrap summary
    if `bootstrap' > 0 {
        di as txt _n "Bootstrap inference (paper Alg 3):"
        di as txt "  Replications: " as res `bootstrap'
        di as txt "  Alpha level: " as res `bsalpha'
    }
    else {
        di as txt _n "{it:Note: Bootstrap inference skipped (bootstrap(0)).}"
        di as txt "{it:  Standard errors per Athey et al. (2025) Algorithm 3 require bootstrap.}"
        di as txt "{it:  To enable inference, re-run without bootstrap(0) (default is 200 reps),}"
        di as txt "{it:  or call: trop_bootstrap, nreps(200)}"
    }

    di as txt "{hline 78}"
end
