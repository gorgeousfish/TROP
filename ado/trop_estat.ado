*! trop_estat
*! Dispatcher for trop postestimation subcommands
*
*  Syntax:
*      estat <subcommand> [, options]
*
*  Supported subcommands:
*      summarize    - Summary of sample structure and treatment allocation
*      vce          - Variance-covariance matrix of estimators
*      sensitivity  - Hyperparameter sensitivity analysis
*      weights      - Weight diagnostics
*      bootstrap    - Bootstrap distribution diagnostics
*      loocv        - LOOCV hyperparameter selection diagnostics
*      factors      - Factor matrix (L) analysis
*      triplerob    - Triple-robustness bias decomposition (paper Theorem 5.1)
*
*  Abbreviations:
*      sum    -> summarize
*      sens   -> sensitivity
*      weight -> weights
*      boot   -> bootstrap
*      triple -> triplerob

program define trop_estat
    version 17

    // Parse the leading subcommand token
    gettoken subcmd rest : 0, parse(" ,")

    // Retain the comma in the remaining string to ensure downstream programs receive a valid option string
    local rest = trim("`rest'")

    // Verify that active estimation results correspond to trop
    if "`e(cmd)'" != "trop" {
        di as error "last estimates not found or not from trop"
        di as error "Run {bf:trop} command before using {bf:estat}"
        exit 301
    }

    // Route execution to the specified subcommand
    if "`subcmd'" == "summarize" | "`subcmd'" == "sum" {
        trop_estat_summarize `rest'
    }
    else if "`subcmd'" == "vce" {
        trop_estat_vce `rest'
    }
    else if "`subcmd'" == "sensitivity" | "`subcmd'" == "sens" {
        trop_estat_sensitivity `rest'
    }
    else if "`subcmd'" == "weights" | "`subcmd'" == "weight" {
        trop_estat_weights `rest'
    }
    else if "`subcmd'" == "bootstrap" | "`subcmd'" == "boot" {
        trop_estat_bootstrap `rest'
    }
    else if "`subcmd'" == "loocv" {
        trop_estat_loocv `rest'
    }
    else if "`subcmd'" == "factors" {
        trop_estat_factors `rest'
    }
    else if "`subcmd'" == "triplerob" | "`subcmd'" == "triple" {
        trop_estat_triplerob `rest'
    }
    else {
        di as error "estat subcommand '{bf:`subcmd'}' not recognized"
        di as error ""
        di as error "Valid subcommands:"
        di as error "  {bf:summarize}   - Sample and treatment structure summary"
        di as error "  {bf:vce}         - Variance-covariance matrix"
        di as error "  {bf:sensitivity} - Hyperparameter sensitivity analysis"
        di as error "  {bf:weights}     - Weight diagnostics"
        di as error "  {bf:bootstrap}   - Bootstrap distribution diagnostics"
        di as error "  {bf:loocv}       - LOOCV hyperparameter selection diagnostics"
        di as error "  {bf:factors}     - Factor matrix (L) SVD analysis"
        di as error "  {bf:triplerob}   - Triple-robustness bias decomposition (Theorem 5.1)"
        exit 199
    }
end
