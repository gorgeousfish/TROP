
/*
  trop_estat_loocv -- Display LOOCV hyperparameter selection diagnostics.

  Reports the regularization parameters selected via leave-one-out
  cross-validation, the LOOCV objective function value, the proportion
  of valid leave-one-out fits, the grid search style, and, when
  applicable, the coordinates of the first failed observation.
*/


program define trop_estat_loocv
    version 17
    
    // --- Verify that trop estimation results are in memory ---
    if "`e(cmd)'" != "trop" {
        di as error "last estimates not found"
        exit 301
    }
    
    // --- Header ---
    di as txt ""
    di as txt "LOOCV Diagnostics"
    di as txt "{hline 61}"
    
    // --- Selected regularization parameters ---
    di as txt "Selected hyperparameters:"
    
    capture confirm scalar e(lambda_time)
    if !_rc {
        di as txt "  lambda_time     = " as res %10.3f e(lambda_time)
    }
    else {
        di as txt "  lambda_time     = " as res "(not available)"
    }
    
    capture confirm scalar e(lambda_unit)
    if !_rc {
        di as txt "  lambda_unit     = " as res %10.3f e(lambda_unit)
    }
    else {
        di as txt "  lambda_unit     = " as res "(not available)"
    }
    
    capture confirm scalar e(lambda_nn)
    if !_rc {
        di as txt "  lambda_nn       = " as res %10.3f e(lambda_nn)
    }
    else {
        di as txt "  lambda_nn       = " as res "(not available)"
    }
    
    di as txt ""
    
    // --- LOOCV performance summary ---
    di as txt "LOOCV performance:"
    
    // Minimized objective function value Q(lambda_hat)
    capture confirm scalar e(loocv_score)
    if !_rc {
        di as txt "  Objective Q(" as res "λ̂" as txt ")  = " as res %10.4f e(loocv_score)
    }
    else {
        di as txt "  Objective Q(" as res "λ̂" as txt ")  = " as res "(not available)"
    }
    
    // Valid fits as a fraction of total leave-one-out iterations
    capture confirm scalar e(loocv_n_valid)
    local has_valid = !_rc
    capture confirm scalar e(loocv_n_attempted)
    local has_attempted = !_rc
    
    if `has_valid' & `has_attempted' {
        local n_valid = e(loocv_n_valid)
        local n_attempted = e(loocv_n_attempted)
        
        if `n_attempted' > 0 {
            local pct_valid = 100 * `n_valid' / `n_attempted'
            di as txt "  Valid fits      = " as res %5.0f `n_valid' ///
               as txt " / " as res %5.0f `n_attempted' ///
               as txt " (" as res %5.1f `pct_valid' as txt "%)"
            
            // Warn if more than 10% of leave-one-out fits failed
            if `pct_valid' < 90 {
                di as err "  Warning: High failure rate (>" as res "10%" ///
                   as err ") may indicate data quality issues"
            }
        }
        else {
            di as txt "  Valid fits      = " as res "(no attempts recorded)"
        }
    }
    else {
        di as txt "  Valid fits      = " as res "(not available)"
    }
    
    // Grid search style (e.g., "auto", "manual")
    if "`e(grid_style)'" != "" {
        di as txt "  Grid style      = " as res "`e(grid_style)'"
    }
    else {
        di as txt "  Grid style      = " as res "(not available)"
    }
    
    // --- First failed leave-one-out observation, if any ---
    capture confirm scalar e(loocv_first_failed_t)
    local has_failed_t = !_rc
    capture confirm scalar e(loocv_first_failed_i)
    local has_failed_i = !_rc
    
    if `has_failed_t' & `has_failed_i' {
        local first_t = e(loocv_first_failed_t)
        local first_i = e(loocv_first_failed_i)
        
        // Non-negative indices indicate at least one leave-one-out fit failed
        if `first_t' >= 0 & `first_i' >= 0 {
            di as txt ""
            di as txt "First failed observation:"
            di as txt "  Time index      = " as res %5.0f `first_t'
            di as txt "  Unit index      = " as res %5.0f `first_i'
        }
    }
    
    di as txt "{hline 61}"
end
