*! estat weights subcommand for trop

/*==============================================================================
  trop_estat_weights

  Display descriptive statistics for TROP weights (Time and Unit).

  Syntax:
    estat weights [, heatmap detailed]

  Options:
    heatmap  -  Display a heatmap of the weight matrix.
    detailed -  Display detailed percentiles of the weight vectors.

  Description:
    Reports summary statistics for the estimated weights.

    For the Two-Step estimator, weights are observation-specific.
    This command reports statistics for the weights associated with the
    first treated observation. Two-Step weights consist of:
      theta (T x 1): Time weights
      omega (N x 1): Unit weights

    For the Joint estimator, weights are global across observations:
      delta_time (T x 1): Global time weights
      delta_unit (N x 1): Global unit weights
==============================================================================*/

program define trop_estat_weights
    version 17
    syntax [, heatmap detailed]
    
    // Check for trop estimation results
    if "`e(cmd)'" != "trop" {
        di as error "last estimates not found"
        exit 301
    }
    
    // Identify estimation method
    local method = "`e(method)'"
    if "`method'" == "" {
        local method = "twostep"
    }
    
    // Execute method-specific display
    if "`method'" == "twostep" {
        _estat_weights_twostep
    }
    else if "`method'" == "joint" {
        _estat_weights_joint
    }
    else {
        _estat_weights_generic
    }
    
    // Display heatmap if requested
    if "`heatmap'" != "" {
        _graph_weight_heatmap
    }
    
    // Display detailed percentiles if requested
    if "`detailed'" != "" {
        _estat_weights_detailed
    }
end

/*==============================================================================
  Display weights for Two-Step estimator

  Displays summary statistics (Mean, Std. Dev., Min, Max) for theta and omega.
  As weights vary by treated observation, statistics are reported for the
  first treated observation.
==============================================================================*/
program define _estat_weights_twostep
    di as txt ""
    di as txt "Weight Distribution (Twostep; Algorithm 2 default)"
    di as txt "{hline 61}"
    di as txt "Note: Weights vary by treated observation. Showing statistics"
    di as txt "      for first treated observation."
    di as txt ""
    
    // Process time weights (theta)
    capture confirm matrix e(theta)
    local has_theta = (_rc == 0)
    
    capture confirm matrix e(omega)
    local has_omega = (_rc == 0)
    
    if `has_theta' {
        tempname theta
        matrix `theta' = e(theta)
        local T = rowsof(`theta')
        
        di as txt "Time weights (theta):"
        mata: _display_weight_vector_stats("e(theta)", "time")
    }
    else {
        di as txt "Time weights (theta): (not available)"
    }
    
    di as txt ""
    
    // Process unit weights (omega)
    if `has_omega' {
        tempname omega
        matrix `omega' = e(omega)
        local N = rowsof(`omega')
        
        di as txt "Unit weights (omega):"
        mata: _display_weight_vector_stats("e(omega)", "unit")
    }
    else {
        di as txt "Unit weights (omega): (not available)"
    }
    
    di as txt "{hline 61}"
end

/*==============================================================================
  Display weights for Joint estimator

  Displays summary statistics for delta_time and delta_unit.
  Weights are global in the Joint estimator.
==============================================================================*/
program define _estat_weights_joint
    di as txt ""
    di as txt "Weight Distribution (Joint; shared-tau extension)"
    di as txt "{hline 61}"
    
    capture confirm matrix e(delta_time)
    local has_delta_time = (_rc == 0)
    
    capture confirm matrix e(delta_unit)
    local has_delta_unit = (_rc == 0)
    
    // Global time weights
    di as txt "Global time weights (delta_time):"
    if `has_delta_time' {
        mata: _display_weight_vector_stats("e(delta_time)", "time")
    }
    else {
        di as txt "  (not available)"
    }
    
    di as txt ""
    
    // Global unit weights
    di as txt "Global unit weights (delta_unit):"
    if `has_delta_unit' {
        mata: _display_weight_vector_stats("e(delta_unit)", "unit")
    }
    else {
        di as txt "  (not available)"
    }
    
    di as txt "{hline 61}"
end

/*==============================================================================
  Generic weight display

  Displays summary statistics for available weight vectors when the
  estimation method is not explicitly identified.
==============================================================================*/
program define _estat_weights_generic
    di as txt ""
    di as txt "Weight Distribution"
    di as txt "{hline 61}"
    
    capture confirm matrix e(theta)
    if !_rc {
        di as txt "Time weights:"
        mata: _display_weight_vector_stats("e(theta)", "time")
        di as txt ""
    }
    
    capture confirm matrix e(omega)
    if !_rc {
        di as txt "Unit weights:"
        mata: _display_weight_vector_stats("e(omega)", "unit")
    }
    
    di as txt "{hline 61}"
end

/*==============================================================================
  Weight Heatmap

  Generates a heatmap of the weight matrix.
  If e(W_mat) is available, it is used directly.
  Otherwise, an outer-product approximation is constructed:
    - Twostep: theta * omega'
    - Joint: delta_time * delta_unit'
  Requires the user-written command -heatplot-.
==============================================================================*/
program define _graph_weight_heatmap
    // Check for heatplot package
    capture which heatplot
    if _rc {
        di as error "heatplot command not found"
        di as error "Install via: {bf:ssc install heatplot}"
        di as error ""
        di as error "Alternative: Use {bf:estat weights} without {bf:heatmap}"
        exit 199
    }
    
    // Check for stored weight matrix
    capture confirm matrix e(W_mat)
    local has_W_mat = (_rc == 0)
    
    if `has_W_mat' {
        heatplot e(W_mat), ///
            title("Weight Matrix") ///
            xlabel("Time period") ///
            ylabel("Unit") ///
            legend(title("Weight")) ///
            scheme(s2color)
    }
    else {
        // Construct weight matrix from available vectors
        local can_reconstruct = 0
        capture confirm matrix e(theta)
        if !_rc {
            capture confirm matrix e(omega)
            if !_rc {
                local can_reconstruct = 1
                tempname W_approx
                mata: st_matrix("`W_approx'", st_matrix("e(theta)") * st_matrix("e(omega)")')
            }
        }
        if !`can_reconstruct' {
            capture confirm matrix e(delta_time)
            if !_rc {
                capture confirm matrix e(delta_unit)
                if !_rc {
                    local can_reconstruct = 1
                    tempname W_approx
                    mata: st_matrix("`W_approx'", st_matrix("e(delta_time)") * st_matrix("e(delta_unit)")')
                }
            }
        }
        if `can_reconstruct' {
            heatplot `W_approx', ///
                title("Weight Matrix (reconstructed)") ///
                xlabel("Time period") ///
                ylabel("Unit") ///
                legend(title("Weight")) ///
                note("Note: Outer product approximation.") ///
                scheme(s2color)
        }
        else {
            di as error "Cannot generate heatmap: weight matrices not available"
            exit 111
        }
    }
end

/*==============================================================================
  Detailed Weight Distribution

  Displays percentiles (10, 25, 50, 75, 90) for weight vectors.
==============================================================================*/
program define _estat_weights_detailed
    di as txt ""
    di as txt "Detailed weight distribution (percentiles):"
    di as txt ""
    
    // Time weights
    capture confirm matrix e(theta)
    if !_rc {
        di as txt "Time weights percentiles:"
        mata: _display_weight_percentiles("e(theta)")
        di as txt ""
    }
    else {
        capture confirm matrix e(delta_time)
        if !_rc {
            di as txt "Time weights percentiles:"
            mata: _display_weight_percentiles("e(delta_time)")
            di as txt ""
        }
    }
    
    // Unit weights
    capture confirm matrix e(omega)
    if !_rc {
        di as txt "Unit weights percentiles:"
        mata: _display_weight_percentiles("e(omega)")
    }
    else {
        capture confirm matrix e(delta_unit)
        if !_rc {
            di as txt "Unit weights percentiles:"
            mata: _display_weight_percentiles("e(delta_unit)")
        }
    }
end

/*==============================================================================
  Mata Helper Functions
==============================================================================*/

version 17
mata:
mata set matastrict on

/*------------------------------------------------------------------------------
  _display_weight_vector_stats()

  Displays summary statistics (Mean, Std. Dev., Min, Max) for a weight vector.
  Also reports the index (time or unit) of the minimum and maximum values.

  Arguments:
    matname  - Name of the Stata matrix (e.g., "e(theta)")
    wtype    - "time" or "unit" (determines index labeling)
------------------------------------------------------------------------------*/
void _display_weight_vector_stats(string scalar matname, string scalar wtype)
{
    real colvector w
    real scalar n, min_val, max_val, mean_val, sd_val
    real scalar min_idx, max_idx, i
    
    w = st_matrix(matname)
    
    // Ensure column orientation
    if (cols(w) > rows(w)) {
        w = w'
    }
    
    n = rows(w)
    min_val = min(w)
    max_val = max(w)
    mean_val = mean(w)
    sd_val = sqrt(variance(w))
    
    // Find index of first occurrence of min and max
    min_idx = .
    max_idx = .
    for (i=1; i<=n; i++) {
        if (w[i] == min_val & min_idx == .) min_idx = i
        if (w[i] == max_val & max_idx == .) max_idx = i
    }
    
    printf("{txt}  Mean     = {res}%10.4f\n", mean_val)
    printf("{txt}  Std.Dev. = {res}%10.4f\n", sd_val)
    
    if (wtype == "time") {
        printf("{txt}  Min      = {res}%10.4f{txt} (t={res}%g{txt})\n", min_val, min_idx)
        printf("{txt}  Max      = {res}%10.4f{txt} (t={res}%g{txt})\n", max_val, max_idx)
    }
    else {
        printf("{txt}  Min      = {res}%10.4f{txt} (i={res}%g{txt})\n", min_val, min_idx)
        printf("{txt}  Max      = {res}%10.4f{txt} (i={res}%g{txt})\n", max_val, max_idx)
    }
}

/*------------------------------------------------------------------------------
  _trop_interpolate_percentile()

  Calculates the p-th percentile of a sorted vector using linear interpolation.
  Index calculation: (n-1)*p.

  Arguments:
    sorted_v - Sorted column vector
    p        - Percentile (0 to 1)

  Returns:
    Interpolated value
------------------------------------------------------------------------------*/
real scalar _trop_interpolate_percentile(real colvector sorted_v, real scalar p)
{
    real scalar n, idx_f, idx_low, idx_high, frac
    
    n = rows(sorted_v)
    if (n == 0) return(.)
    if (n == 1) return(sorted_v[1])
    
    idx_f = (n - 1) * p
    idx_low = floor(idx_f)
    idx_high = ceil(idx_f)
    
    idx_low = max((0, min((n - 1, idx_low))))
    idx_high = max((0, min((n - 1, idx_high))))
    
    if (idx_low == idx_high) {
        return(sorted_v[idx_low + 1])
    }
    else {
        frac = idx_f - idx_low
        return(sorted_v[idx_low + 1] * (1 - frac) + sorted_v[idx_high + 1] * frac)
    }
}

/*------------------------------------------------------------------------------
  _display_weight_percentiles()

  Displays the 10th, 25th, 50th, 75th, and 90th percentiles of a weight vector.

  Arguments:
    matname - Name of the Stata matrix
------------------------------------------------------------------------------*/
void _display_weight_percentiles(string scalar matname)
{
    real colvector w, sorted_w
    real scalar n, p10, p25, p50, p75, p90
    
    w = st_matrix(matname)
    
    // Ensure column orientation
    if (cols(w) > rows(w)) {
        w = w'
    }
    
    n = rows(w)
    sorted_w = sort(w, 1)
    
    p10 = _trop_interpolate_percentile(sorted_w, 0.10)
    p25 = _trop_interpolate_percentile(sorted_w, 0.25)
    p50 = _trop_interpolate_percentile(sorted_w, 0.50)
    p75 = _trop_interpolate_percentile(sorted_w, 0.75)
    p90 = _trop_interpolate_percentile(sorted_w, 0.90)
    
    printf("{txt}  10%%      = {res}%10.4f\n", p10)
    printf("{txt}  25%%      = {res}%10.4f\n", p25)
    printf("{txt}  50%%      = {res}%10.4f\n", p50)
    printf("{txt}  75%%      = {res}%10.4f\n", p75)
    printf("{txt}  90%%      = {res}%10.4f\n", p90)
}

end
