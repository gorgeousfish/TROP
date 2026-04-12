/*──────────────────────────────────────────────────────────────────────────────
  trop_estat_helpers.mata

  Helper functions for estat subcommands: weight diagnostics and SVD
  analysis of the low-rank factor matrix.

  The TROP estimator assigns exponential-decay weights to both units and
  time periods.  These helpers quantify the effective concentration of
  those weights via Shannon entropy, Kish's effective sample size, and a
  top-k concentration index.  A separate routine performs SVD on the
  estimated factor matrix mu to report its effective rank, singular-value
  spectrum, and matrix norms (Frobenius and nuclear).

  Contents
    _compute_entropy()            Shannon entropy of a weight vector
    _compute_ess()                effective sample size (Kish, 1965)
    _compute_concentration()      top-k cumulative weight fraction
    compute_weight_stats()        aggregate weight diagnostics
    compute_bootstrap_stats()     bootstrap distribution summary
    _trop_estat_factors_svd()     SVD of the factor matrix
──────────────────────────────────────────────────────────────────────────────*/

version 17
mata:
mata set matastrict on

/*──────────────────────────────────────────────────────────────────────────────
  _compute_entropy()

  Shannon entropy of a normalised weight vector:

      H(w) = -sum_{i: w_i > 0} w_i * ln(w_i)

  with the convention 0 * ln(0) = 0.

  Bounds:  H = ln(n) when w_i = 1/n for all i  (uniform);
           H = 0     when a single w_k = 1      (degenerate).

  Arguments
    w   N x 1 column vector with sum(w) = 1

  Returns
    real scalar in [0, ln(n)]
──────────────────────────────────────────────────────────────────────────────*/

real scalar _compute_entropy(real colvector w)
{
    real scalar entropy
    real colvector w_safe

    w_safe = w :+ (w :== 0)
    entropy = -sum(w :* ln(w_safe))

    return(entropy)
}

/*──────────────────────────────────────────────────────────────────────────────
  _compute_ess()

  Effective sample size (Kish, 1965) for a normalised weight vector:

      ESS = 1 / sum(w_i^2)     when sum(w) = 1

  Bounds:  ESS = n when w_i = 1/n  (uniform);
           ESS = 1 when a single w_k = 1  (degenerate).

  Arguments
    w   N x 1 column vector with sum(w) = 1

  Returns
    real scalar in [1, n]
──────────────────────────────────────────────────────────────────────────────*/

real scalar _compute_ess(real colvector w)
{
    real scalar ess, sum_w_squared

    sum_w_squared = sum(w :^ 2)

    if (sum_w_squared < 1e-16) {
        ess = 1
    }
    else {
        ess = 1 / sum_w_squared
    }

    return(ess)
}

/*──────────────────────────────────────────────────────────────────────────────
  _compute_concentration()

  Concentration index: the smallest fraction k/n of units whose
  cumulative weight (sorted descending) reaches at least 50%.

  Arguments
    w               N x 1 column vector with sum(w) = 1
    concentration   (output) scalar k/n
    top_k           (output) scalar k
──────────────────────────────────────────────────────────────────────────────*/

void _compute_concentration(real colvector w,
                            real scalar concentration,
                            real scalar top_k)
{
    real colvector w_sorted, cumsum_w
    real scalar n

    n = rows(w)

    w_sorted = sort(w, -1)
    cumsum_w = runningsum(w_sorted)

    top_k = sum(cumsum_w :< 0.5) + 1

    if (top_k < 1) top_k = 1
    if (top_k > n) top_k = n

    concentration = top_k / n
}

/*──────────────────────────────────────────────────────────────────────────────
  struct weight_stats

  Container for weight-vector diagnostics.
──────────────────────────────────────────────────────────────────────────────*/

struct weight_stats {
    real scalar n
    real scalar min_val
    real scalar max_val
    real scalar mean_val
    real scalar entropy
    real scalar ess
    real scalar concentration
    real scalar top_k
}

/*──────────────────────────────────────────────────────────────────────────────
  compute_weight_stats()

  Computes descriptive statistics, Shannon entropy, effective sample size,
  and concentration index for a weight vector.

  Arguments
    w      N x 1 weight vector (normalised)
    type   string label ("time" or "unit")

  Returns
    struct weight_stats
──────────────────────────────────────────────────────────────────────────────*/

struct weight_stats scalar compute_weight_stats(real colvector w,
                                                string scalar type)
{
    struct weight_stats scalar stats
    real scalar concentration_val, top_k_val

    stats.n = rows(w)
    stats.min_val = min(w)
    stats.max_val = max(w)
    stats.mean_val = mean(w)

    stats.entropy = _compute_entropy(w)
    stats.ess = _compute_ess(w)

    _compute_concentration(w, concentration_val, top_k_val)
    stats.concentration = concentration_val
    stats.top_k = top_k_val

    return(stats)
}

/*──────────────────────────────────────────────────────────────────────────────
  struct bootstrap_stats

  Container for bootstrap distribution summary statistics.
──────────────────────────────────────────────────────────────────────────────*/

struct bootstrap_stats {
    real scalar mean_val
    real scalar sd
    real scalar skewness
    real scalar kurtosis
    real scalar min_val
    real scalar max_val
    real scalar n_converged
    string scalar skew_label
    string scalar kurt_label
}

/*──────────────────────────────────────────────────────────────────────────────
  compute_bootstrap_stats()

  Distributional diagnostics of bootstrap ATT estimates: mean, standard
  deviation, skewness, and kurtosis.

      skewness = E[(X - mu)^3] / sigma^3
      kurtosis = E[(X - mu)^4] / sigma^4

  Missing values (from non-converged iterations) are excluded.

  Arguments
    att_boot    B x 1 vector of bootstrap ATT estimates
    converged   optional B x 1 convergence indicator

  Returns
    struct bootstrap_stats
──────────────────────────────────────────────────────────────────────────────*/

struct bootstrap_stats scalar compute_bootstrap_stats(
    real colvector att_boot,
    | real colvector converged
)
{
    struct bootstrap_stats scalar stats
    real colvector z, att_valid
    real scalar has_converged, n_total

    has_converged = (args() >= 2 && rows(converged) > 0)

    n_total = rows(att_boot)
    att_valid = select(att_boot, att_boot :< .)

    stats.mean_val = mean(att_valid)
    stats.sd = sqrt(variance(att_valid))
    stats.min_val = min(att_valid)
    stats.max_val = max(att_valid)

    if (has_converged) {
        stats.n_converged = sum(converged)
    }
    else {
        stats.n_converged = rows(att_valid)
    }

    if (stats.sd > 0) {
        z = (att_valid :- stats.mean_val) / stats.sd
        stats.skewness = mean(z :^ 3)
        stats.kurtosis = mean(z :^ 4)
    }
    else {
        stats.skewness = 0
        stats.kurtosis = .
    }

    if (abs(stats.skewness) < 0.5) {
        stats.skew_label = "(approximately symmetric)"
    }
    else if (stats.skewness > 0) {
        stats.skew_label = "(right-skewed)"
    }
    else {
        stats.skew_label = "(left-skewed)"
    }

    if (stats.kurtosis >= .) {
        stats.kurt_label = "(degenerate: all estimates identical)"
    }
    else if (stats.kurtosis >= 2.5 && stats.kurtosis <= 3.5) {
        stats.kurt_label = "(approximately normal)"
    }
    else if (stats.kurtosis > 3.5) {
        stats.kurt_label = "(heavy-tailed)"
    }
    else {
        stats.kurt_label = "(light-tailed)"
    }

    return(stats)
}

/*──────────────────────────────────────────────────────────────────────────────
  _trop_estat_factors_svd()

  Singular value decomposition of a factor matrix stored in e().

  Displays: singular values, variance explained, matrix norms, and
  element-level statistics.

      L = U * diag(sigma) * V'
      variance share:   sigma_i^2 / sum(sigma^2)
      Frobenius norm:   ||L||_F = sqrt(sum(sigma^2))
      nuclear norm:     ||L||_* = sum(sigma_i)

  The nuclear norm is the penalty term in the TROP objective; its
  magnitude relative to the Frobenius norm indicates how much
  regularisation is active.

  If the matrix is numerically zero, abbreviated output is printed.
  When e(effective_rank) is unavailable, the continuous effective rank
  sum(sigma) / sigma_1 is computed as a fallback.

  Arguments
    matname   name of a Stata matrix stored in e()
──────────────────────────────────────────────────────────────────────────────*/

void _trop_estat_factors_svd(string scalar matname)
{
    real matrix L, L_work, Vt
    real vector s
    real scalar i, total_var, tol, n_sv, effective_rank
    real scalar frob_norm, nuclear_norm, max_abs, min_abs
    real scalar T, N, transposed

    L = st_matrix(matname)
    T = rows(L)
    N = cols(L)

    frob_norm = sqrt(sum(L:^2))
    if (frob_norm < 1e-30) {
        printf("{txt}Singular value decomposition:\n")
        printf("{txt}  Effective rank  = {res}%8.2f\n", 0)
        printf("{txt}  Top singular values:\n")
        printf("{txt}  (factor matrix is effectively zero, no variance to decompose)\n")
        printf("\n{txt}Matrix norms:\n")
        printf("{txt}  ||L||_F (Frobenius) = {res}%10.3f\n", frob_norm)
        printf("{txt}  ||L||_* (Nuclear)   = {res}%10.3f\n", 0)
        printf("\n{txt}Element statistics:\n")
        printf("{txt}  max|L_it|           = {res}%10.3f\n", max(abs(L)))
        printf("{txt}  min|L_it|           = {res}%10.3f\n", min(abs(L)))
        return
    }

    /* _svd() requires rows >= cols; transpose if needed */
    transposed = 0
    if (T < N) {
        L_work = L'
        transposed = 1
    }
    else {
        L_work = L
    }

    s = J(0, 1, .)
    Vt = J(0, 0, .)
    _svd(L_work, s, Vt)
    if (length(s) == 0 || hasmissing(s)) {
        printf("{err}SVD decomposition failed. Factor matrix may be degenerate.{txt}\n")
        return
    }

    tol = 1e-10
    total_var = sum(s:^2)
    n_sv = length(s)

    /* Retrieve stored effective rank; fall back to sum(sigma)/sigma_1 */
    effective_rank = .
    if (length(st_numscalar("e(effective_rank)")) == 1) {
        effective_rank = st_numscalar("e(effective_rank)")
    }
    if (effective_rank >= . | effective_rank == .) {
        if (length(s) > 0 && s[1] > 0) {
            effective_rank = sum(s) / s[1]
        }
        else {
            effective_rank = 0
        }
    }

    printf("{txt}Singular value decomposition:\n")
    printf("{txt}  Effective rank  = {res}%8.2f\n", effective_rank)
    printf("{txt}  Top singular values:\n")

    for (i = 1; i <= min((5, n_sv)); i++) {
        if (s[i] < tol) {
            printf("{txt}    σ%g = {res}%9.3f{txt} (< tol, effectively zero)\n",
                   i, s[i])
        }
        else {
            printf("{txt}    σ%g = {res}%9.3f{txt} (explains {res}%5.1f%%{txt} variance)\n",
                   i, s[i], 100 * s[i]^2 / total_var)
        }
    }

    nuclear_norm = sum(s)

    printf("\n{txt}Matrix norms:\n")
    printf("{txt}  ||L||_F (Frobenius) = {res}%10.3f\n", frob_norm)
    printf("{txt}  ||L||_* (Nuclear)   = {res}%10.3f\n", nuclear_norm)

    max_abs = max(abs(L))
    min_abs = min(abs(L))

    printf("\n{txt}Element statistics:\n")
    printf("{txt}  max|L_it|           = {res}%10.3f\n", max_abs)
    printf("{txt}  min|L_it|           = {res}%10.3f\n", min_abs)
}

end
