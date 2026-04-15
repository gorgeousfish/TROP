{smcl}
{vieweralsosee "trop estat" "help trop_estat"}{...}
{vieweralsosee "trop predict" "help trop_predict"}{...}
{viewerjumpto "Syntax" "trop##syntax"}{...}
{viewerjumpto "Description" "trop##description"}{...}
{viewerjumpto "Options" "trop##options"}{...}
{viewerjumpto "Remarks" "trop##remarks"}{...}
{viewerjumpto "Examples" "trop##examples"}{...}
{viewerjumpto "Stored results" "trop##results"}{...}
{viewerjumpto "Methods and formulas" "trop##methods"}{...}
{viewerjumpto "References" "trop##references"}{...}
{viewerjumpto "Author" "trop##author"}{...}
{viewerjumpto "Installation" "trop##installation"}{...}
{viewerjumpto "Also see" "trop##alsosee"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{cmd:trop} {hline 2}}Triply Robust Panel Estimator{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:trop} {depvar} {it:treatvar} {ifin}{cmd:,} {opth panelvar(varname)} {opth timevar(varname)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth panelvar(varname)}}panel unit identifier{p_end}
{synopt:{opth timevar(varname)}}time period identifier{p_end}

{syntab:Estimation Method}
{synopt:{opt method(twostep|joint)}}estimation method; default is {cmd:twostep}{p_end}

{syntab:Lambda Grid Settings}
{synopt:{opt grid_style(default|extended)}}preset lambda grid; default is {cmd:default}{p_end}
{synopt:{opth lambda_time_grid(numlist)}}custom lambda_time grid values{p_end}
{synopt:{opth lambda_unit_grid(numlist)}}custom lambda_unit grid values{p_end}
{synopt:{opth lambda_nn_grid(numlist)}}custom lambda_nn grid values{p_end}

{syntab:LOOCV Control}
{synopt:{opth fixedlambda(numlist)}}skip LOOCV; use fixed lambda values (3 values){p_end}
{synopt:{opt max_loocv_samples(#)}}maximum LOOCV subsamples; default is {cmd:0}{p_end}
{synopt:{opt tol(#)}}convergence tolerance; default is {cmd:1e-6}{p_end}
{synopt:{opt maxiter(#)}}maximum iterations; default is {cmd:100}{p_end}

{syntab:Bootstrap Inference}
{synopt:{opt bootstrap(#)}}bootstrap replications; default is {cmd:0} (disabled){p_end}

{syntab:Other}
{synopt:{opt seed(#)}}random number seed; default is {cmd:42}{p_end}
{synopt:{opt level(#)}}confidence level; default is {cmd:c(level)}{p_end}
{synopt:{opt verbose}}display detailed progress information{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:trop} implements the Triply Robust Panel (TROP) estimator proposed by
Athey, Imbens, Qu, and Viviano (2025) for estimating the average treatment
effect on the treated (ATT) in panel data settings.  The estimator combines
unit weights, time weights, and a nuclear-norm regularized low-rank factor
model to predict counterfactual outcomes for treated unit-time pairs.

{pstd}
The triple robustness property (Theorem 5.1) ensures that the bias of the
estimator is bounded by the product of three components: unit imbalance,
time imbalance, and misspecification of the regression adjustment.  The
estimator is consistent if any one of the three components is zero
(Corollary 1).

{pstd}
Two estimation methods are available:

{phang}
{opt twostep} (default) implements Algorithm 2 of the paper, which estimates
individual treatment effects for each treated unit-time pair, allowing for
heterogeneous effects across units and time periods.

{phang}
{opt joint} implements the approach described in Remark 6.1, which assumes
homogeneous treatment effects and estimates a single scalar tau via weighted
least squares.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth panelvar(varname)} specifies the variable that identifies the panel
units (e.g., individuals, firms, countries).

{phang}
{opth timevar(varname)} specifies the variable that identifies the time
periods.

{dlgtab:Estimation Method}

{phang}
{opt method(twostep|joint)} specifies the estimation method.

{pmore}
{opt twostep} (default) uses per-observation weights following Algorithm 2.
This allows for heterogeneous treatment effects and is the recommended
approach.

{pmore}
{opt joint} uses global weights following Remark 6.1. This assumes
homogeneous treatment effects but provides higher precision when the
assumption holds.

{dlgtab:Lambda Grid Settings}

{phang}
{opt grid_style(default|extended)} specifies a preset lambda grid for
LOOCV hyperparameter selection.

{pmore}
{opt default} uses 180 combinations (6x6x5) and is recommended for
most applications.

{pmore}
{opt extended} uses 4,032 combinations (14x16x18) covering all optimal values
from Table 2 of the paper. This provides a finer search but is computationally
more expensive.

{phang}
{opth lambda_time_grid(numlist)} specifies custom values for lambda_time.
Overrides {opt grid_style()} for this parameter.

{phang}
{opth lambda_unit_grid(numlist)} specifies custom values for lambda_unit.
Overrides {opt grid_style()} for this parameter.

{phang}
{opth lambda_nn_grid(numlist)} specifies custom values for lambda_nn (nuclear
norm penalty). Overrides {opt grid_style()} for this parameter.

{pmore}
{bf:Infinity values:} Use {cmd:.} (missing) or {cmd:1e10} to represent
infinity in custom grids. For lambda_time and lambda_unit, infinity means
uniform weights. For lambda_nn, infinity disables the low-rank factor model.

{dlgtab:LOOCV Control}

{phang}
{opth fixedlambda(numlist)} specifies fixed values for (lambda_time,
lambda_unit, lambda_nn), bypassing LOOCV hyperparameter selection entirely.
Exactly 3 non-negative values must be provided.

{pmore}
Example: {cmd:fixedlambda(0.5 1.0 0.1)} sets lambda_time=0.5,
lambda_unit=1.0, lambda_nn=0.1. When this option is specified,
{opt grid_style()} and {opt max_loocv_samples()} are ignored.
{opt tol()} and {opt maxiter()} still apply to the alternating minimization
in the estimation step.

{phang}
{opt max_loocv_samples(#)} specifies the maximum number of control
observations to use in LOOCV. Default is {cmd:0}, which uses all control
observations (no subsampling), corresponding to the full criterion in
Eq. (5). Set to a positive integer (e.g., 500) to subsample for faster
computation on very large panels.

{phang}
{opt tol(#)} specifies the convergence tolerance for the alternating
minimization algorithm. Default is {cmd:1e-6}.

{phang}
{opt maxiter(#)} specifies the maximum number of iterations for the
alternating minimization algorithm. Default is {cmd:100}.

{dlgtab:Bootstrap Inference}

{phang}
{opt bootstrap(#)} specifies the number of bootstrap replications for
variance estimation following Algorithm 3 of the paper. Default is {cmd:0}
(disabled). Set to a positive integer (e.g., 200) to enable bootstrap
inference. The bootstrap constructs each replicate by sampling N_0 control
units with replacement and N_1 treated units with replacement separately.

{dlgtab:Other}

{phang}
{opt seed(#)} specifies the random number seed for reproducibility.
Default is {cmd:42}.

{phang}
{opt level(#)} specifies the confidence level for confidence intervals.
Default is {cmd:c(level)}, typically 95.

{phang}
{opt verbose} displays detailed progress information including LOOCV
diagnostics, convergence status, and timing information.


{marker remarks}{...}
{title:Remarks}

{pstd}
{bf:Weight structures}

{pstd}
The two estimation methods use fundamentally different weight structures,
following Algorithm 2 and Remark 6.1 of the paper respectively:

{phang}
{bf:method(twostep)} — Per-observation weights (Eq. 3, Algorithm 2):{break}
Each treated observation (i,t) receives its own set of weights. Time weights
measure distance from the specific target period:
theta_s = exp(-lambda_time * |t - s|). Unit weights measure the RMSE
between each control unit j and the specific target unit i over common
control periods, excluding the target period t:
omega_j = exp(-lambda_unit * dist_{-t}(j, i)). This permits heterogeneous
treatment effects across units and periods.

{phang}
{bf:method(joint)} — Global weights (Remark 6.1):{break}
A single set of weights is shared across all treated observations. Time
weights measure distance to the center of the treated block:
delta_time[t] = exp(-lambda_time * |t - center|). Unit weights measure
the RMSE from each unit's pre-treatment trajectory to the average treated
trajectory: delta_unit[j] = exp(-lambda_unit * RMSE(j, avg_treated)).
This assumes homogeneous treatment effects but is computationally more
efficient.

{pstd}
{bf:When to use each method:}

{phang2}Use {opt twostep} (default) when treatment effects may vary across
units and time periods, or when the panel has complex treatment assignment
patterns.{p_end}

{phang2}Use {opt joint} when treatment effects are expected to be
homogeneous, when computational speed is important, or when the number of
treated observations is large.{p_end}

{pstd}
{bf:Interpretation of e(mu)}

{pstd}
The content of {cmd:e(mu)} differs between estimation methods:

{phang}
{bf:method(twostep):} {cmd:e(mu)} is set to {cmd:.} (missing). The twostep
model is Y(0) = alpha_i + beta_t + L_{it} without an explicit global
intercept. Unit fixed effects are stored in {cmd:e(alpha)}.

{phang}
{bf:method(joint):} {cmd:e(mu)} returns the global intercept mu (a scalar).
The joint model is Y(0) = mu + alpha_i + beta_t + L_{it} with
identification constraints alpha_1 = beta_1 = 0.

{pstd}
Both parameterizations yield mathematically equivalent counterfactual
predictions.

{pstd}
{bf:Lambda grid comparison}

{col 5}Grid Style{col 20}lambda_time{col 35}lambda_unit{col 50}lambda_nn{col 62}Total
{col 5}{hline 65}
{col 5}default{col 20}6 values{col 35}6 values{col 50}5 values{col 62}180
{col 5}extended{col 20}14 values{col 35}16 values{col 50}18 values{col 62}4,032
{col 5}{hline 65}

{pstd}
The {opt extended} grid covers all optimal values from Table 2 of the paper:

{col 5}Dataset{col 25}lambda_unit{col 38}lambda_time{col 51}lambda_nn
{col 5}{hline 60}
{col 5}CPS log-wage{col 25}0{col 38}0.1{col 51}0.9
{col 5}CPS urate{col 25}1.6{col 38}0.35{col 51}0.011
{col 5}PWT{col 25}0.3{col 38}0.4{col 51}0.006
{col 5}Germany{col 25}1.2{col 38}0.2{col 51}0.011
{col 5}Basque{col 25}0{col 38}0.35{col 51}0.006
{col 5}Smoking{col 25}0.25{col 38}0.4{col 51}0.011
{col 5}Boatlift{col 25}0.2{col 38}0.2{col 51}0.151
{col 5}{hline 60}


{marker examples}{...}
{title:Examples}

{pstd}
{bf:Setup: Generate test panel data}

{phang2}{cmd:. clear all}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. set obs 1000}{p_end}
{phang2}{cmd:. gen id = ceil(_n/10)}{p_end}
{phang2}{cmd:. bysort id: gen t = _n}{p_end}
{phang2}{cmd:. gen y = rnormal() + 0.5*(id>80)*(t>7)}{p_end}
{phang2}{cmd:. gen d = (id > 80 & t > 7)}{p_end}

{pstd}
{bf:Basic usage (twostep method, default)}

{phang2}{cmd:. trop y d, panelvar(id) timevar(t) seed(42)}{p_end}

{pstd}
{bf:Joint method}

{phang2}{cmd:. trop y d, panelvar(id) timevar(t) method(joint) seed(42)}{p_end}

{pstd}
{bf:With bootstrap inference}

{phang2}{cmd:. trop y d, panelvar(id) timevar(t) bootstrap(200) seed(42)}{p_end}

{pstd}
{bf:Using extended grid for finer search}

{phang2}{cmd:. trop y d, panelvar(id) timevar(t) grid_style(extended) seed(42)}{p_end}

{pstd}
{bf:Custom lambda grid}

{phang2}{cmd:. trop y d, panelvar(id) timevar(t) lambda_nn_grid(0 0.1 1e10) seed(42)}{p_end}

{pstd}
{bf:Fixed lambda (skip LOOCV)}

{phang2}{cmd:. trop y d, panelvar(id) timevar(t) fixedlambda(0.5 1.0 0.1) seed(42)}{p_end}

{pstd}
{bf:Verbose output}

{phang2}{cmd:. trop y d, panelvar(id) timevar(t) seed(42) verbose}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:trop} stores the following in {cmd:e()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:e(att)}}ATT point estimate{p_end}
{synopt:{cmd:e(se)}}bootstrap standard error{p_end}
{synopt:{cmd:e(t)}}t statistic (att/se){p_end}
{synopt:{cmd:e(ci_lower)}}confidence interval lower bound{p_end}
{synopt:{cmd:e(ci_upper)}}confidence interval upper bound{p_end}
{synopt:{cmd:e(pvalue)}}two-sided p-value{p_end}
{synopt:{cmd:e(df_r)}}degrees of freedom for t-distribution{p_end}
{synopt:{cmd:e(ci_lower_percentile)}}bootstrap percentile CI lower bound{p_end}
{synopt:{cmd:e(ci_upper_percentile)}}bootstrap percentile CI upper bound{p_end}
{synopt:{cmd:e(mu)}}global intercept ({cmd:.} for twostep, scalar for joint){p_end}
{synopt:{cmd:e(lambda_time)}}selected or fixed lambda_time{p_end}
{synopt:{cmd:e(lambda_unit)}}selected or fixed lambda_unit{p_end}
{synopt:{cmd:e(lambda_nn)}}selected or fixed lambda_nn{p_end}
{synopt:{cmd:e(loocv_score)}}LOOCV objective Q(lambda_hat){p_end}
{synopt:{cmd:e(loocv_n_valid)}}number of successful LOOCV fits{p_end}
{synopt:{cmd:e(loocv_n_attempted)}}total LOOCV fit attempts{p_end}
{synopt:{cmd:e(loocv_n_control_total)}}total control observations for LOOCV{p_end}
{synopt:{cmd:e(loocv_fail_rate)}}LOOCV failure rate (0 to 1){p_end}
{synopt:{cmd:e(loocv_subsampled)}}1 if LOOCV used subsampling, 0 otherwise{p_end}
{synopt:{cmd:e(loocv_used)}}1 if LOOCV was performed, 0 if skipped{p_end}
{synopt:{cmd:e(loocv_first_failed_t)}}first failed observation time (-1 if none){p_end}
{synopt:{cmd:e(loocv_first_failed_i)}}first failed observation unit (-1 if none){p_end}
{synopt:{cmd:e(max_loocv_samples)}}maximum LOOCV subsamples requested{p_end}
{synopt:{cmd:e(n_lambda_time)}}number of lambda_time grid values{p_end}
{synopt:{cmd:e(n_lambda_unit)}}number of lambda_unit grid values{p_end}
{synopt:{cmd:e(n_lambda_nn)}}number of lambda_nn grid values{p_end}
{synopt:{cmd:e(n_grid_combinations)}}total grid combinations{p_end}
{synopt:{cmd:e(n_grid_per_cycle)}}grid points per coordinate descent cycle{p_end}
{synopt:{cmd:e(effective_rank)}}effective rank of L (sum(s)/s[1]){p_end}
{synopt:{cmd:e(n_iterations)}}number of iterations{p_end}
{synopt:{cmd:e(converged)}}convergence indicator (1=yes, 0=no){p_end}
{synopt:{cmd:e(n_obs_estimated)}}successfully estimated observations (twostep){p_end}
{synopt:{cmd:e(n_obs_failed)}}failed observations (twostep){p_end}
{synopt:{cmd:e(N_units)}}number of panel units N{p_end}
{synopt:{cmd:e(N_periods)}}number of time periods T{p_end}
{synopt:{cmd:e(N_treated_units)}}number of ever-treated units{p_end}
{synopt:{cmd:e(N_obs)}}total number of observations{p_end}
{synopt:{cmd:e(N_treat)}}number of treated unit-period pairs (W=1){p_end}
{synopt:{cmd:e(N_control)}}number of control observations (W=0){p_end}
{synopt:{cmd:e(N_control_units)}}number of never-treated units{p_end}
{synopt:{cmd:e(T_treat_periods)}}number of periods with treatment{p_end}
{synopt:{cmd:e(bootstrap_reps)}}bootstrap replications requested{p_end}
{synopt:{cmd:e(n_bootstrap_valid)}}successful bootstrap iterations{p_end}
{synopt:{cmd:e(alpha_level)}}significance level for CI{p_end}
{synopt:{cmd:e(level)}}confidence level (e.g., 95){p_end}
{synopt:{cmd:e(seed)}}random number seed used{p_end}
{synopt:{cmd:e(balanced)}}1 if panel is balanced, 0 otherwise{p_end}
{synopt:{cmd:e(miss_rate)}}fraction of missing observations{p_end}
{synopt:{cmd:e(min_pre_treated)}}minimum pre-treatment periods for treated units{p_end}
{synopt:{cmd:e(min_valid_pairs)}}minimum common control periods across pairs{p_end}
{synopt:{cmd:e(has_switching)}}1 if treatment switching detected{p_end}
{synopt:{cmd:e(max_switches)}}maximum treatment switches observed{p_end}
{synopt:{cmd:e(data_validated)}}1 if data validation passed{p_end}
{synopt:{cmd:e(time_min)}}minimum time period value{p_end}
{synopt:{cmd:e(time_max)}}maximum time period value{p_end}
{synopt:{cmd:e(time_range)}}range of time periods{p_end}
{synopt:{cmd:e(n_pre_periods)}}number of pre-treatment periods{p_end}
{synopt:{cmd:e(n_post_periods)}}number of post-treatment periods{p_end}

{p2col 5 28 32 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector (1x1, ATT){p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix (1x1; bootstrap only){p_end}
{synopt:{cmd:e(alpha)}}unit fixed effects (N x 1){p_end}
{synopt:{cmd:e(beta)}}time fixed effects (T x 1){p_end}
{synopt:{cmd:e(factor_matrix)}}low-rank factor matrix L (T x N){p_end}
{synopt:{cmd:e(tau)}}individual treatment effects (twostep only){p_end}
{synopt:{cmd:e(theta)}}time weights (T x 1; twostep only){p_end}
{synopt:{cmd:e(omega)}}unit weights (N x 1; twostep only){p_end}
{synopt:{cmd:e(delta_time)}}global time weights (T x 1; joint only){p_end}
{synopt:{cmd:e(delta_unit)}}global unit weights (N x 1; joint only){p_end}
{synopt:{cmd:e(bootstrap_estimates)}}bootstrap ATT estimates (B x 1){p_end}
{synopt:{cmd:e(lambda_time_grid)}}lambda_time grid values searched{p_end}
{synopt:{cmd:e(lambda_unit_grid)}}lambda_unit grid values searched{p_end}
{synopt:{cmd:e(lambda_nn_grid)}}lambda_nn grid values searched{p_end}
{synopt:{cmd:e(lambda_grid)}}Cartesian product of lambda grids (K x 3){p_end}
{synopt:{cmd:e(cv_curve)}}LOOCV scores at grid points (K x 4){p_end}

{p2col 5 28 32 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}"trop"{p_end}
{synopt:{cmd:e(cmdline)}}full command line as typed{p_end}
{synopt:{cmd:e(method)}}"twostep" or "joint"{p_end}
{synopt:{cmd:e(depvar)}}dependent variable name{p_end}
{synopt:{cmd:e(treatvar)}}treatment variable name{p_end}
{synopt:{cmd:e(panelvar)}}panel variable name{p_end}
{synopt:{cmd:e(timevar)}}time variable name{p_end}
{synopt:{cmd:e(grid_style)}}grid style used{p_end}
{synopt:{cmd:e(treatment_pattern)}}treatment pattern detected{p_end}
{synopt:{cmd:e(vcetype)}}"Bootstrap" (with bootstrap){p_end}
{synopt:{cmd:e(estat_cmd)}}"trop_estat"{p_end}
{synopt:{cmd:e(title)}}"TROP Estimator"{p_end}
{synopt:{cmd:e(predict)}}"trop_p"{p_end}

{p2col 5 28 32 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


{marker methods}{...}
{title:Methods and formulas}

{pstd}
The TROP estimator predicts counterfactual outcomes for treated unit-time
pairs using a working model Y_{it}(0) = alpha_i + beta_t + L_{it}, where
alpha_i are unit fixed effects, beta_t are time fixed effects, and L is a
low-rank factor component.

{pstd}
For each treated observation (i,t), the parameters are estimated by solving
a weighted nuclear-norm penalized regression (Eq. 2):

{p 8 8 2}
(alpha_hat, beta_hat, L_hat) = argmin_{alpha,beta,L} sum_{j,s}
theta_s^{i,t} * omega_j^{i,t} * (1-W_{js}) *
(Y_{js} - alpha_j - beta_s - L_{js})^2 + lambda_nn * ||L||_*

{pstd}
where ||L||_* denotes the nuclear norm of L.  The weights are defined as
exponential distance-decay functions (Eq. 3):

{p 8 8 2}
theta_s^{i,t} = exp(-lambda_time * |t - s|)

{p 8 8 2}
omega_j^{i,t} = exp(-lambda_unit * dist_{-t}(j, i))

{pstd}
The unit distance is the RMSE between the outcome trajectories of units j
and i over common control periods, excluding the target period t:

{p 8 8 2}
dist_{-t}(j, i) = sqrt( sum_u 1{u!=t}(1-W_{iu})(1-W_{ju})(Y_{iu}-Y_{ju})^2
/ sum_u 1{u!=t}(1-W_{iu})(1-W_{ju}) )

{pstd}
The treatment effect for each treated observation (i,t) is estimated as:

{p 8 8 2}
tau_hat_{it} = Y_{it} - alpha_hat_i - beta_hat_t - L_hat_{it}

{pstd}
The ATT is the average of individual treatment effects over all treated
unit-time pairs (Algorithm 2):

{p 8 8 2}
tau_hat = (1 / sum_{i,t} W_{it}) * sum_{i,t} W_{it} * tau_hat_{it}

{pstd}
The triple robustness property (Theorem 5.1) states that the bias satisfies:

{p 8 8 2}
|E[tau_hat - tau | L]| <= ||Delta^u||_2 * ||Delta^t||_2 * ||B||_*

{pstd}
where Delta^u is unit imbalance (the discrepancy between the weighted
average of control unit loadings and the treated unit loading), Delta^t is
time imbalance (the analogous discrepancy for time factor loadings), and B
captures regression adjustment misspecification.  The estimator is
consistent if any one of the three components is zero (Corollary 1):
(a) balance over unit loadings, (b) balance over factor loadings, or
(c) correct regression adjustment specification.

{pstd}
Tuning parameters (lambda_time, lambda_unit, lambda_nn) are selected via
leave-one-out cross-validation (LOOCV) minimizing (Eq. 5):

{p 8 8 2}
Q(lambda) = sum_{i,t} (1 - W_{it}) * (tau_hat_{it}(lambda))^2

{pstd}
The grid search uses coordinate descent (footnote 2 of the paper): each
parameter is optimized in turn while holding the other two at their current
optimal values.

{pstd}
Bootstrap variance estimation follows Algorithm 3: units are resampled
with replacement within the treated and control groups separately, and
the TROP estimator is recomputed on each bootstrap sample.  The bootstrap
variance is:

{p 8 8 2}
V_hat = (1/B) * sum_{b=1}^{B} (tau_hat^{(b)} - tau_bar)^2

{pstd}
where tau_hat^{(b)} is the ATT estimate from the b-th bootstrap sample
and tau_bar is the mean of the bootstrap estimates.


{marker references}{...}
{title:References}

{phang}
Athey, S., G. W. Imbens, Z. Qu, and D. Viviano. 2025.
Triply robust panel estimators.
{it:arXiv preprint arXiv:2508.21536}.
{p_end}


{marker author}{...}
{title:Author}

{pstd}
Xuanyu Cai{break}
City University of Macau{break}
xuanyuCAI@outlook.com

{pstd}
Wenli Xu{break}
City University of Macau{break}
wlxu@cityu.edu.mo


{marker installation}{...}
{title:Installation}

{pstd}
{bf:From GitHub:}

{phang2}{cmd:. net install trop, from("https://raw.githubusercontent.com/gorgeousfish/TROP/main") replace}{p_end}

{pstd}
{bf:Local installation:}

{phang2}{cmd:. net install trop, from("/path/to/trop_stata") replace}{p_end}

{pstd}
{bf:Verify installation:}

{phang2}{cmd:. trop, version}{p_end}
{phang2}{cmd:. trop_check}{p_end}

{pstd}
The package includes a precompiled plugin for macOS ARM64 (Apple Silicon).
Other platforms (macOS Intel, Linux x64, Windows x64) can be built from the
Rust source code.  See the project repository for build instructions.


{marker alsosee}{...}
{title:Also see}

{psee}
Online: {helpb trop}, {helpb trop_estat}, {helpb trop_predict}
{p_end}
