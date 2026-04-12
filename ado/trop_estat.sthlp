{smcl}
{vieweralsosee "trop" "help trop"}{...}
{vieweralsosee "trop estat" "help trop_estat"}{...}
{vieweralsosee "trop predict" "help trop_predict"}{...}
{viewerjumpto "Syntax" "trop_estat##syntax"}{...}
{viewerjumpto "Description" "trop_estat##description"}{...}
{viewerjumpto "Subcommands" "trop_estat##subcommands"}{...}
{viewerjumpto "Examples" "trop_estat##examples"}{...}
{viewerjumpto "Stored results" "trop_estat##results"}{...}
{viewerjumpto "Methods and formulas" "trop_estat##methods"}{...}
{viewerjumpto "References" "trop_estat##references"}{...}
{viewerjumpto "Author" "trop_estat##author"}{...}
{title:Title}

{p2colset 5 24 26 2}{...}
{p2col:{cmd:estat} (after {cmd:trop}) {hline 2}}Postestimation statistics for trop{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:estat} {it:subcommand} [{cmd:,} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr:subcommand}
{synoptline}
{synopt:{opt sum:marize}}sample and treatment structure summary{p_end}
{synopt:{opt weight:s}}weight diagnostics{p_end}
{synopt:{opt loocv}}LOOCV hyperparameter selection diagnostics{p_end}
{synopt:{opt factors}}factor matrix (L) SVD analysis{p_end}
{synopt:{opt boot:strap}}bootstrap distribution diagnostics{p_end}
{synopt:{opt sens:itivity}}hyperparameter sensitivity analysis{p_end}
{synopt:{opt vce}}variance-covariance matrix{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
Abbreviations are shown by underlining.


{marker description}{...}
{title:Description}

{pstd}
{cmd:estat} displays postestimation diagnostics after {helpb trop}.  Seven
subcommands are provided for inspecting the estimation sample, weight
distributions, hyperparameter selection, factor structure, bootstrap
inference, sensitivity of the LOOCV objective, and the variance-covariance
matrix.

{pstd}
All subcommands require that {cmd:trop} has been executed previously.  If no
{cmd:trop} estimates are found in {cmd:e()}, error code 301 is returned.


{marker subcommands}{...}
{title:Subcommands}

{dlgtab:estat summarize}

{p 8 16 2}
{cmd:estat summarize} [{cmd:,} {opt detailed}]

{pstd}
Displays a summary of the estimation sample structure and treatment
distribution, including the number of units (N), periods (T), total
observations, treated and control observation counts, treatment adoption
pattern, outcome variable descriptive statistics, and estimation method.

{pstd}
Options:

{phang}
{opt detailed} displays a guide for tabulating the unit-by-period treatment
distribution in the original dataset.

{dlgtab:estat weights}

{p 8 16 2}
{cmd:estat weights} [{cmd:,} {opt heatmap} {opt detailed}]

{pstd}
Displays weight diagnostics.  Output differs by estimation method:

{phang}
{bf:Twostep method:} A note is displayed indicating that weights vary by
treated observation.  Summary statistics (mean, standard deviation, minimum,
maximum) are reported for the time weights (theta, T x 1) and unit weights
(omega, N x 1) associated with the first treated observation.

{phang}
{bf:Joint method:} Summary statistics are reported for the global time
weights (delta_time, T x 1) and unit weights (delta_unit, N x 1).

{pstd}
Options:

{phang}
{opt heatmap} generates a heatmap of the weight matrix.  If {cmd:e(W_mat)} is
available it is used directly; otherwise an outer-product approximation is
constructed from the weight vectors.  Requires the user-written command
{cmd:heatplot} (install via {cmd:ssc install heatplot}).

{phang}
{opt detailed} displays the 10th, 25th, 50th, 75th, and 90th percentiles of
the weight vectors.

{dlgtab:estat loocv}

{p 8 16 2}
{cmd:estat loocv}

{pstd}
Displays LOOCV hyperparameter selection diagnostics, including:

{phang2}- Selected hyperparameters (lambda_time, lambda_unit, lambda_nn){p_end}
{phang2}- LOOCV objective value Q(lambda_hat){p_end}
{phang2}- Valid fits count and percentage{p_end}
{phang2}- Grid style used{p_end}
{phang2}- First failed observation coordinates (if any){p_end}
{phang2}- Warning if failure rate exceeds 10%{p_end}

{dlgtab:estat factors}

{p 8 16 2}
{cmd:estat factors}

{pstd}
Displays an analysis of the low-rank factor matrix L via singular value
decomposition (SVD), including:

{phang2}- Matrix dimensions (T x N){p_end}
{phang2}- Effective rank (singular values above tolerance){p_end}
{phang2}- Top singular values with variance explained{p_end}
{phang2}- Matrix norms (Frobenius, Nuclear){p_end}

{dlgtab:estat bootstrap}

{p 8 16 2}
{cmd:estat bootstrap} [{cmd:,} {opt graph}]

{pstd}
Displays bootstrap distribution diagnostics.  This subcommand requires that
the {opt bootstrap()} option was specified in the original {cmd:trop} command.

{pstd}
Output includes:

{phang2}- Number of bootstrap samples and valid samples{p_end}
{phang2}- ATT distribution statistics (mean, standard deviation, median, IQR){p_end}
{phang2}- Percentiles (2.5%, 5%, 25%, 50%, 75%, 95%, 97.5%){p_end}
{phang2}- Normality test (Jarque-Bera statistic and p-value){p_end}

{pstd}
Options:

{phang}
{opt graph} plots a histogram of the bootstrap distribution with a normal
density overlay, the point estimate, and confidence interval reference lines.

{dlgtab:estat sensitivity}

{p 8 16 2}
{cmd:estat sensitivity} [{cmd:,} {opt graph} {opt table(#)}]

{pstd}
Displays hyperparameter sensitivity analysis based on the LOOCV grid search
results.  Reports the search range for each regularization parameter, the
selected optimum, the LOOCV score, and boundary diagnostics.

{pstd}
When full grid evaluation data are available (more than half of the Cartesian
product has valid CV scores), a ranked table of top grid points by CV loss is
displayed together with sensitivity metrics (coefficient of variation, loss
range, relative range).

{pstd}
Options:

{phang}
{opt graph} plots the CV loss curve.  Available only when full grid evaluation
data exist and a single regularization parameter varies.

{phang}
{opt table(#)} specifies the number of top grid points to display; default is
{cmd:table(10)}.

{dlgtab:estat vce}

{p 8 16 2}
{cmd:estat vce} [{cmd:,} {opt correlation}]

{pstd}
Displays the variance-covariance matrix of the estimates.  Because the TROP
estimator produces a scalar ATT, the matrix {cmd:e(V)} is 1 x 1.  The
standard error, variance estimation method, and bootstrap replication count
(if applicable) are also reported.

{pstd}
Options:

{phang}
{opt correlation} displays the correlation matrix instead of the
variance-covariance matrix.  For the 1 x 1 case this is trivially 1.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. clear all}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. set obs 1000}{p_end}
{phang2}{cmd:. gen id = ceil(_n/10)}{p_end}
{phang2}{cmd:. bysort id: gen t = _n}{p_end}
{phang2}{cmd:. gen y = rnormal() + 0.5*(id>80)*(t>7)}{p_end}
{phang2}{cmd:. gen d = (id > 80 & t > 7)}{p_end}

{pstd}Run trop with bootstrap{p_end}
{phang2}{cmd:. trop y d, panelvar(id) timevar(t) seed(42) bootstrap(200)}{p_end}

{pstd}Sample summary{p_end}
{phang2}{cmd:. estat summarize}{p_end}

{pstd}LOOCV diagnostics{p_end}
{phang2}{cmd:. estat loocv}{p_end}

{pstd}Weight diagnostics{p_end}
{phang2}{cmd:. estat weights}{p_end}

{pstd}Weight diagnostics with percentiles{p_end}
{phang2}{cmd:. estat weights, detailed}{p_end}

{pstd}Factor matrix analysis{p_end}
{phang2}{cmd:. estat factors}{p_end}

{pstd}Bootstrap distribution{p_end}
{phang2}{cmd:. estat bootstrap}{p_end}

{pstd}Bootstrap distribution with histogram{p_end}
{phang2}{cmd:. estat bootstrap, graph}{p_end}

{pstd}Hyperparameter sensitivity{p_end}
{phang2}{cmd:. estat sensitivity}{p_end}

{pstd}Variance-covariance matrix{p_end}
{phang2}{cmd:. estat vce}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
The {cmd:estat} subcommands do not store additional results in {cmd:e()}.
All displayed statistics are computed from the {cmd:e()} results stored by
{helpb trop}.  See {helpb trop##results:trop stored results} for the complete
list of scalars, macros, and matrices available after estimation.


{marker methods}{...}
{title:Methods and formulas}

{pstd}
{bf:estat weights}

{pstd}
For the twostep method, time weights theta_s and unit weights omega_j are
computed as exponential-decay functions of distance (Eq. 3 of Athey et al.,
2025):

{p 8 8 2}
theta_s = exp(-lambda_time * |t - s|)

{p 8 8 2}
omega_j = exp(-lambda_unit * dist(j, i))

{pstd}
where dist(j, i) is the root-mean-squared difference in outcomes between
units j and i over shared control periods.  For the joint method, the
analogous global weights are delta_time (T x 1) and delta_unit (N x 1).

{pstd}
{bf:estat loocv}

{pstd}
The LOOCV objective Q(lambda) is defined as (Eq. 5):

{p 8 8 2}
Q(lambda) = sum_{i,t} (1 - W_{it}) * (tau_hat_{it}(lambda))^2

{pstd}
where tau_hat_{it}(lambda) is the estimated pseudo-treatment effect for
control observation (i, t).  Hyperparameters are selected via coordinate
descent (cycling) over the grid of (lambda_time, lambda_unit, lambda_nn).

{pstd}
{bf:estat factors}

{pstd}
The effective rank of the factor matrix L is computed via SVD as
sum(s) / s[1], where s is the vector of singular values.  This continuous
measure indicates the degree to which the low-rank component concentrates on
a small number of factors.

{pstd}
{bf:estat bootstrap}

{pstd}
Bootstrap variance is estimated following Algorithm 3 of Athey et al. (2025).
The cluster bootstrap resamples N_0 control units and N_1 treated units with
replacement, re-estimating the TROP estimator for each replicate.  The
bootstrap variance is:

{p 8 8 2}
V_hat = (1/B) * sum_{b=1}^{B} (tau_hat^(b) - tau_bar)^2

{pstd}
Normality of the bootstrap distribution is assessed via the Jarque-Bera
statistic: JB = (B/6) * [S^2 + (1/4)(K - 3)^2], where S is skewness and K
is kurtosis.  Under the null of normality, JB follows a chi-squared
distribution with 2 degrees of freedom.

{pstd}
{bf:estat sensitivity}

{pstd}
Sensitivity metrics include the coefficient of variation of CV loss across
grid points, the loss range, and the relative range.  Boundary diagnostics
check whether the selected optimum lies at the edge of the search grid.


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


{title:Also see}

{psee}
Online: {helpb trop}, {helpb trop_estat}, {helpb trop_predict}
{p_end}
