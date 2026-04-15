# trop

**Triply Robust Panel Estimator for Stata**

[![Stata 17+](https://img.shields.io/badge/Stata-17%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 1.0.0](https://img.shields.io/badge/Version-1.0.0-green.svg)]()

![trop](image/image.png)

## Overview

`trop` implements the **Triply RObust Panel (TROP) estimator** proposed by Athey, Imbens, Qu, and Viviano (2025) for Stata. The estimator combines three components — **unit weights**, **time weights**, and a **nuclear-norm-regularized low-rank regression adjustment** — to estimate average treatment effects on the treated (ATT) in panel data with potentially complex assignment patterns.

In semi-synthetic simulations calibrated to seven real datasets (Table 1 of the paper), TROP achieves the **lowest RMSE in 20 out of 21 specifications**, outperforming DID, SC, SDID, MC, and DIFP estimators across diverse data-generating processes:

| Data set | Outcome | Treatment | N | T | TROP | SDID | SC | DID | MC | DIFP |
|----------|---------|-----------|---|---|------|------|----|-----|-----|------|
| CPS | log-wage | min wage | 50 | 40 | **1.00** | 1.14 | 1.44 | 1.91 | 1.26 | 1.22 |
| CPS | urate | min wage | 50 | 40 | **1.00** | 1.05 | 1.11 | 1.89 | 1.10 | 1.09 |
| PWT | log-GDP | democracy | 111 | 48 | **1.00** | 1.44 | 1.59 | 7.85 | 1.76 | 1.54 |
| Germany | GDP | random | 17 | 44 | **1.00** | 1.46 | 2.82 | 3.58 | 1.56 | 2.46 |
| Basque | GDP | random | 18 | 43 | **1.00** | 1.02 | 4.55 | 9.11 | 1.70 | 2.47 |
| Smoking | packs pc | random | 39 | 31 | **1.00** | 1.22 | 1.48 | 2.16 | 1.14 | 1.45 |
| Boatlift | log-wage | random | 44 | 19 | **1.00** | 1.34 | 1.62 | 1.35 | 1.04 | 1.62 |

<sub>Normalized RMSE from Table 1 of Athey et al. (2025). Full results cover 21 specifications across 7 datasets.</sub>

**Features:**

- **Triple Robust Estimation** — Asymptotically unbiased if *any one of* unit weights, time weights, or the regression adjustment removes biases (Theorem 5.1)
- **Two Estimation Methods** — Twostep (per-observation, heterogeneous effects; Algorithm 2) and Joint (weighted least squares, homogeneous effect)
- **Leave-One-Out Cross-Validation** — Data-driven selection of tuning parameters via coordinate-cycling LOOCV (Algorithm 1)
- **Bootstrap Inference** — Stratified unit block bootstrap for variance estimation and confidence intervals (Algorithm 3)
- **General Assignment Patterns** — Handles staggered adoption, switching treatments, and arbitrary binary treatment matrices
- **Post-Estimation Diagnostics** — 7 `estat` subcommands and 8 `predict` types for comprehensive analysis
- **High-Performance Backend** — Core computation in Rust via compiled plugin; no external dependencies

## Key Concepts

### Triple Robustness

The TROP estimator combines three components, each targeting a different source of confounding:

| Component | Role | Controlled by |
|-----------|------|--------------|
| **Unit weights** $\omega_j$ | Upweight control units similar to treated units | $\lambda_{\text{unit}}$ |
| **Time weights** $\theta_s$ | Upweight time periods close to the treatment period | $\lambda_{\text{time}}$ |
| **Low-rank factor model** $\mathbf{L}$ | Capture unobserved interactive fixed effects | $\lambda_{nn}$ |

The key insight is that the estimator's bias is bounded by the **product** of three imbalance terms (Theorem 5.1). If any single component successfully removes bias, the overall bias vanishes — hence *triple* robustness. This multiplicative bound is strictly tighter than the additive bounds of DID, SC, or SDID.

### Special Cases

The TROP framework nests existing estimators:

| Setting | Recovers |
|---------|----------|
| $\lambda_{nn} = \infty$, uniform weights | **DID / TWFE** |
| Uniform weights, $\lambda_{nn} < \infty$ | **Matrix Completion** |
| $\lambda_{nn} = \infty$, specific unit/time weights | **SC / SDID** |

## Requirements

- Stata 17.0 or later
- No additional Stata packages required
- Precompiled plugin included for macOS ARM64 (Apple Silicon); other platforms can be built from the Rust source

## Installation

### Install the package

```stata
net install trop, from("https://raw.githubusercontent.com/gorgeousfish/TROP/main") replace
```

This automatically installs:
- All commands and help files
- Pre-compiled Mata library
- macOS ARM64 plugin

### Verify Installation

```stata
trop, version
trop_check
```

## Quick Start with Examples

All examples below use the **CPS log-wage dataset** — 50 US states × 40 years (1979–2018) of state-level log wages, where `d` flags state-years in which a minimum wage increase was in effect. This is one of the seven benchmark datasets from Athey et al. (2025).

Download the dataset first:

```stata
net get trop                 /* downloads datasets to the current directory */
use cps_logwage.dta, clear
```

Or load directly from GitHub:

```stata
use "https://raw.githubusercontent.com/gorgeousfish/TROP/main/data/cps_logwage.dta", clear
```

**Dataset:** N = 50, T = 40, 2,000 observations. `y` = log state-level wage; `d` = minimum wage treatment (8 treated state-year cells, 0.4%); `id` = state identifier; `t` = year.

### Example 1: Fixed Hyperparameters (Twostep)

Using the paper's recommended values for CPS log-wage (Table S.1 of Athey et al. 2025):

```stata
use "https://raw.githubusercontent.com/gorgeousfish/TROP/main/data/cps_logwage.dta", clear
trop y d, panelvar(id) timevar(t) fixedlambda(0.1 0 0.9)
```

Output:

```
------------------------------------------------------------------------------
TROP Estimation Results
------------------------------------------------------------------------------
Method:            twostep
Grid style:        default (17 grid points/cycle, coordinate descent)

Panel dimensions:  N = 50, T = 40
Observations:      2000
Treated:           8 ( 0.4%)

Fixed hyperparameters (LOOCV skipped):
  lambda_time =   0.1000
  lambda_unit =   0.0000
  lambda_nn   =   0.9000

Treatment Effect (ATT):
  tau     =     0.031406
```

### Example 2: LOOCV-Selected Hyperparameters

```stata
use "https://raw.githubusercontent.com/gorgeousfish/TROP/main/data/cps_logwage.dta", clear
trop y d, panelvar(id) timevar(t)
```

Output:

```
------------------------------------------------------------------------------
TROP Estimation Results
------------------------------------------------------------------------------
Method:            twostep
Grid style:        default (17 grid points/cycle, coordinate descent)

Panel dimensions:  N = 50, T = 40
Observations:      2000
Treated:           8 ( 0.4%)

Selected hyperparameters (via LOOCV):
  lambda_time =   0.5000
  lambda_unit =   5.0000
  lambda_nn   =   0.1000
  Q(lambda_hat) =   3.456566

Treatment Effect (ATT):
  tau     =     0.033572
```

**Note:** LOOCV on N = 50, T = 40 typically takes 20–40 minutes. For faster evaluation use `fixedlambda()` or limit search with `max_loocv_samples(500)`.

### Example 3: Bootstrap Inference

```stata
use "https://raw.githubusercontent.com/gorgeousfish/TROP/main/data/cps_logwage.dta", clear
trop y d, panelvar(id) timevar(t) ///
    fixedlambda(0.1 0 0.9) bootstrap(200) seed(42)
```

Output:

```
Treatment Effect (ATT):
  tau     =     0.031406
  SE     =     0.015716
  t      =       1.9984
  p-value=       0.0858
  95% CI = [   -0.005756,     0.068568]
```

### Example 4: Joint Method

```stata
use "https://raw.githubusercontent.com/gorgeousfish/TROP/main/data/cps_logwage.dta", clear
trop y d, panelvar(id) timevar(t) ///
    method(joint) fixedlambda(0.1 0 0.9)
```

Output:

```
Treatment Effect (ATT):
  tau     =     0.031406

Global intercept:
  mu     =     5.154320
```

### Example 5: Post-Estimation Workflow

```stata
use "https://raw.githubusercontent.com/gorgeousfish/TROP/main/data/cps_logwage.dta", clear
trop y d, panelvar(id) timevar(t) fixedlambda(0.1 0 0.9)

* Predict counterfactual outcomes and treatment effects
predict y0_hat, y0
predict te_hat, te

* Diagnostics
estat summarize
estat weights
estat factors
```

`estat summarize` output:

```
-----------------------------------------------------------------
Estimation sample summary
-----------------------------------------------------------------
  Number of observations:        2000    (balanced panel)
  Number of units (N):             50
  Number of periods (T):           40
  Missing rate:                   0.0%

Treatment structure:
  Treated observations:             8    (  0.4%)
  Control observations:          1992    ( 99.6%)
  Treated units:                    8    ( 16.0%)
  Treated periods:                  1    (  2.5%)
  Pattern:                   multiple_treated_simultaneous
-----------------------------------------------------------------
```

To also inspect LOOCV diagnostics, run without `fixedlambda()`:

```stata
use "https://raw.githubusercontent.com/gorgeousfish/TROP/main/data/cps_logwage.dta", clear
trop y d, panelvar(id) timevar(t)
estat loocv
```

### Example 6: Standalone Bootstrap (Post-Estimation)

```stata
use "https://raw.githubusercontent.com/gorgeousfish/TROP/main/data/cps_logwage.dta", clear

* Estimate without bootstrap first (faster iteration)
trop y d, panelvar(id) timevar(t) fixedlambda(0.1 0 0.9)

* Then add bootstrap inference separately
trop_bootstrap, nreps(200) seed(42)
```

Output:

```
------------------------------------------------------------
TROP Bootstrap Inference Results
------------------------------------------------------------
ATT estimate:       0.031406
Bootstrap SE:       0.015716
95% CI:       [   -0.005756,     0.068568]
p-value:              0.0858

Bootstrap reps:        200
Valid reps:            200
------------------------------------------------------------
```

### Example 7: PWT Log-GDP Panel (111 Countries × 48 Years)

For a large panel, use the Penn World Tables democracy dataset:

```stata
use "https://raw.githubusercontent.com/gorgeousfish/TROP/main/data/pwt_loggdp.dta", clear

* Paper's hyperparameters for PWT; add maxiter(200) for large panels
trop y d, panelvar(id) timevar(t) fixedlambda(0.4 0.3 0.006) maxiter(200)
```

Output:

```
------------------------------------------------------------------------------
TROP Estimation Results
------------------------------------------------------------------------------
Method:            twostep
Grid style:        default (17 grid points/cycle, coordinate descent)

Panel dimensions:  N = 111, T = 48
Observations:      5328
Treated:           29 ( 0.5%)

Fixed hyperparameters (LOOCV skipped):
  lambda_time =   0.4000
  lambda_unit =   0.3000
  lambda_nn   =   0.0060

Treatment Effect (ATT):
  tau     =    -0.013552
```

**Available datasets** (download via `net get trop`):

| File | Description | N | T |
|------|-------------|---|---|
| `cps_logwage.dta` | CPS state-level log wage (min wage treatment) | 50 | 40 |
| `cps_urate.dta` | CPS state-level unemployment rate (min wage treatment) | 50 | 40 |
| `pwt_loggdp.dta` | Penn World Tables log GDP (democracy transition) | 111 | 48 |
| `germany_gdp.dta` | Abadie & Gardeazabal (2003) West Germany GDP | 17 | 44 |
| `basque_gdp.dta` | Abadie (2003) Basque Country GDP | 18 | 43 |
| `smoking_packs.dta` | California Prop 99 cigarette consumption | 39 | 31 |

> `germany_gdp`, `basque_gdp`, and `smoking_packs` have `d = 0` throughout — they are raw outcome panels designed for semi-synthetic simulation (as in Table 1 of Athey et al. 2025). Assign a treatment indicator before running `trop`.

## Recommended Workflow

```
Step 1: Validate     →  trop_validate depvar treatvar, panelvar() timevar()
Step 2: Estimate     →  trop depvar treatvar, panelvar() timevar()
Step 3: Diagnose     →  estat summarize / estat weights / estat loocv / estat factors
Step 4: Inference    →  trop_bootstrap, nreps(1000)    — or —    bootstrap() in Step 2
```

**Step 1 — Validate data structure.** `trop_validate` checks panel balance, treatment patterns, missing data, and minimum sample requirements before estimation. This step is also performed automatically inside `trop`.

```stata
trop_validate y d, panelvar(id) timevar(t)
```

**Step 2 — Estimate ATT.** Choose `method(twostep)` (default) for heterogeneous effects, or `method(joint)` for faster homogeneous estimation. LOOCV selects hyperparameters automatically.

| Method | Use when | Speed |
|--------|----------|-------|
| `twostep` | Heterogeneous effects; general assignment | Slower |
| `joint` | Homogeneous effect; simultaneous adoption only | Faster |

**Step 3 — Diagnose results.**

| Subcommand | What to check |
|------------|---------------|
| `estat summarize` | Panel dimensions, treatment pattern, balance |
| `estat weights` | Whether weights concentrate on few units/periods |
| `estat loocv` | Convergence, failure rate, selected lambdas |
| `estat factors` | Effective rank, explained variance by top singular values |

**Step 4 — Conduct inference.** Use `bootstrap()` inline or `trop_bootstrap` post-estimation. The bootstrap resamples units within treatment strata (Algorithm 3).

```stata
* Inline (re-runs estimation + bootstrap together)
trop y d, panelvar(id) timevar(t) bootstrap(1000) seed(42)

* Post-estimation (uses stored estimation results)
trop_bootstrap, nreps(1000) seed(42)
```

## Tutorial

An interactive Jupyter notebook tutorial is included as ancillary material:

```stata
net get trop
```

The downloaded file `10_trop_stata.ipynb` covers data generation, estimation with both methods, all `estat` diagnostics, prediction, and a real-data CPS example.

## Commands

| Command          | Description                                       |
| ---------------- | ------------------------------------------------- |
| `trop`           | Main estimation command (twostep or joint method)  |
| `trop_bootstrap` | Standalone bootstrap inference (post-estimation)   |
| `trop_validate`  | Panel data structure validation                    |
| `trop_check`     | Environment and installation check                 |

Post-estimation (available after `trop`):

| Command          | Description                                       |
| ---------------- | ------------------------------------------------- |
| `estat`          | Diagnostics dispatcher (7 subcommands; see below) |
| `predict`        | Prediction dispatcher (8 types; see below)        |

## Options

### trop Options

**Required:**

| Option              | Description                      |
| ------------------- | -------------------------------- |
| `panelvar(varname)` | Unit (panel) identifier variable |
| `timevar(varname)`  | Time identifier variable         |

**Optional:**

| Option                       | Description                                                     | Default    |
| ---------------------------- | --------------------------------------------------------------- | ---------- |
| `method(string)`             | Estimation method: `twostep` or `joint`                         | `twostep`  |
| `grid_style(string)`         | Lambda grid style: `default` or `extended`                      | `default`  |
| `lambda_time_grid(numlist)`  | User-specified grid for lambda_time                             | auto       |
| `lambda_unit_grid(numlist)`  | User-specified grid for lambda_unit                             | auto       |
| `lambda_nn_grid(numlist)`    | User-specified grid for lambda_nn                               | auto       |
| `fixedlambda(numlist)`       | Fix (lambda_time lambda_unit lambda_nn); skip LOOCV             | —          |
| `max_loocv_samples(integer)` | LOOCV subsample size (0 = full sample)                          | `0`        |
| `tol(real)`                  | Convergence tolerance for iterative estimation                  | `1e-6`     |
| `maxiter(integer)`           | Maximum number of iterations                                    | `100`      |
| `bootstrap(integer)`         | Number of bootstrap replications (0 = none)                     | `0`        |
| `seed(integer)`              | Random number generator seed                                    | `42`       |
| `level(cilevel)`             | Confidence level for intervals                                  | `c(level)` |
| `verbose`                    | Display detailed diagnostic output                              | off        |

**Grid styles:**
- `default` — 6 × 6 × 5 = 180 grid combinations, 17 evaluations per coordinate-descent cycle
- `extended` — 14 × 16 × 18 = 4,032 combinations, 48 evaluations per cycle (finer search, slower)

### trop_bootstrap Options

| Option              | Description                                                     | Default    |
| ------------------- | --------------------------------------------------------------- | ---------- |
| `nreps(integer)`    | Number of bootstrap replications                                | `1000`     |
| `level(real)`       | Confidence level in percent (10–99.99)                          | `c(level)` |
| `seed(integer)`     | Random number generator seed                                    | `42`       |
| `maxiter(integer)`  | Maximum iterations per replication                              | `100`      |
| `tol(real)`         | Convergence tolerance per replication                           | `1e-6`     |
| `verbose`           | Display progress information                                    | off        |

## Stored Results

### Scalars

*Core estimates:*

| Scalar          | Description                                    |
| --------------- | ---------------------------------------------- |
| `e(att)`        | ATT point estimate                             |
| `e(se)`         | Bootstrap standard error                       |
| `e(t)`          | t statistic (att/se)                           |
| `e(pvalue)`     | Two-sided p-value                              |
| `e(ci_lower)`   | Lower bound of confidence interval             |
| `e(ci_upper)`   | Upper bound of confidence interval             |
| `e(mu)`         | Global intercept (joint only; missing for twostep) |

*Tuning parameters:*

| Scalar           | Description                                   |
| ---------------- | --------------------------------------------- |
| `e(lambda_time)` | Selected lambda_time                          |
| `e(lambda_unit)` | Selected lambda_unit                          |
| `e(lambda_nn)`   | Selected lambda_nn                            |
| `e(loocv_score)` | Optimal LOOCV score Q(lambda_hat)             |

*Sample information:*

| Scalar                | Description                              |
| --------------------- | ---------------------------------------- |
| `e(N_units)`          | Number of units                          |
| `e(N_periods)`        | Number of time periods                   |
| `e(N_obs)`            | Total observations                       |
| `e(N_treat)`          | Number of treated observations           |
| `e(N_control)`        | Number of control observations           |
| `e(N_treated_units)`  | Number of treated units                  |
| `e(N_control_units)`  | Number of control units                  |
| `e(T_treat_periods)`  | Number of treatment periods              |
| `e(bootstrap_reps)`   | Number of bootstrap replications         |

*Convergence:*

| Scalar             | Description                                |
| ------------------ | ------------------------------------------ |
| `e(n_iterations)`  | Number of iterations                       |
| `e(converged)`     | Convergence indicator (1/0)                |
| `e(n_obs_estimated)` | Successfully estimated observations (twostep only) |
| `e(n_obs_failed)`  | Failed observations (twostep, if > 0)      |

*LOOCV diagnostics:*

| Scalar                     | Description                             |
| -------------------------- | --------------------------------------- |
| `e(loocv_n_valid)`         | Number of valid LOOCV evaluations       |
| `e(loocv_n_attempted)`     | Number of attempted LOOCV evaluations   |
| `e(loocv_n_control_total)` | Total control observations for LOOCV    |
| `e(loocv_fail_rate)`       | LOOCV failure rate                      |
| `e(loocv_subsampled)`      | Whether LOOCV subsampling was used (1/0)|
| `e(max_loocv_samples)`     | LOOCV subsample size setting            |
| `e(loocv_used)`            | Whether LOOCV was performed (1/0)       |
| `e(seed)`                  | RNG seed used                           |

*Grid information:*

| Scalar                   | Description                                        |
| ------------------------ | -------------------------------------------------- |
| `e(n_lambda_time)`       | Number of lambda_time grid values                  |
| `e(n_lambda_unit)`       | Number of lambda_unit grid values                  |
| `e(n_lambda_nn)`         | Number of lambda_nn grid values                    |
| `e(n_grid_combinations)` | Total Cartesian grid combinations                  |
| `e(n_grid_per_cycle)`    | Grid evaluations per coordinate-descent cycle      |

*Other:*

| Scalar                    | Description                             |
| ------------------------- | --------------------------------------- |
| `e(balanced)`             | Balanced panel indicator (1/0)          |
| `e(miss_rate)`            | Missing data rate                       |
| `e(alpha_level)`          | Significance level for CI               |
| `e(effective_rank)`       | Effective rank of factor matrix         |
| `e(n_bootstrap_valid)`    | Number of valid bootstrap replications  |
| `e(data_validated)`       | Data validation indicator (1/0)         |

### Macros

| Macro                  | Description                              |
| ---------------------- | ---------------------------------------- |
| `e(cmd)`               | `"trop"`                                 |
| `e(cmdline)`           | Full command line as issued              |
| `e(method)`            | `"twostep"` or `"joint"`                |
| `e(grid_style)`        | `"default"`, `"extended"`, or `"custom"` |
| `e(depvar)`            | Dependent variable name                  |
| `e(treatvar)`          | Treatment variable name                  |
| `e(panelvar)`          | Panel variable name                      |
| `e(timevar)`           | Time variable name                       |
| `e(vcetype)`           | `"Bootstrap"` or `""`                    |
| `e(estat_cmd)`         | `"trop_estat"`                           |
| `e(treatment_pattern)` | Treatment assignment pattern description |

### Matrices

| Matrix                   | Description                                          |
| ------------------------ | ---------------------------------------------------- |
| `e(b)`                   | Coefficient vector (1×1, ATT)                        |
| `e(V)`                   | Variance-covariance matrix (1×1; requires bootstrap) |
| `e(alpha)`               | Unit fixed effects (N×1)                             |
| `e(beta)`                | Time fixed effects (T×1)                             |
| `e(factor_matrix)`       | Low-rank factor matrix L (T×N)                       |
| `e(tau)`                 | Individual treatment effects (twostep only)          |
| `e(bootstrap_estimates)` | Bootstrap distribution (B×1; requires bootstrap)     |
| `e(theta)`               | Time weights (twostep only)                          |
| `e(omega)`               | Unit weights (twostep only)                          |
| `e(delta_time)`          | Time weights (joint only)                            |
| `e(delta_unit)`          | Unit weights (joint only)                            |
| `e(lambda_time_grid)`    | Lambda time grid values                              |
| `e(lambda_unit_grid)`    | Lambda unit grid values                              |
| `e(lambda_nn_grid)`      | Lambda nuclear norm grid values                      |

## Post-Estimation

### estat Subcommands

| Subcommand          | Abbreviation | Description                               |
| ------------------- | ------------ | ----------------------------------------- |
| `estat summarize`   | `sum`        | Sample structure and treatment allocation  |
| `estat vce`         |              | Variance-covariance matrix display         |
| `estat sensitivity` | `sens`       | Hyperparameter sensitivity analysis        |
| `estat weights`     | `weight`     | Unit and time weight diagnostics           |
| `estat bootstrap`   | `boot`       | Bootstrap distribution diagnostics         |
| `estat loocv`       |              | LOOCV hyperparameter selection diagnostics |
| `estat factors`     |              | Factor matrix (L) SVD analysis             |

### predict Types

After `trop`, use `predict newvar, type` to generate predictions:

| Type        | Description                                           |
| ----------- | ----------------------------------------------------- |
| `y0`        | Counterfactual outcome Y(0) **[default]**             |
| `y1`        | Counterfactual outcome Y(1)                           |
| `te`        | Treatment effect (treated obs only)                   |
| `residuals` | Residuals Y - Y(0)                                    |
| `alpha`     | Unit fixed effects                                    |
| `beta`      | Time fixed effects                                    |
| `mu`        | Global intercept (joint only)                         |
| `xb`        | Linear prediction (equivalent to `y0`)                |

## Methodology

### The TROP Estimator

The TROP estimator models the potential control outcome as $Y_{it}(0) = \alpha_i + \beta_t + L_{it} + \epsilon_{it}$, where $\alpha_i$ are unit fixed effects, $\beta_t$ are time fixed effects, $L_{it}$ is a low-rank factor component, and $\epsilon_{it}$ is idiosyncratic noise.

For each treated unit–time pair $(i,t)$, the estimator predicts the counterfactual outcome by solving a weighted nuclear-norm penalized regression (Eq. 2 of the paper):

$$(\hat{\alpha}, \hat{\beta}, \hat{\mathbf{L}}) = \arg\min_{\alpha, \beta, \mathbf{L}} \sum_{j=1}^{N} \sum_{s=1}^{T} \theta_s^{i,t} \omega_j^{i,t} (1-W_{js})(Y_{js} - \alpha_j - \beta_s - L_{js})^2 + \lambda_{nn} \|\mathbf{L}\|_*$$

where the weights exhibit exponential decay (Eq. 3):

$$\theta_s^{i,t} = \exp(-\lambda_{\text{time}} \cdot |t - s|), \qquad \omega_j^{i,t} = \exp(-\lambda_{\text{unit}} \cdot \text{dist}_{-t}^{\text{unit}}(j, i))$$

The unit distance measures RMSE of outcome differences over shared control periods:

$$\text{dist}_{-t}^{\text{unit}}(j,i) = \left(\frac{\sum_{u} \mathbf{1}\{u \neq t\}(1-W_{iu})(1-W_{ju})(Y_{iu}-Y_{ju})^2}{\sum_{u} \mathbf{1}\{u \neq t\}(1-W_{iu})(1-W_{ju})}\right)^{1/2}$$

The treatment effect is then $\hat{\tau}_{it} = Y_{it} - \hat{\alpha}_i - \hat{\beta}_t - \hat{L}_{it}$.

This formulation encompasses DID, SC, MC, and SDID all as special cases. For $\lambda_{nn} = \infty$ and $\omega_j = \theta_s = 1$, we recover the DID/TWFE estimator. For $\omega_j = \theta_s = 1$ and $\lambda_{nn} < \infty$, we recover the MC estimator. For $\lambda_{nn} = \infty$ with specific unit and time weights, we recover SC and SDID.

### Triple Robustness Property

The bias satisfies a multiplicative bound (Theorem 5.1):

$$\left|\mathbb{E}[\hat{\tau} - \tau \mid \mathbf{L}]\right| \leq \|\Delta^{\mathbf{u}}(\omega, \Gamma)\|_2 \times \|\Delta^{\mathbf{t}}(\theta, \Lambda)\|_2 \times \|B\|_*$$

where $\Delta^{\mathbf{u}}$ is unit imbalance, $\Delta^{\mathbf{t}}$ is time imbalance, and $B$ captures regression adjustment misspecification. The estimator is consistent if **any one** of the three terms is negligible (Corollary 1):

1. Balance over unit factor loadings ($\|\Delta^{\mathbf{u}}\|_2 \approx 0$)
2. Balance over time factor loadings ($\|\Delta^{\mathbf{t}}\|_2 \approx 0$)
3. Correct regression adjustment specification ($\|B\|_* \approx 0$)

This *multiplicative* bias bound is strictly tighter than the *additive* bounds that govern DID, SC, and SDID, giving TROP stronger robustness properties.

### Tuning Parameter Selection

The triplet $(\lambda_{\text{time}}, \lambda_{\text{unit}}, \lambda_{nn})$ is selected via leave-one-out cross-validation minimizing (Eq. 5):

$$Q(\lambda) = \sum_{i=1}^{N} \sum_{t=1}^{T} (1 - W_{it})(\hat{\tau}_{it}(\lambda))^2$$

This is equivalent to choosing the tuning parameters with the smallest out-of-sample squared error for predicting the potential outcome under control on the control observations. The grid search uses coordinate descent (Algorithm 1): each parameter is optimized in turn while holding the other two at their most recently selected values, then the cycle repeats until convergence.

### Estimation Methods

- **Twostep** (Algorithm 2): Per-observation estimation allowing heterogeneous treatment effects. For each treated pair $(i,t)$, the model is fitted as if $(i,t)$ were the only treated observation. The ATT is $\hat{\tau} = \frac{1}{\sum_{i,t} W_{it}} \sum_{i,t} W_{it} \hat{\tau}_{it}$.
- **Joint** (Remark 6.1): Weighted least squares with a single scalar treatment effect $\tau$, assuming homogeneous effects across all treated unit–time pairs. Uses global weights shared across all treated observations. Computationally more efficient when the homogeneity assumption holds.

### Bootstrap Inference

Variance estimation follows Algorithm 3: stratified unit block bootstrap that separately resamples treated and control units with replacement. For each replication $b = 1, \ldots, B$, the full estimation procedure (including LOOCV if applicable) is repeated to obtain $\hat{\tau}^{(b)}$, and the bootstrap variance is:

$$\hat{V}_{\tau} = \frac{1}{B} \sum_{b=1}^{B} (\hat{\tau}^{(b)} - \bar{\hat{\tau}})^2$$

where $\bar{\hat{\tau}} = \frac{1}{B}\sum_{b=1}^B \hat{\tau}^{(b)}$ is the mean of bootstrap estimates.

## Architecture

The package uses a four-layer design for performance and numerical accuracy:

```
┌─────────────────────────────────────────────────┐
│  Stata User Interface  (trop.ado)               │
│  - Syntax parsing, option handling              │
├─────────────────────────────────────────────────┤
│  Mata Interface Layer                           │
│  - Input validation and data conversion         │
│  - e() result storage                           │
├─────────────────────────────────────────────────┤
│  C Bridge Plugin                                │
│  - Pointer conversion, error code mapping       │
├─────────────────────────────────────────────────┤
│  Rust Core                                      │
│  - Distance matrices, weight computation        │
│  - LOOCV grid search, SVD estimation            │
│  - Stratified unit block bootstrap              │
└─────────────────────────────────────────────────┘
```

All numerical computation (LOOCV, SVD, bootstrap) is performed in Rust for speed and precision. The Mata layer handles data validation and result storage. Users do not need the Rust toolchain — the pre-compiled plugin is included.

## References

Athey, S., Imbens, G., Qu, Z., & Viviano, D. (2025). Triply robust panel estimators. arXiv preprint arXiv:2508.21536.

## Authors

**Stata Implementation:**

- **Xuanyu Cai**, City University of Macau
  Email: [xuanyuCAI@outlook.com](mailto:xuanyuCAI@outlook.com)
- **Wenli Xu**, City University of Macau
  Email: [wlxu@cityu.edu.mo](mailto:wlxu@cityu.edu.mo)

**Methodology:**

- **Susan Athey**, Stanford University
- **Guido Imbens**, Stanford University
- **Zhaonan Qu**, Columbia University
- **Davide Viviano**, Harvard University

## License

AGPL-3.0. See [LICENSE](LICENSE) for details.

## Citation

If you use this package in your research, please cite both the methodology paper and the Stata implementation:

**APA Format:**

> Cai, X., & Xu, W. (2025). *trop: Stata module for Triply Robust Panel estimation* [Computer software]. GitHub. https://github.com/gorgeousfish/TROP
>
> Athey, S., Imbens, G., Qu, Z., & Viviano, D. (2025). Triply robust panel estimators. arXiv preprint arXiv:2508.21536.

**BibTeX:**

```bibtex
@software{trop2025stata,
  title={trop: Stata module for Triply Robust Panel estimation},
  author={Xuanyu Cai and Wenli Xu},
  year={2025},
  version={1.0.0},
  url={https://github.com/gorgeousfish/TROP}
}

@article{athey2025triply,
  title={Triply robust panel estimators},
  author={Athey, Susan and Imbens, Guido and Qu, Zhaonan and Viviano, Davide},
  journal={arXiv preprint arXiv:2508.21536},
  year={2025}
}
```

## See Also

- Original paper by Athey, Imbens, Qu & Viviano: https://arxiv.org/abs/2508.21536
- Related Stata packages: [`sdid`](https://github.com/Daniel-Pailanir/sdid) (Synthetic DID), [`diddesign`](https://github.com/gorgeousfish/diddesign) (Double DID)
