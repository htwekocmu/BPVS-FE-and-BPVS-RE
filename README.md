# BPVS

**Bayesian Panel Variable Selection with Regularized Horseshoe Priors**

`BPVS` implements Bayesian Panel Variable Selection in fixed-effects (unit
dummies, `BPVS-FE`) and random-effects (random intercepts, `BPVS-RE`)
specifications. A regularized horseshoe (continuous spike-and-slab) prior is
placed on the slope coefficients; the package computes posterior inclusion
probabilities (PIPs), variable rankings, WAIC, Bayesian R-squared, and Bayesian
and frequentist Hausman comparisons between the two panel specifications.
Estimation is implemented with [`brms`](https://paul-buerkner.github.io/brms/)
on [`Stan`](https://mc-stan.org/).

The methods are described in

> Pastpipatkul, P. and Ko, H. (2025). *Bayesian Panel Variable Selection with
> the Regularized Horseshoe Prior: Fixed-Effects and Random-Effects
> Specifications, Estimation, and Empirical Validation.*

## Installation

You need R (>= 4.0). The package is installed from GitHub:

```r
install.packages("remotes")
remotes::install_github("htwekocmu/BPVS-FE-and-BPVS-RE", ref = "add-r-package-and-shiny")
library(BPVS)
```

You also need a Stan backend for `brms` (either `rstan` or, preferably,
`cmdstanr`):

```r
install.packages("rstan")          # option 1 (simplest)
# or
install.packages("cmdstanr")       # option 2 (recommended)
cmdstanr::install_cmdstan()
```

## Quick start

### 1. Fit the model on your own panel

```r
library(BPVS)

# Read your panel data: one row per (id, time), with numeric predictors
dat <- read.csv("your_panel.csv")

# Outcome, id, and time columns:
res <- bpvs_pipeline(Growth ~ IE + IH + Debt + FDI + INF + POP + TOP + UR + EXR + GR,
                     data    = dat,
                     effects = c("fe", "re"),   # fit both specifications
                     outcome_var = "Growth",
                     id_var  = "id",
                     time_var = "year")

# Inspect results
res$rank_fe          # BPVS-FE ranking (PIP, selection)
res$rank_re          # BPVS-RE ranking
res$comparison$bayesian_hausman   # posterior of beta_FE - beta_RE
res$comparison$hausman_plm        # frequentist Hausman test
```

### 2. Try the bundled example data

```r
data(developing)
res <- bpvs_pipeline(Growth ~ IE + IH + Debt + FDI + INF + POP + TOP + UR + EXR + GR,
                     data    = developing,
                     effects = c("fe", "re"),
                     outcome_var = "Growth", id_var = "id", time_var = "year")
```

### 3. Launch the Shiny app

```r
BPVS::run_bpvs_app()
```

The app lets you upload a CSV, choose the outcome/id/time columns and the
specifications, and explore rankings, PIPs, coefficient plots, and Hausman
diagnostics interactively.

### 4. Simulation / validation

```r
set.seed(456)
sim <- simulate_panel_data(N_units = 50, T_periods = 10,
                           beta_true = c(1.5, -0.8, 0, 0.5, 0, 0),
                           sigma = 0.5, tau_unit = 0.4)
res  <- bpvs_pipeline(Income ~ FDI + INF + POP + TOP + UR + EXR, sim,
                      effects = c("fe", "re"), outcome_var = "Income")
```

## Key functions

| Function               | Purpose                                              |
|------------------------|------------------------------------------------------|
| `fit_bpvs_fe()`        | BPVS with fixed effects (unit dummies)               |
| `fit_bpvs_re()`        | BPVS with random intercepts                          |
| `bpvs_pipeline()`      | Fit + rank + compare in one call                     |
| `compute_pip()`        | Posterior inclusion probabilities                    |
| `rank_variables()`     | Rank predictors, label Selected/Weak/Excluded        |
| `compare_bpvs()`       | FE vs RE comparison + Hausman diagnostics            |
| `bayesian_hausman()`   | Posterior of beta_FE - beta_RE per variable          |
| `hausman_plm()`        | Frequentist Hausman test via `plm`                   |
| `simulate_panel_data()`| Generate synthetic panel with known coefficients     |
| `plot_inclusion()`     | PIP lollipop chart                                   |
| `plot_coefficients()`  | Coefficient forest plot with 95% CIs                 |
| `plot_pip_compare()`   | FE vs RE PIP scatter                                 |
| `run_bpvs_app()`       | Launch the Shiny app                                 |

## Model summary

Likelihood:

```
y_it ~ N(mu_it, sigma^2)
mu_it = alpha_i + sum_j beta_j x_itj
```

- **BPVS-RE**: `alpha_i ~ N(0, tau_unit^2)` (random intercepts).
- **BPVS-FE**: `alpha_i` free, estimated as unshrunk unit dummies; the
  horseshoe prior is applied only to the slope block.

Prior on slopes (regularized horseshoe):

```
beta_j ~ N(0, tau^2 c^2 lambda_j^2 / (c^2 + tau^2 lambda_j^2))
lambda_j ~ Cauchy+(0, 1)     # local shrinkage
c^2 ~ Inv-Gamma(nu/2, nu s^2/2)
tau_0 = p0/(p - p0) * sigma / sqrt(n)   # auto-scaled global shrinkage
```

Prior on variance components: `sigma ~ student_t(3, 0, 10)`,
`tau_unit ~ exponential(1)`.

Inclusion: `PIP_j = Pr(|beta_j| > 0.01 * sd(y) | data)`; a variable is
*Selected* if `PIP_j > 0.5`.

## License

GPL-3.
