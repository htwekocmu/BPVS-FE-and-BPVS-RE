# ======================================================================
# bpvs_model.R -- Bayesian Panel Variable Selection (BPVS) Model
# Pathairat Pastpipatkul and Htwe Ko (2025)
# Revised: FE (unit dummies) and RE (random intercept) variants
#
# PROBABILITY MODEL:
#
#   Likelihood (Panel):
#     y_it ~ N(mu_it, sigma^2)
#     mu_it = alpha_i + sum_{j=1}^p beta_j * x_{itj}
#
#   BPVS-RE (random intercepts):
#     alpha_i ~ N(0, tau_unit^2)             # random effects (country effects)
#     fitted with:  y ~ x1 + ... + xp + (1 | id)
#
#   BPVS-FE (fixed effects, unit dummies):
#     alpha_i treated as free unit-specific intercepts
#     fitted with nonlinear brms formula:
#       y ~ alpha + eta
#       alpha ~ 0 + id                       # unshrunk unit dummies
#       eta  ~ 0 + x1 + ... + xp             # horseshoe only on slopes
#
#   Prior on slopes (Regularized Horseshoe - continuous spike-and-slab):
#     beta_j ~ N(0, tau^2 * lambda_j^2)
#     lambda_j ~ Cauchy+(0, 1)               # local shrinkage
#     tau ~ Cauchy+(0, tau_0)                # global shrinkage
#     tau_0 = p0/(p - p0) * sigma / sqrt(N)  # auto-scaled
#
#   Prior on variance components:
#     sigma ~ student_t(3, 0, 10)
#     tau_unit ~ exponential(1)
#
# IDENTIFICATION:
#   - FE: unit effects identified via dummy coefficients (proper wide prior)
#   - RE: alpha_i identified via proper prior variance tau_unit
#   - Horseshoe identified when p << N or sufficient MCMC iterations
#   - PIP_j = Pr(|beta_j| > eps | data), eps = 0.01 * sd(y)
# ======================================================================

#' Select the Stan backend
#'
#' Returns \code{"cmdstanr"} when the \pkg{cmdstanr} package is installed and
#' its CmdStan installation is present, otherwise \code{"rstan"}.
#'
#' @return A length-one character vector, either \code{"cmdstanr"} or
#'   \code{"rstan"}.
#' @export
bpvs_backend <- function() {
  has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
  if (has_cmdstanr) {
    ok <- tryCatch(
      cmdstanr::cmdstan_version() >= "2.26",
      error = function(e) FALSE
    )
    if (isTRUE(ok)) return("cmdstanr")
  }
  "rstan"
}

#' Extract slope terms from a linear formula
#'
#' Pulls the right-hand-side term labels of a formula, i.e. the names of the
#' predictors that will receive the horseshoe prior.
#'
#' @param formula An object of class \code{formula}.
#' @return A character vector of slope term labels.
#' @keywords internal
#' @noRd
extract_slope_terms <- function(formula) {
  t <- terms(formula)
  attr(t, "term.labels")
}

#' Fit BPVS with fixed effects (unit dummies)
#'
#' Estimates the Bayesian Panel Variable Selection model in which unit-specific
#' intercepts are free parameters implemented as unshrunk dummy coefficients and
#' the regularized horseshoe prior is applied only to the slope block.
#'
#' The model is
#' \deqn{y_{it} \sim N(\alpha_i + \mathbf{x}_{it}^\top \boldsymbol{\beta},
#' \sigma^2),}{y_it ~ N(alpha_i + x_it' beta, sigma^2),}
#' fitted through the nonlinear \pkg{brms} formula
#' \code{y ~ alpha + eta}, \code{alpha ~ 0 + id}, \code{eta ~ 0 + x1 + ... + xp}.
#'
#' @param formula A linear formula of the form \code{y ~ x1 + x2 + ... + xp}.
#' @param data A data frame containing the variables in \code{formula} and the
#'   unit identifier \code{id_var}.
#' @param id_var Name of the column holding the cross-sectional identifier.
#' @param chains,iter,warmup,cores,seed Number of chains, iterations per chain,
#'   warm-up draws, cores, and RNG seed passed to \code{\link[brms]{brm}}.
#' @param adapt_delta,max_treedepth NUTS sampler controls passed to
#'   \code{\link[brms]{brm}}.
#' @param refresh MCMC progress print interval passed to \code{\link[brms]{brm}}.
#' @param backend Stan backend, see \code{\link{bpvs_backend}}. Auto-detected
#'   when \code{NULL}.
#' @param ... Additional arguments passed to \code{\link[brms]{brm}}.
#'
#' @return A fitted \pkg{brms} model object (class \code{brmsfit}) with
#'   attributes \code{bpvs_effects}, \code{bpvs_slopes}, and
#'   \code{bpvs_id_var} attached for downstream use by
#'   \code{\link{compute_pip}}.
#' @export
#' @seealso \code{\link{fit_bpvs_re}} \code{\link{bpvs_pipeline}}
fit_bpvs_fe <- function(formula, data, id_var = "id",
                        chains = 4, iter = 4000, warmup = 2000, cores = 4,
                        seed = 123, adapt_delta = 0.95, max_treedepth = 12,
                        refresh = 500, backend = NULL, ...) {

  if (!id_var %in% names(data)) {
    stop("id_var '", id_var, "' not found in data.")
  }
  data <- as.data.frame(data)
  data[[id_var]] <- factor(data[[id_var]])

  y_name <- all.vars(formula)[1]
  slopes <- extract_slope_terms(formula)
  if (length(slopes) == 0) stop("No predictors on RHS of formula.")

  nl_formula <- bf(
    as.formula(paste(y_name, "~ alpha + eta")),
    as.formula(paste("alpha ~ 0 +", id_var)),
    as.formula(paste("eta ~ 0 +", paste(slopes, collapse = " + "))),
    nl = TRUE
  )

  priors <- c(
    set_prior("horseshoe(df = 1, scale_global = 2, autoscale = TRUE)",
              nlpar = "eta", class = "b"),
    set_prior("normal(0, 10)", nlpar = "alpha", class = "b")
  )

  if (is.null(backend)) backend <- bpvs_backend()

  fit <- brm(
    formula    = nl_formula,
    data       = data,
    family     = gaussian(),
    prior      = priors,
    chains     = chains,
    iter       = iter,
    warmup     = warmup,
    cores      = cores,
    seed       = seed,
    backend    = backend,
    refresh    = refresh,
    save_pars  = save_pars(all = TRUE),
    control    = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    ...
  )
  attr(fit, "bpvs_effects") <- "fe"
  attr(fit, "bpvs_slopes")  <- slopes
  attr(fit, "bpvs_id_var")  <- id_var
  return(fit)
}

#' Fit BPVS with random intercepts (BPVS-RE)
#'
#' Estimates the Bayesian Panel Variable Selection model in which unit-specific
#' intercepts are exchangeable random effects, \eqn{\alpha_i \sim N(0,
#' \tau_{\text{unit}}^2)}{alpha_i ~ N(0, tau_unit^2)}, and the regularized
#' horseshoe prior is applied to the slopes.
#'
#' The model is estimated as \code{y ~ x1 + ... + xp + (1 | id)}.
#'
#' @inheritParams fit_bpvs_fe
#'
#' @return A fitted \pkg{brms} model object (class \code{brmsfit}) with
#'   attributes \code{bpvs_effects}, \code{bpvs_slopes}, and
#'   \code{bpvs_id_var} attached.
#' @export
#' @seealso \code{\link{fit_bpvs_fe}} \code{\link{bpvs_pipeline}}
fit_bpvs_re <- function(formula, data, id_var = "id",
                        chains = 4, iter = 4000, warmup = 2000, cores = 4,
                        seed = 123, adapt_delta = 0.95, max_treedepth = 12,
                        refresh = 500, backend = NULL, ...) {

  if (!id_var %in% names(data)) {
    stop("id_var '", id_var, "' not found in data.")
  }
  data <- as.data.frame(data)
  data[[id_var]] <- factor(data[[id_var]])

  y_name <- all.vars(formula)[1]
  slopes <- extract_slope_terms(formula)
  re_formula <- as.formula(paste(
    y_name, "~", paste(slopes, collapse = " + "), "+ (1 |", id_var, ")"
  ))

  if (is.null(backend)) backend <- bpvs_backend()

  fit <- brm(
    formula    = re_formula,
    data       = data,
    family     = gaussian(),
    prior      = c(set_prior("horseshoe(df = 1, scale_global = 2, autoscale = TRUE)",
                             class = "b")),
    chains     = chains,
    iter       = iter,
    warmup     = warmup,
    cores      = cores,
    seed       = seed,
    backend    = backend,
    refresh    = refresh,
    save_pars  = save_pars(all = TRUE),
    control    = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    ...
  )
  attr(fit, "bpvs_effects") <- "re"
  attr(fit, "bpvs_slopes")  <- slopes
  attr(fit, "bpvs_id_var")  <- id_var
  return(fit)
}

#' Fit a BPVS model (dispatcher)
#'
#' Convenience dispatcher that calls \code{\link{fit_bpvs_fe}} or
#' \code{\link{fit_bpvs_re}} depending on \code{effects}.
#'
#' @inheritParams fit_bpvs_fe
#' @param effects One of \code{"fe"} or \code{"re"}.
#'
#' @return A fitted \pkg{brms} model object.
#' @export
fit_bpvs <- function(formula, data, effects = c("re", "fe"), id_var = "id", ...) {
  effects <- match.arg(effects)
  if (effects == "fe") {
    fit_bpvs_fe(formula, data, id_var = id_var, ...)
  } else {
    fit_bpvs_re(formula, data, id_var = id_var, ...)
  }
}

#' Compute posterior inclusion probabilities
#'
#' Estimates \eqn{\mathrm{PIP}_j = \Pr(|\beta_j| > \varepsilon \mid data)}{PIP_j
#' = Pr(|beta_j| > eps | data)} as the posterior frequency of draws satisfying
#' \eqn{|\beta_j^{(s)}| > \varepsilon}{|beta_j^(s)| > eps}, where
#' \eqn{\varepsilon = \kappa \cdot \operatorname{sd}(y)}{eps = kappa * sd(y)}
#' with \code{effect_threshold = kappa}.
#'
#' Column prefixes differ by specification: fixed effects store slopes as
#' \code{b_eta_<var>}, random effects as \code{b_<var>}.
#'
#' @param fit A fitted BPVS model (output of \code{\link{fit_bpvs_fe}} or
#'   \code{\link{fit_bpvs_re}}).
#' @param data The data frame used to fit the model.
#' @param outcome_var Name of the outcome column.
#' @param effect_threshold Multiplier \eqn{\kappa} of \code{sd(y)} defining the
#'   "effectively zero" radius. Defaults to \code{0.01}.
#'
#' @return A data frame with columns \code{Variable}, \code{Estimate}
#'   (posterior median), \code{CI_lower}, \code{CI_upper}, and \code{PIP}.
#' @export
compute_pip <- function(fit, data, outcome_var = "Income",
                        effect_threshold = 0.01) {
  post <- as_draws_df(fit)
  eps <- effect_threshold * sd(data[[outcome_var]], na.rm = TRUE)

  effects <- attr(fit, "bpvs_effects") %||% "re"
  pattern <- if (effects == "fe") "^b_eta_" else "^b_"

  b_names <- grep(pattern, names(post), value = TRUE)
  b_names <- setdiff(b_names, "b_Intercept")

  pip <- sapply(b_names, function(col) {
    mean(abs(post[[col]]) > eps, na.rm = TRUE)
  })
  med   <- sapply(b_names, function(col) median(post[[col]], na.rm = TRUE))
  q025  <- sapply(b_names, function(col) quantile(post[[col]], 0.025, na.rm = TRUE))
  q975  <- sapply(b_names, function(col) quantile(post[[col]], 0.975, na.rm = TRUE))

  var_names <- gsub("^b_eta_|^b_", "", b_names)

  data.frame(
    Variable   = var_names,
    Estimate   = round(med, 6),
    CI_lower   = round(q025, 6),
    CI_upper   = round(q975, 6),
    PIP        = round(pip, 6),
    stringsAsFactors = FALSE
  )
}

#' Rank predictors by posterior inclusion probability
#'
#' Sorts a PIP table in decreasing order of \code{PIP} and assigns each
#' variable a decision label: \code{Selected} (\code{PIP > 0.5}),
#' \code{Weak} (\code{0.1 < PIP <= 0.5}), or \code{Excluded}
#' (\code{PIP <= 0.1}).
#'
#' @param pip_df A PIP table as returned by \code{\link{compute_pip}}.
#'
#' @return A data frame with columns \code{Variable}, \code{Estimate},
#'   \code{CI_lower}, \code{CI_upper}, \code{PIP}, \code{Rank}, and
#'   \code{Selection}.
#' @export
rank_variables <- function(pip_df) {
  pip_df %>%
    arrange(desc(PIP)) %>%
    mutate(
      Rank      = seq_along(PIP),
      Selection = case_when(
        PIP > 0.5  ~ "Selected",
        PIP > 0.1  ~ "Weak",
        TRUE       ~ "Excluded"
      )
    )
}
