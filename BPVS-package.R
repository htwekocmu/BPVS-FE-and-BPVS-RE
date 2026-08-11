#' BPVS: Bayesian Panel Variable Selection
#'
#' Bayesian Panel Variable Selection with regularized horseshoe priors.
#' Fixed-effects (unit dummies) and random-effects (random intercepts)
#' specifications with posterior inclusion probabilities, WAIC, Bayes R2,
#' and Bayesian / frequentist Hausman comparisons.
#'
#' The two workhorse fitting functions are \code{\link{fit_bpvs_fe}} (fixed
#' effects, unit dummies) and \code{\link{fit_bpvs_re}} (random intercepts).
#' \code{\link{bpvs_pipeline}} fits both, computes PIPs, ranks the predictors,
#' and compares the two specifications in a single call.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom magrittr %>%
#' @importFrom stats as.formula terms sd median quantile na.omit gaussian rnorm
#' @importFrom utils data read.csv write.csv
#' @importFrom brms brm bf set_prior save_pars bayes_R2
#' @importFrom posterior as_draws_df
#' @importFrom loo waic
#' @importFrom dplyr arrange mutate select filter group_by summarise case_when n desc
#' @importFrom rlang sym .data :=
#' @import ggplot2
NULL

utils::globalVariables(c("Variable", "PIP", "Estimate", "CI_lower", "CI_upper",
                         "Rank", "Selection", "Delta", "Delta_median",
                         "Delta_q025", "Delta_q975", "Pr_delta_not0",
                         "RE_Consistent", "FE_Estimate", "FE_PIP",
                         "FE_Selection", "RE_Estimate", "RE_PIP", "RE_Selection",
                         "PIP_FE", "PIP_RE", "True", "Estimated", "Recovered",
                         "xmin", "xmax"))
