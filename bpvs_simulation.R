# ======================================================================
# bpvs_simulation.R -- synthetic panel data generation and recovery
# Pathairat Pastpipatkul and Htwe Ko (2025)
# ======================================================================

#' Simulate panel data with known coefficients
#'
#' Generates a balanced panel with unit-specific random intercepts and
#' correlated Gaussian predictors, so that the true coefficients (including
#' zero entries) are known. Useful for validating parameter recovery and
#' selection performance.
#'
#' The design matrix is drawn from
#' \deqn{\mathbf{x}_{it} \sim N(\mathbf{0}, \mathbf{\Sigma}),
#' \quad \mathbf{\Sigma}_{ij} = \rho^{|i-j|},}{x_it ~ N(0, Sigma),
#' Sigma_ij = rho^|i-j|,}
#' i.e. a Toeplitz autocorrelation structure, and the outcome is
#' \deqn{y_{it} = \alpha_i + \mathbf{x}_{it}^\top \boldsymbol{\beta}
#' + \varepsilon_{it},}{y_it = alpha_i + x_it' beta + eps_it,}
#' with \eqn{\alpha_i \sim N(0, \tau_{\text{unit}}^2)}{alpha_i ~ N(0,
#' tau_unit^2)} and \eqn{\varepsilon_{it} \sim N(0, \sigma^2)}{eps_it ~ N(0,
#' sigma^2)}.
#'
#' @param N_units Number of cross-sectional units.
#' @param T_periods Number of time periods per unit.
#' @param beta_true Numeric vector of true slope coefficients.
#' @param sigma Idiosyncratic error standard deviation.
#' @param tau_unit Standard deviation of the unit random intercepts.
#' @param rho Autocorrelation of the Toeplitz predictor covariance.
#' @param predictors Optional character vector of predictor names; defaults to
#'   \code{paste0("x", 1:p)}.
#'
#' @return A data frame with columns \code{id}, \code{year}, \code{Income}
#'   (the outcome), and one column per predictor.
#' @export
simulate_panel_data <- function(N_units = 50, T_periods = 10,
                                beta_true = c(1.5, -0.8, 0, 0.5, 0, 0),
                                sigma = 0.5, tau_unit = 0.4, rho = 0.3,
                                predictors = NULL) {
  p <- length(beta_true)
  N <- N_units * T_periods

  if (is.null(predictors)) predictors <- paste0("x", 1:p)

  X <- MASS::mvrnorm(N, mu = rep(0, p),
                     Sigma = outer(1:p, 1:p, function(i, j) rho^abs(i - j)))
  colnames(X) <- predictors

  unit_id <- rep(1:N_units, each = T_periods)
  alpha <- rnorm(N_units, 0, tau_unit)

  y <- numeric(N)
  for (i in 1:N) {
    y[i] <- alpha[unit_id[i]] + sum(X[i, ] * beta_true) + rnorm(1, 0, sigma)
  }

  data.frame(
    id     = factor(unit_id),
    year   = rep(1:T_periods, times = N_units),
    Income = y,
    X
  )
}

#' Build a recovery summary from a fitted BPVS model
#'
#' Compares the posterior medians and credible intervals of the slopes with the
#' true coefficients used to generate the data, and reports posterior inclusion
#' probabilities.
#'
#' @param fit A fitted BPVS model (FE or RE).
#' @param effects_label Label for the model, e.g. \code{"FE"} or \code{"RE"}.
#' @param true_beta Named numeric vector of true coefficients.
#' @param data The (simulated) data frame used to fit the model.
#' @param outcome_var Name of the outcome column.
#' @param effect_threshold Multiplier of \code{sd(y)} defining "effectively
#'   zero" (see \code{\link{compute_pip}}).
#'
#' @return A data frame with columns \code{Model}, \code{Variable}, \code{True},
#'   \code{Estimated}, \code{CI_lower}, \code{CI_upper}, \code{CI_95},
#'   \code{PIP}, \code{Selected_PIP}, \code{True_NonZero}, and
#'   \code{Recovery_Error}.
#' @export
recovery_summary <- function(fit, effects_label = "FE", true_beta,
                             data, outcome_var = "Income",
                             effect_threshold = 0.01) {
  post <- as_draws_df(fit)
  pattern <- if (effects_label == "FE") "^b_eta_" else "^b_"
  beta_cols <- setdiff(grep(pattern, names(post), value = TRUE), "b_Intercept")
  var_names <- gsub("^b_eta_|^b_", "", beta_cols)

  beta_med  <- sapply(beta_cols, function(col) median(post[[col]]))
  beta_q025 <- sapply(beta_cols, function(col) quantile(post[[col]], 0.025))
  beta_q975 <- sapply(beta_cols, function(col) quantile(post[[col]], 0.975))

  eps <- effect_threshold * sd(data[[outcome_var]])
  pip <- sapply(beta_cols, function(col) mean(abs(post[[col]]) > eps))

  data.frame(
    Model          = effects_label,
    Variable       = var_names,
    True           = true_beta[var_names],
    Estimated      = round(beta_med, 4),
    CI_lower       = round(beta_q025, 4),
    CI_upper       = round(beta_q975, 4),
    CI_95          = sprintf("[%.4f, %.4f]", beta_q025, beta_q975),
    PIP            = round(pip, 4),
    Selected_PIP   = ifelse(pip > 0.5, "Yes", "No"),
    True_NonZero   = ifelse(true_beta[var_names] != 0, "Yes", "No"),
    Recovery_Error = round(beta_med - true_beta[var_names], 4),
    stringsAsFactors = FALSE
  )
}

#' Recovery diagnostics
#'
#' Computes selection accuracy, RMSE, and misclassification rate of a recovery
#' summary.
#'
#' @param summ A recovery summary from \code{\link{recovery_summary}}.
#'
#' @return A data frame with columns \code{Model}, \code{Correct_Selection},
#'   \code{RMSE}, and \code{Misclassification}.
#' @export
recovery_metrics <- function(summ) {
  correct_sel <- summ$Selected_PIP == ifelse(summ$True_NonZero == "Yes", "Yes", "No")
  rmse <- sqrt(mean((summ$Estimated - summ$True)^2))
  misclass <- mean((summ$True != 0) != (summ$PIP > 0.5))
  data.frame(
    Model = unique(summ$Model),
    Correct_Selection = sprintf("%d/%d (%.0f%%)", sum(correct_sel), nrow(summ),
                                100 * mean(correct_sel)),
    RMSE = round(rmse, 4),
    Misclassification = round(misclass, 4)
  )
}
