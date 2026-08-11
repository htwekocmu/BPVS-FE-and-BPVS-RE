# ======================================================================
# bpvs_diagnostics.R -- fit statistics, FE vs RE comparison, Hausman
# Pathairat Pastpipatkul and Htwe Ko (2025)
# ======================================================================

#' Null-coalescing operator
#'
#' Returns \code{a} if it is not \code{NULL}, otherwise \code{b}.
#'
#' @param a First object.
#' @param b Fallback object.
#' @return Either \code{a} or \code{b}.
#' @name null-coalesce
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Extract model fit statistics (WAIC and Bayesian R-squared)
#'
#' Computes the Widely Applicable Information Criterion and its standard error
#' via \code{\link[loo]{waic}} and the Bayesian R-squared via
#' \code{\link[brms]{bayes_R2}}.
#'
#' @param fit A fitted BPVS model.
#'
#' @return A named list with elements \code{WAIC}, \code{WAIC_SE}, and
#'   \code{Bayes_R2}.
#' @export
extract_model_stats <- function(fit) {
  w <- tryCatch(waic(fit), error = function(e) NULL)
  r2 <- tryCatch(bayes_R2(fit)[1], error = function(e) NA)
  if (!is.null(w)) {
    waic_val <- w$estimates["waic", "Estimate"]
    waic_se  <- w$estimates["waic", "SE"]
  } else {
    waic_val <- NA; waic_se <- NA
  }
  list(
    WAIC       = round(waic_val, 2),
    WAIC_SE    = round(waic_se, 2),
    Bayes_R2   = round(r2, 4)
  )
}

#' Pretty-print BPVS results
#'
#' Prints a variable ranking table together with model fit statistics and a
#' variable-selection summary to the console.
#'
#' @param ranking A ranking table from \code{\link{rank_variables}}.
#' @param stats A named list from \code{\link{extract_model_stats}}.
#' @param title A title string, e.g. \code{"BPVS-FE"}.
#' @return Invisibly returns \code{NULL}; called for its console output.
#' @export
print_bpvs_results <- function(ranking, stats, title = "BPVS") {
  cat("\n============================================================\n")
  cat(paste(title, "-- BAYESIAN PANEL VARIABLE SELECTION RESULTS\n"))
  cat("============================================================\n")

  cat("\n--- Variable Ranking (by Posterior Inclusion Probability) ---\n")
  print(ranking %>% select(Rank, Variable, Estimate, PIP, Selection),
        row.names = FALSE)

  cat(sprintf("\n--- Model Fit ---\n"))
  cat(sprintf("  WAIC:       %.2f (se = %.2f)\n", stats$WAIC, stats$WAIC_SE))
  cat(sprintf("  Bayes R2:   %.4f\n", stats$Bayes_R2))

  selected <- ranking$Variable[ranking$Selection == "Selected"]
  weak     <- ranking$Variable[ranking$Selection == "Weak"]

  cat("\n--- Variable Selection Summary ---\n")
  if (length(selected) > 0) {
    cat("  Selected (PIP > 0.50):", paste(selected, collapse = ", "), "\n")
  } else {
    cat("  Selected (PIP > 0.50): (none)\n")
  }
  if (length(weak) > 0) {
    cat("  Weak     (PIP > 0.10):", paste(weak, collapse = ", "), "\n")
  }
  excluded <- ranking$Variable[ranking$Selection == "Excluded"]
  if (length(excluded) > 0) {
    cat("  Excluded (PIP <= 0.10):", paste(excluded, collapse = ", "), "\n")
  }
  cat("============================================================\n")
}

#' Bayesian Hausman comparison of FE and RE slope posteriors
#'
#' For each predictor, computes the posterior distribution of
#' \eqn{\Delta_j = \beta_j^{FE} - \beta_j^{RE}}{\Delta_j = beta_j^FE -
#' beta_j^RE} using paired posterior draws from the two fitted models. If the
#' 95\% credible interval of \eqn{\Delta_j} excludes zero, the random-effects
#' estimate of that coefficient is flagged as inconsistent.
#'
#' @param fit_fe A fitted BPVS-FE model.
#' @param fit_re A fitted BPVS-RE model.
#' @param data The data frame used to fit the models.
#' @param outcome_var Name of the outcome column.
#' @param effect_threshold Multiplier of \code{sd(y)} defining "effectively
#'   zero" (see \code{\link{compute_pip}}).
#'
#' @return A data frame with one row per variable and columns
#'   \code{Variable}, \code{Delta_median}, \code{Delta_q025},
#'   \code{Delta_q975}, \code{Pr_delta_not0}, and \code{RE_Consistent}.
#' @export
bayesian_hausman <- function(fit_fe, fit_re, data, outcome_var = "Income",
                             effect_threshold = 0.01) {
  slopes <- attr(fit_fe, "bpvs_slopes") %||% attr(fit_re, "bpvs_slopes")

  post_fe <- as_draws_df(fit_fe)
  post_re <- as_draws_df(fit_re)

  delta_list <- lapply(slopes, function(v) {
    fe_col <- paste0("b_eta_", v)
    re_col <- paste0("b_", v)
    delta <- post_fe[[fe_col]] - post_re[[re_col]]
    data.frame(
      Variable = v,
      Delta    = delta
    )
  })
  delta_df <- do.call(rbind, delta_list)

  delta_df %>%
    group_by(Variable) %>%
    summarise(
      Delta_median  = median(Delta),
      Delta_q025    = quantile(Delta, 0.025),
      Delta_q975    = quantile(Delta, 0.975),
      Pr_delta_not0 = mean(abs(Delta) > effect_threshold * sd(data[[outcome_var]], na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      RE_Consistent = ifelse(Delta_q025 > 0 | Delta_q975 < 0, "No (RE biased)", "Yes")
    )
}

#' Frequentist Hausman test via \pkg{plm}
#'
#' Computes the classical Hausman specification test comparing the fixed-effects
#' (within) and random-effects GLS estimators.
#'
#' @param formula A linear formula.
#' @param data A data frame containing the formula variables.
#' @param id_var Name of the cross-sectional identifier.
#' @param time_var Name of the time identifier.
#'
#' @return A list with elements \code{statistic} and \code{p.value}, or
#'   \code{NULL} if the test cannot be computed.
#' @export
hausman_plm <- function(formula, data, id_var = "id", time_var = "year") {
  if (!requireNamespace("plm", quietly = TRUE)) return(NULL)
  pdf <- tryCatch(
    plm::pdata.frame(as.data.frame(data), index = c(id_var, time_var)),
    error = function(e) NULL
  )
  if (is.null(pdf)) return(NULL)
  fe <- tryCatch(plm::plm(formula, data = pdf, model = "within"),
                 error = function(e) NULL)
  re <- tryCatch(plm::plm(formula, data = pdf, model = "random"),
                 error = function(e) NULL)
  if (is.null(fe) || is.null(re)) return(NULL)
  h <- tryCatch(plm::phtest(fe, re), error = function(e) NULL)
  if (is.null(h)) return(NULL)
  list(
    statistic = unname(h$statistic),
    p.value   = h$p.value
  )
}

#' Compare BPVS-FE and BPVS-RE
#'
#' Builds a combined table of FE and RE estimates and PIPs, appends the
#' Bayesian Hausman results, and computes the frequentist Hausman test.
#'
#' @inheritParams bayesian_hausman
#' @param formula A linear formula (needed for the frequentist Hausman test).
#' @param id_var Name of the cross-sectional identifier column.
#' @param time_var Name of the time identifier column.
#'
#' @return A list with elements \code{comparison}, \code{stats_fe},
#'   \code{stats_re}, \code{bayesian_hausman}, and \code{hausman_plm}.
#' @export
compare_bpvs <- function(fit_fe, fit_re, data, outcome_var = "Income",
                         formula = NULL, id_var = "id", time_var = "year",
                         effect_threshold = 0.01) {

  fe_pip <- compute_pip(fit_fe, data, outcome_var, effect_threshold)
  re_pip <- compute_pip(fit_re, data, outcome_var, effect_threshold)

  fe_rank <- rank_variables(fe_pip) %>%
    select(Variable, FE_Estimate = Estimate, FE_PIP = PIP, FE_Selection = Selection)
  re_rank <- rank_variables(re_pip) %>%
    select(Variable, RE_Estimate = Estimate, RE_PIP = PIP, RE_Selection = Selection)

  comparison <- merge(fe_rank, re_rank, by = "Variable", all = TRUE) %>%
    arrange(desc(FE_PIP))

  stats_fe <- extract_model_stats(fit_fe)
  stats_re <- extract_model_stats(fit_re)

  bh <- bayesian_hausman(fit_fe, fit_re, data, outcome_var, effect_threshold)
  comparison <- merge(comparison, bh, by = "Variable", all.x = TRUE)

  haus <- hausman_plm(formula, data, id_var, time_var)

  list(
    comparison = comparison,
    stats_fe   = stats_fe,
    stats_re   = stats_re,
    bayesian_hausman = bh,
    hausman_plm = haus
  )
}

#' Print the FE vs RE comparison
#'
#' Prints model fit, the comparison table, the Bayesian Hausman table, and the
#' frequentist Hausman test result to the console.
#'
#' @param cmp An object returned by \code{\link{compare_bpvs}}.
#' @return Invisibly returns \code{NULL}; called for its console output.
#' @export
print_compare <- function(cmp) {
  cat("\n============================================================\n")
  cat("BPVS -- FE vs RE COMPARISON\n")
  cat("============================================================\n")

  cat("\n--- Model Fit ---\n")
  cat(sprintf("  FE: WAIC = %.2f (se %.2f), Bayes R2 = %.4f\n",
              cmp$stats_fe$WAIC, cmp$stats_fe$WAIC_SE, cmp$stats_fe$Bayes_R2))
  cat(sprintf("  RE: WAIC = %.2f (se %.2f), Bayes R2 = %.4f\n",
              cmp$stats_re$WAIC, cmp$stats_re$WAIC_SE, cmp$stats_re$Bayes_R2))

  cat("\n--- FE vs RE: Estimates and PIPs ---\n")
  print(cmp$comparison %>%
          select(Variable, FE_Estimate, RE_Estimate, FE_PIP, RE_PIP,
                 FE_Selection, RE_Selection),
        row.names = FALSE)

  cat("\n--- Bayesian Hausman (posterior of beta_FE - beta_RE) ---\n")
  print(cmp$bayesian_hausman %>%
          select(Variable, Delta_median, Delta_q025, Delta_q975,
                 Pr_delta_not0, RE_Consistent),
        row.names = FALSE)
  n_biased <- sum(cmp$bayesian_hausman$RE_Consistent == "No (RE biased)")
  if (n_biased > 0) {
    cat(sprintf("  %d variable(s) show FE/RE disagreement -> RE potentially biased there.\n", n_biased))
  } else {
    cat("  No FE/RE disagreement -> RE consistent with FE.\n")
  }

  if (!is.null(cmp$hausman_plm)) {
    cat("\n--- Frequentist Hausman (plm, auxiliary) ---\n")
    cat(sprintf("  Statistic = %.3f, p-value = %.4f\n",
                cmp$hausman_plm$statistic, cmp$hausman_plm$p.value))
    if (cmp$hausman_plm$p.value < 0.05) {
      cat("  p < 0.05 -> RE inconsistent; prefer BPVS-FE.\n")
    } else {
      cat("  p >= 0.05 -> no evidence against RE; BPVS-RE acceptable.\n")
    }
  }
  cat("============================================================\n")
}

#' Fit, rank, and compare BPVS models (pipeline)
#'
#' One-call driver: fits BPVS-FE and/or BPVS-RE, computes posterior inclusion
#' probabilities, ranks the predictors, extracts fit statistics, prints the
#' results, and (when both specifications are fitted) computes the Bayesian and
#' frequentist Hausman comparison.
#'
#' @inheritParams fit_bpvs_fe
#' @param effects Character vector of specifications to fit; any subset of
#'   \code{c("fe", "re")}.
#' @param outcome_var Name of the outcome column.
#' @param time_var Name of the time identifier (for the frequentist Hausman
#'   test).
#'
#' @return An invisible list with elements (as available):
#'   \code{fit_fe}, \code{fit_re}, \code{rank_fe}, \code{rank_re},
#'   \code{pip_fe}, \code{pip_re}, \code{stats_fe}, \code{stats_re}, and
#'   \code{comparison}.
#' @export
#' @examples
#' \dontrun{
#' # Quick demo on synthetic data (small MCMC, for illustration only)
#' set.seed(456)
#' dat <- simulate_panel_data(N_units = 20, T_periods = 6,
#'                            beta_true = c(1.5, -0.8, 0),
#'                            sigma = 0.8, tau_unit = 0.4)
#' fmla <- Income ~ FDI + INF + POP
#' res <- bpvs_pipeline(fmla, dat, effects = c("fe", "re"),
#'                      outcome_var = "Income", chains = 2,
#'                      iter = 500, warmup = 250, cores = 1, seed = 1)
#' }
bpvs_pipeline <- function(formula, data, effects = c("fe", "re"),
                          outcome_var = "Income", id_var = "id",
                          time_var = "year", chains = 4, iter = 4000,
                          warmup = 2000, cores = 4, seed = 123,
                          adapt_delta = 0.95, max_treedepth = 12, ...) {

  effects <- match.arg(effects, several.ok = TRUE)

  results <- list()

  if ("fe" %in% effects) {
    cat("\n--- Fitting BPVS-FE (fixed effects, unit dummies) ---\n")
    fit_fe <- fit_bpvs_fe(formula, data, id_var = id_var, chains = chains,
                          iter = iter, warmup = warmup, cores = cores,
                          seed = seed, adapt_delta = adapt_delta,
                          max_treedepth = max_treedepth, ...)
    fe_pip  <- compute_pip(fit_fe, data, outcome_var)
    fe_rank <- rank_variables(fe_pip)
    fe_stat <- extract_model_stats(fit_fe)
    print_bpvs_results(fe_rank, fe_stat, title = "BPVS-FE")
    results$fit_fe   <- fit_fe
    results$rank_fe  <- fe_rank
    results$pip_fe   <- fe_pip
    results$stats_fe <- fe_stat
  }

  if ("re" %in% effects) {
    cat("\n--- Fitting BPVS-RE (random intercepts) ---\n")
    fit_re <- fit_bpvs_re(formula, data, id_var = id_var, chains = chains,
                          iter = iter, warmup = warmup, cores = cores,
                          seed = seed, adapt_delta = adapt_delta,
                          max_treedepth = max_treedepth, ...)
    re_pip  <- compute_pip(fit_re, data, outcome_var)
    re_rank <- rank_variables(re_pip)
    re_stat <- extract_model_stats(fit_re)
    print_bpvs_results(re_rank, re_stat, title = "BPVS-RE")
    results$fit_re   <- fit_re
    results$rank_re  <- re_rank
    results$pip_re   <- re_pip
    results$stats_re <- re_stat
  }

  if (all(c("fe", "re") %in% effects)) {
    cmp <- compare_bpvs(results$fit_fe, results$fit_re, data, outcome_var,
                        formula = formula, id_var = id_var, time_var = time_var)
    print_compare(cmp)
    results$comparison <- cmp
  }

  invisible(results)
}
