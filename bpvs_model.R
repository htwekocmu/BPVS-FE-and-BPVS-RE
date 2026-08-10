# ======================================================================
# bpvs_model.R — Bayesian Panel Variable Selection (BPVS) Model
# Pathairat Pastpipatkul and Htwe Ko (2025)
# Revised: FE (unit dummies) and RE (random intercept) variants
# ======================================================================
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
#   Prior on slopes (Regularized Horseshoe — continuous spike-and-slab):
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
#
# OUTPUTS:
#   fit_bpvs_fe()        — fits BPVS-FE via brms (unit dummies + horseshoe)
#   fit_bpvs_re()        — fits BPVS-RE via brms (random intercept + horseshoe)
#   fit_bpvs()           — dispatcher: effects = "fe" | "re"
#   compute_pip()        — computes posterior inclusion probabilities
#   rank_variables()     — ranks predictors by PIP with selection labels
#   extract_model_stats()— WAIC and Bayesian R2
#   compare_bpvs()       — FE vs RE comparison + Bayesian/Frequentist Hausman
#   bpvs_pipeline()      — fit + compute + rank + print in one call
# ======================================================================

required_packages <- c("brms", "tidyverse", "posterior", "bayesplot", "loo", "plm")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

backend <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"

# -----------------------------------------------------------------------
# extract_slope_terms: pull predictor names from a linear formula RHS
# -----------------------------------------------------------------------
extract_slope_terms <- function(formula) {
  t <- terms(formula)
  attr(t, "term.labels")
}

# -----------------------------------------------------------------------
# fit_bpvs_fe: Fit BPVS with FIXED EFFECTS (unit dummies)
#   y ~ alpha + eta ;  alpha ~ 0 + id ;  eta ~ 0 + x1 + ... + xp
#   Horseshoe applied ONLY to slope block (eta); dummies unshrunk.
# -----------------------------------------------------------------------
fit_bpvs_fe <- function(formula, data, id_var = "id",
                        chains = 4, iter = 4000, warmup = 2000, cores = 4,
                        seed = 123, adapt_delta = 0.95, max_treedepth = 12,
                        refresh = 500, ...) {

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

# -----------------------------------------------------------------------
# fit_bpvs_re: Fit BPVS with RANDOM EFFECTS (random intercepts)
# -----------------------------------------------------------------------
fit_bpvs_re <- function(formula, data, id_var = "id",
                        chains = 4, iter = 4000, warmup = 2000, cores = 4,
                        seed = 123, adapt_delta = 0.95, max_treedepth = 12,
                        refresh = 500, ...) {

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

# -----------------------------------------------------------------------
# fit_bpvs: dispatcher
# -----------------------------------------------------------------------
fit_bpvs <- function(formula, data, effects = c("re", "fe"), id_var = "id", ...) {
  effects <- match.arg(effects)
  if (effects == "fe") {
    fit_bpvs_fe(formula, data, id_var = id_var, ...)
  } else {
    fit_bpvs_re(formula, data, id_var = id_var, ...)
  }
}

# -----------------------------------------------------------------------
# compute_pip: Compute Posterior Inclusion Probabilities
#   PIP_j = Pr(|beta_j| > eps | data), eps = effect_threshold * sd(y)
#   FE: slope block columns are b_eta_<var>
#   RE: slope columns are b_<var>
# -----------------------------------------------------------------------
compute_pip <- function(fit, data, outcome_var = "Income",
                        effect_threshold = 0.01) {
  post <- as_draws_df(fit)
  eps <- effect_threshold * sd(data[[outcome_var]], na.rm = TRUE)

  effects <- attr(fit, "bpvs_effects") %||% "re"
  pattern <- if (effects == "fe") "^b_eta_" else "^b_"

  b_names <- grep(pattern, names(post), value = TRUE)
  # RE: drop intercept; FE: eta block has no intercept
  b_names <- setdiff(b_names, "b_Intercept")

  pip <- sapply(b_names, function(col) {
    mean(abs(post[[col]]) > eps, na.rm = TRUE)
  })
  med   <- sapply(b_names, function(col) median(post[[col]], na.rm = TRUE))
  q025  <- sapply(b_names, function(col) quantile(post[[col]], 0.025, na.rm = TRUE))
  q975  <- sapply(b_names, function(col) quantile(post[[col]], 0.975, na.rm = TRUE))

  # Clean variable names: strip "b_eta_" / "b_" prefix
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

# -----------------------------------------------------------------------
# rank_variables: Rank predictors by PIP with selection labels
# -----------------------------------------------------------------------
rank_variables <- function(pip_df) {
  pip_df %>%
    arrange(desc(PIP)) %>%
    mutate(
      Rank      = 1:n(),
      Selection = case_when(
        PIP > 0.5  ~ "Selected",
        PIP > 0.1  ~ "Weak",
        TRUE       ~ "Excluded"
      )
    )
}

# -----------------------------------------------------------------------
# extract_model_stats: Extract WAIC and Bayesian R2
# -----------------------------------------------------------------------
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

# -----------------------------------------------------------------------
# print_bpvs_results: Pretty-print BPVS output
# -----------------------------------------------------------------------
print_bpvs_results <- function(ranking, stats, title = "BPVS") {
  cat("\n============================================================\n")
  cat(paste(title, "— BAYESIAN PANEL VARIABLE SELECTION RESULTS\n"))
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

# -----------------------------------------------------------------------
# bayesian_hausman: compare FE vs RE slope posteriors
#   For each variable, posterior of (beta_FE - beta_RE);
#   if the 95% interval excludes 0 -> RE inconsistent for that variable.
# -----------------------------------------------------------------------
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

# -----------------------------------------------------------------------
# hausman_plm: Frequentist Hausman test via plm (auxiliary)
# -----------------------------------------------------------------------
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

# -----------------------------------------------------------------------
# compare_bpvs: FE vs RE comparison table + Hausman checks
# -----------------------------------------------------------------------
compare_bpvs <- function(fit_fe, fit_re, data, outcome_var = "Income",
                         formula = NULL, id_var = "id", time_var = "year",
                         effect_threshold = 0.01) {

  fe_pip <- compute_pip(fit_fe, data, outcome_var, effect_threshold)
  re_pip <- compute_pip(fit_re, data, outcome_var, effect_threshold)

  fe_rank <- rank_variables(fe_pip) %>% select(Variable, FE_Estimate = Estimate, FE_PIP = PIP,
                                               FE_Selection = Selection)
  re_rank <- rank_variables(re_pip) %>% select(Variable, RE_Estimate = Estimate, RE_PIP = PIP,
                                               RE_Selection = Selection)

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

# -----------------------------------------------------------------------
# print_compare: print FE vs RE comparison
# -----------------------------------------------------------------------
print_compare <- function(cmp) {
  cat("\n============================================================\n")
  cat("BPVS — FE vs RE COMPARISON\n")
  cat("============================================================\n")

  cat("\n--- Model Fit ---\n")
  cat(sprintf("  FE: WAIC = %.2f (se %.2f), Bayes R2 = %.4f\n",
              cmp$stats_fe$WAIC, cmp$stats_fe$WAIC_SE, cmp$stats_fe$Bayes_R2))
  cat(sprintf("  RE: WAIC = %.2f (se %.2f), Bayes R2 = %.4f\n",
              cmp$stats_re$WAIC, cmp$stats_re$WAIC_SE, cmp$stats_re$Bayes_R2))

  cat("\n--- FE vs RE: Estimates and PIPs ---\n")
  print(cmp$comparison %>% select(Variable, FE_Estimate, RE_Estimate,
                                  FE_PIP, RE_PIP, FE_Selection, RE_Selection),
        row.names = FALSE)

  cat("\n--- Bayesian Hausman (posterior of beta_FE - beta_RE) ---\n")
  print(cmp$bayesian_hausman %>% select(Variable, Delta_median, Delta_q025,
                                        Delta_q975, Pr_delta_not0, RE_Consistent),
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

# -----------------------------------------------------------------------
# bpvs_pipeline: fit (FE and/or RE) + compute + rank + print
# -----------------------------------------------------------------------
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

`%||%` <- function(a, b) if (is.null(a)) b else a

