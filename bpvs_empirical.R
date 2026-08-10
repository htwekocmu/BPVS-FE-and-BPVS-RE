# ======================================================================
# bpvs_empirical.R — BPVS Real Data Application (FE and RE)
# Pathairat Pastpipatkul and Htwe Ko (2025)
# ======================================================================
#
# USAGE:
#   source("bpvs_model.R")
#   source("bpvs_empirical.R")
#
# OUTPUTS:
#   Console — FE & RE rankings, PIPs, selection, comparison + Hausman
#   bpvs_emp_results/ — PNG figures + saved models (.rds)
# ======================================================================

cat("============================================================\n")
cat("BPVS — EMPIRICAL APPLICATION (FE and RE)\n")
cat("============================================================\n")
cat("\nSources: bpvs_model.R\n\n")

source("bpvs_model.R")
library(patchwork)
library(ggrepel)

set.seed(2025123)

# -------------------------------------------------------------------
# 1. Load and prepare data
# -------------------------------------------------------------------

# --- EDIT THIS PATH to point to your CSV file ---
data_path <- "C:/Users/User/OneDrive/Desktop/developing.csv"

# --- EDIT THIS if your outcome column has a different name ---
outcome_var <- "Growth"
time_var    <- "year"
id_var      <- "id"

if (!file.exists(data_path)) {
  stop("Data file not found: ", data_path, "\n",
       "Please update 'data_path' in bpvs_empirical.R to point to your CSV.")
}

data_raw <- read.csv(data_path)
cat(sprintf("Loaded: %s\n", data_path))
cat(sprintf("Dimensions: %d rows x %d columns\n", nrow(data_raw), ncol(data_raw)))

# Drop character identifiers, ensure panel structure
data <- data_raw %>%
  dplyr::select(-any_of("Country")) %>%
  mutate(!!sym(id_var) := factor(.data[[id_var]]))

cat(sprintf("Panel: %d countries, %d time periods, %d total observations\n",
            n_distinct(data[[id_var]]), n_distinct(data[[time_var]]), nrow(data)))

if (!outcome_var %in% names(data)) {
  stop("Outcome variable '", outcome_var, "' not found. Adjust 'outcome_var'.")
}

# Identify numeric predictors (exclude id, year, outcome)
predictors <- setdiff(names(data), c(id_var, time_var, outcome_var))
cat("Predictors:", paste(predictors, collapse = ", "), "\n")

# Scale predictors for horseshoe prior
data_scaled <- data
for (v in predictors) {
  m <- mean(data_scaled[[v]], na.rm = TRUE)
  s <- sd(data_scaled[[v]], na.rm = TRUE)
  data_scaled[[v]] <- (data_scaled[[v]] - m) / s
}

# -------------------------------------------------------------------
# 2. Build formula
# -------------------------------------------------------------------
formula <- as.formula(paste(outcome_var, "~", paste(predictors, collapse = " + ")))
cat("Model formula:\n")
print(formula)

# -------------------------------------------------------------------
# 3. Fit BPVS-FE and BPVS-RE
# -------------------------------------------------------------------
cat("\nFitting BPVS-FE and BPVS-RE (horseshoe prior)...\n")

res <- bpvs_pipeline(
  formula       = formula,
  data          = data_scaled,
  effects       = c("fe", "re"),
  outcome_var   = outcome_var,
  id_var        = id_var,
  time_var      = time_var,
  chains        = 4,
  iter          = 4000,
  warmup        = 2000,
  cores         = 4,
  seed          = 789,
  adapt_delta   = 0.95,
  max_treedepth = 12,
  refresh       = 1000
)

fit_fe  <- res$fit_fe
fit_re  <- res$fit_re
rank_fe <- res$rank_fe
rank_re <- res$rank_re
cmp     <- res$comparison

# -------------------------------------------------------------------
# 4. Save models and results
# -------------------------------------------------------------------
output_dir <- "bpvs_emp_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

saveRDS(fit_fe, file.path(output_dir, "bpvs_fit_fe.rds"))
saveRDS(fit_re, file.path(output_dir, "bpvs_fit_re.rds"))
write.csv(rank_fe, file.path(output_dir, "bpvs_ranking_fe.csv"), row.names = FALSE)
write.csv(rank_re, file.path(output_dir, "bpvs_ranking_re.csv"), row.names = FALSE)
write.csv(cmp$comparison, file.path(output_dir, "bpvs_comparison.csv"), row.names = FALSE)
write.csv(cmp$bayesian_hausman, file.path(output_dir, "bpvs_bayes_hausman.csv"), row.names = FALSE)
if (!is.null(cmp$hausman_plm)) {
  write.csv(as.data.frame(cmp$hausman_plm), file.path(output_dir, "bpvs_hausman_plm.csv"),
            row.names = FALSE)
}

# -------------------------------------------------------------------
# 5. Visualizations
# -------------------------------------------------------------------

# 5a. Inclusion probability bar chart (FE vs RE)
make_inclusion_plot <- function(ranking, title, filename) {
  ranking$Variable <- factor(ranking$Variable, levels = rev(ranking$Variable))
  p <- ggplot(ranking, aes(x = Variable, y = PIP)) +
    geom_segment(aes(xend = Variable, yend = 0), color = "gray60", linewidth = 0.5) +
    geom_point(aes(color = Selection), size = 4) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_hline(yintercept = 0.1, linetype = "dotted", color = "orange", alpha = 0.5) +
    scale_color_manual(values = c("Selected" = "#1b9e77", "Weak" = "#d95f02",
                                  "Excluded" = "#7570b3")) +
    labs(title = title,
         subtitle = "Dashed red = 0.5 (Selected), dotted orange = 0.1 (Weak signal)",
         x = "Predictor", y = "PIP = Pr(|beta| > eps | data)",
         color = "Decision") +
    ylim(0, 1.05) + coord_flip() +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(file.path(output_dir, filename), p, width = 8, height = max(4, 0.5 * nrow(ranking)))
  p
}

p_incl_fe <- make_inclusion_plot(rank_fe, "BPVS-FE: Posterior Inclusion Probabilities",
                                 "inclusion_probabilities_FE.png")
p_incl_re <- make_inclusion_plot(rank_re, "BPVS-RE: Posterior Inclusion Probabilities",
                                 "inclusion_probabilities_RE.png")

# 5b. Coefficient plot with 95% CIs
make_coef_plot <- function(ranking, title, filename) {
  ranking$Variable <- factor(ranking$Variable, levels = ranking$Variable)
  p <- ggplot(ranking, aes(x = Estimate, y = Variable)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper, color = Selection),
                    size = 0.9, fatten = 2.5) +
    scale_color_manual(values = c("Selected" = "#1b9e77", "Weak" = "#d95f02",
                                  "Excluded" = "#7570b3")) +
    labs(title = title,
         subtitle = "Horseshoe prior shrinks irrelevant coefficients toward zero",
         x = "Posterior Median (95% CI)", y = "Predictor", color = "Decision") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(file.path(output_dir, filename), p, width = 8, height = max(4, 0.6 * nrow(ranking)))
  p
}

p_coef_fe <- make_coef_plot(rank_fe, "BPVS-FE: Coefficient Estimates (95% CI)",
                            "coefficient_estimates_FE.png")
p_coef_re <- make_coef_plot(rank_re, "BPVS-RE: Coefficient Estimates (95% CI)",
                            "coefficient_estimates_RE.png")

# 5c. PIP comparison FE vs RE
pip_compare <- merge(
  rank_fe %>% select(Variable, PIP_FE = PIP),
  rank_re %>% select(Variable, PIP_RE = PIP),
  by = "Variable", all = TRUE
)
p_pip_compare <- ggplot(pip_compare, aes(x = PIP_FE, y = PIP_RE, label = Variable)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 4, color = "steelblue", alpha = 0.8) +
  geom_text_repel(size = 3.5) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "BPVS: Posterior Inclusion Probabilities — FE vs RE",
       x = "BPVS-FE PIP", y = "BPVS-RE PIP") +
  theme_minimal()
ggsave(file.path(output_dir, "pip_FE_vs_RE.png"), p_pip_compare, width = 6, height = 6)

# 5d. Trace plots for convergence (slopes)
beta_params_fe <- grep("^b_eta_", variables(fit_fe), value = TRUE)
if (length(beta_params_fe) > 12) beta_params_fe <- beta_params_fe[1:12]
p_trace_fe <- mcmc_trace(fit_fe, pars = beta_params_fe, facet_args = list(ncol = 3))
ggsave(file.path(output_dir, "trace_plots_FE.png"), p_trace_fe,
       width = 12, height = 3 * ceiling(length(beta_params_fe) / 3))

beta_params_re <- setdiff(grep("^b_", variables(fit_re), value = TRUE), "b_Intercept")
if (length(beta_params_re) > 12) beta_params_re <- beta_params_re[1:12]
p_trace_re <- mcmc_trace(fit_re, pars = beta_params_re, facet_args = list(ncol = 3))
ggsave(file.path(output_dir, "trace_plots_RE.png"), p_trace_re,
       width = 12, height = 3 * ceiling(length(beta_params_re) / 3))

# Display
print(p_incl_fe)
print(p_incl_re)
print(p_coef_fe)
print(p_coef_re)
print(p_pip_compare)

cat(sprintf("\nResults saved to: %s/\n", output_dir))
cat("============================================================\n")
cat("EMPIRICAL APPLICATION COMPLETE\n")
cat("============================================================\n")
