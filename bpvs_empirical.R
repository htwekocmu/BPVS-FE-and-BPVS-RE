# ======================================================================
# bpvs_empirical.R -- BPVS real-data application (FE and RE)
# Pathairat Pastpipatkul and Htwe Ko (2025)
# ======================================================================

#' Prepare a panel data frame for BPVS estimation
#'
#' Internal helper that drops non-numeric identifier columns, converts the unit
#' identifier to a factor, and scales all numeric predictors to zero mean and
#' unit variance (important for the horseshoe auto-scaling).
#'
#' @param data_raw A data frame.
#' @param outcome_var Name of the outcome column.
#' @param id_var Name of the unit identifier column.
#' @param time_var Name of the time identifier column.
#' @param drop_cols Character vector of additional columns to drop (e.g.
#'   country names).
#'
#' @return A list with elements \code{data} (the scaled data frame),
#'   \code{predictors} (character vector of scaled predictor names), and
#'   \code{formula} (the model formula).
#' @keywords internal
#' @noRd
prepare_panel_data <- function(data_raw, outcome_var = "Growth",
                               id_var = "id", time_var = "year",
                               drop_cols = NULL) {
  if (!is.data.frame(data_raw)) data_raw <- as.data.frame(data_raw)

  data <- data_raw %>%
    dplyr::select(-dplyr::any_of(c("Country", drop_cols))) %>%
    dplyr::mutate(!!rlang::sym(id_var) := factor(.data[[id_var]]))

  if (!outcome_var %in% names(data)) {
    stop("Outcome variable '", outcome_var, "' not found in data.")
  }

  predictors <- setdiff(names(data), c(id_var, time_var, outcome_var))
  if (length(predictors) == 0) stop("No predictors found in data.")

  for (v in predictors) {
    m <- mean(data[[v]], na.rm = TRUE)
    s <- stats::sd(data[[v]], na.rm = TRUE)
    data[[v]] <- (data[[v]] - m) / s
  }

  formula <- as.formula(paste(outcome_var, "~", paste(predictors, collapse = " + ")))

  list(data = data, predictors = predictors, formula = formula)
}

#' Run the BPVS empirical application end-to-end
#'
#' Mirrors the script \file{bpvs_empirical.R}: loads and scales the panel data,
#' fits BPVS-FE and BPVS-RE, computes PIPs and rankings, saves the fitted
#' models and result tables, and writes PNG figures (inclusion probabilities,
#' coefficient estimates, FE-vs-RE PIP scatter) to an output directory.
#'
#' @param data A data frame (alternative to \code{data_path}).
#' @param data_path Path to a CSV file containing the panel data. Ignored if
#'   \code{data} is supplied.
#' @param outcome_var Name of the outcome column.
#' @param id_var Name of the unit identifier column.
#' @param time_var Name of the time identifier column.
#' @param output_dir Directory where fitted models, tables, and figures are
#'   saved.
#' @param effects Which specifications to fit, any subset of \code{c("fe", "re")}.
#' @param ... Additional arguments passed to \code{\link{bpvs_pipeline}} (e.g.
#'   \code{chains}, \code{iter}, \code{warmup}, \code{seed},
#'   \code{adapt_delta}).
#'
#' @return Invisibly, the results list from \code{\link{bpvs_pipeline}}.
#' @export
bpvs_empirical <- function(data = NULL, data_path = NULL,
                           outcome_var = "Growth", id_var = "id",
                           time_var = "year",
                           output_dir = "bpvs_emp_results",
                           effects = c("fe", "re"), ...) {

  if (is.null(data)) {
    if (is.null(data_path) || !file.exists(data_path)) {
      stop("Provide either 'data' or an existing 'data_path'.")
    }
    data <- utils::read.csv(data_path)
    cat(sprintf("Loaded: %s\n", data_path))
    cat(sprintf("Dimensions: %d rows x %d columns\n", nrow(data), ncol(data)))
  }

  prep <- prepare_panel_data(data, outcome_var, id_var, time_var)
  data_scaled <- prep$data
  predictors  <- prep$predictors
  formula     <- prep$formula

  cat(sprintf("Panel: %d units, %d time periods, %d observations\n",
              dplyr::n_distinct(data_scaled[[id_var]]),
              dplyr::n_distinct(data_scaled[[time_var]]), nrow(data_scaled)))
  cat("Predictors:", paste(predictors, collapse = ", "), "\n")
  cat("Model formula:\n"); print(formula)

  res <- bpvs_pipeline(formula = formula, data = data_scaled,
                       effects = effects, outcome_var = outcome_var,
                       id_var = id_var, time_var = time_var, ...)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (!is.null(res$fit_fe)) {
    saveRDS(res$fit_fe, file.path(output_dir, "bpvs_fit_fe.rds"))
    write.csv(res$rank_fe, file.path(output_dir, "bpvs_ranking_fe.csv"), row.names = FALSE)
    p1 <- plot_inclusion(res$rank_fe, "BPVS-FE: Posterior Inclusion Probabilities")
    ggplot2::ggsave(file.path(output_dir, "inclusion_probabilities_FE.png"), p1,
                    width = 8, height = max(4, 0.5 * nrow(res$rank_fe)))
    p2 <- plot_coefficients(res$rank_fe, "BPVS-FE: Coefficient Estimates (95% CI)")
    ggplot2::ggsave(file.path(output_dir, "coefficient_estimates_FE.png"), p2,
                    width = 8, height = max(4, 0.6 * nrow(res$rank_fe)))
  }
  if (!is.null(res$fit_re)) {
    saveRDS(res$fit_re, file.path(output_dir, "bpvs_fit_re.rds"))
    write.csv(res$rank_re, file.path(output_dir, "bpvs_ranking_re.csv"), row.names = FALSE)
    p3 <- plot_inclusion(res$rank_re, "BPVS-RE: Posterior Inclusion Probabilities")
    ggplot2::ggsave(file.path(output_dir, "inclusion_probabilities_RE.png"), p3,
                    width = 8, height = max(4, 0.5 * nrow(res$rank_re)))
    p4 <- plot_coefficients(res$rank_re, "BPVS-RE: Coefficient Estimates (95% CI)")
    ggplot2::ggsave(file.path(output_dir, "coefficient_estimates_RE.png"), p4,
                    width = 8, height = max(4, 0.6 * nrow(res$rank_re)))
  }

  if (!is.null(res$comparison)) {
    write.csv(res$comparison$comparison,
              file.path(output_dir, "bpvs_comparison.csv"), row.names = FALSE)
    write.csv(res$comparison$bayesian_hausman,
              file.path(output_dir, "bpvs_bayes_hausman.csv"), row.names = FALSE)
    if (!is.null(res$comparison$hausman_plm)) {
      write.csv(as.data.frame(res$comparison$hausman_plm),
                file.path(output_dir, "bpvs_hausman_plm.csv"), row.names = FALSE)
    }
    p5 <- plot_pip_compare(res$rank_fe, res$rank_re)
    ggplot2::ggsave(file.path(output_dir, "pip_FE_vs_RE.png"), p5, width = 6, height = 6)
  }

  cat(sprintf("\nResults saved to: %s/\n", output_dir))
  invisible(res)
}
