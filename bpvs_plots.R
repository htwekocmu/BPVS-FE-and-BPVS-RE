# ======================================================================
# bpvs_plots.R -- visualizations for BPVS results
# Pathairat Pastpipatkul and Htwe Ko (2025)
# ======================================================================

#' Inclusion probability bar chart
#'
#' Draws a lollipop chart of posterior inclusion probabilities with reference
#' lines at \code{0.5} (Selected) and \code{0.1} (Weak signal).
#'
#' @param ranking A ranking table from \code{\link{rank_variables}}.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#'
#' @return A \pkg{ggplot2} object.
#' @export
plot_inclusion <- function(ranking, title = "Posterior Inclusion Probabilities",
                           subtitle = "Dashed red = 0.5 (Selected), dotted orange = 0.1 (Weak signal)") {
  ranking$Variable <- factor(ranking$Variable, levels = rev(ranking$Variable))
  ggplot(ranking, aes(x = Variable, y = PIP)) +
    geom_segment(aes(xend = Variable, yend = 0), color = "gray60", linewidth = 0.5) +
    geom_point(aes(color = Selection), size = 4) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_hline(yintercept = 0.1, linetype = "dotted", color = "orange", alpha = 0.5) +
    scale_color_manual(values = c("Selected" = "#1b9e77", "Weak" = "#d95f02",
                                  "Excluded" = "#7570b3")) +
    labs(title = title, subtitle = subtitle,
         x = "Predictor", y = "PIP = Pr(|beta| > eps | data)",
         color = "Decision") +
    ylim(0, 1.05) + coord_flip() +
    theme_minimal() + theme(legend.position = "bottom")
}

#' Coefficient plot with 95 percent credible intervals
#'
#' Draws a forest plot of posterior medians with 95\% credible intervals,
#' colored by the selection decision.
#'
#' @param ranking A ranking table from \code{\link{rank_variables}}.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#'
#' @return A \pkg{ggplot2} object.
#' @export
plot_coefficients <- function(ranking, title = "Coefficient Estimates (95% CI)",
                              subtitle = "Horseshoe prior shrinks irrelevant coefficients toward zero") {
  ranking$Variable <- factor(ranking$Variable, levels = ranking$Variable)
  ggplot(ranking, aes(x = Estimate, y = Variable)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper, color = Selection),
                    size = 0.9, fatten = 2.5) +
    scale_color_manual(values = c("Selected" = "#1b9e77", "Weak" = "#d95f02",
                                  "Excluded" = "#7570b3")) +
    labs(title = title, subtitle = subtitle,
         x = "Posterior Median (95% CI)", y = "Predictor", color = "Decision") +
    theme_minimal() + theme(legend.position = "bottom")
}

#' PIP scatterplot: FE vs RE
#'
#' Plots BPVS-FE inclusion probabilities against BPVS-RE inclusion
#' probabilities, with variable labels. Points near the diagonal indicate that
#' the two specifications agree on the inclusion status of a variable.
#'
#' @param rank_fe A ranking table for BPVS-FE.
#' @param rank_re A ranking table for BPVS-RE.
#' @param title Plot title.
#'
#' @return A \pkg{ggplot2} object.
#' @export
plot_pip_compare <- function(rank_fe, rank_re, title = "Posterior Inclusion Probabilities -- FE vs RE") {
  pip_compare <- merge(
    rank_fe %>% select(Variable, PIP_FE = PIP),
    rank_re %>% select(Variable, PIP_RE = PIP),
    by = "Variable", all = TRUE
  )
  ggplot(pip_compare, aes(x = PIP_FE, y = PIP_RE, label = Variable)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 4, color = "steelblue", alpha = 0.8) +
    ggrepel::geom_text_repel(size = 3.5) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title, x = "BPVS-FE PIP", y = "BPVS-RE PIP") +
    theme_minimal()
}

#' Parameter recovery plot
#'
#' Scatter of true against estimated coefficients with 95\% credible intervals
#' and the 45-degree line, used to visualize simulation recovery.
#'
#' @param summ A recovery summary from \code{\link{simulate_panel_data}} style
#'   output (columns \code{True}, \code{Estimated}, \code{CI_lower},
#'   \code{CI_upper}, \code{Variable}).
#' @param title Plot title.
#'
#' @return A \pkg{ggplot2} object.
#' @export
plot_recovery <- function(summ, title = "Parameter Recovery") {
  summ$Recovered <- ifelse(abs(summ$True) > 0, "Non-zero", "Zero")
  ggplot(summ, aes(x = True, y = Estimated, color = Recovered)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper), size = 0.8) +
    geom_point(size = 3.5) +
    ggrepel::geom_text_repel(aes(label = Variable), size = 3.5) +
    scale_color_manual(values = c("Non-zero" = "#d95f02", "Zero" = "#1b9e77")) +
    coord_fixed(xlim = c(-1.5, 2), ylim = c(-1.5, 2)) +
    labs(title = title,
         x = "True Coefficient", y = "Posterior Median (95% CI)",
         color = "True Status") +
    theme_minimal() + theme(legend.position = "bottom")
}
