# ======================================================================
# bpvs_shiny.R -- Shiny app launcher for BPVS
# Pathairat Pastpipatkul and Htwe Ko (2025)
# ======================================================================

#' Launch the BPVS Shiny application
#'
#' Starts the interactive Shiny app bundled with the package. The app lets the
#' user upload a panel data set (or use the bundled \code{\link{developing}}
#' example), choose the outcome / id / time columns, select the specifications
#' (FE and/or RE), and explore the rankings, PIPs, coefficient plots, and
#' Hausman diagnostics.
#'
#' The app is stored in \code{inst/shiny/BPVS} and is located via
#' \code{system.file("shiny", "BPVS", package = "BPVS")}.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @return The return value of \code{shiny::runApp()}; called for its side
#'   effect of launching the app.
#' @export
run_bpvs_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The 'shiny' package is required to run the BPVS app. ",
         "Install it with install.packages('shiny').")
  }
  app_dir <- system.file("shiny", "BPVS", package = "BPVS")
  if (!nzchar(app_dir)) {
    stop("Could not locate the BPVS Shiny app. Reinstall the package with ",
         "the app directory intact.")
  }
  shiny::runApp(app_dir, ...)
}
