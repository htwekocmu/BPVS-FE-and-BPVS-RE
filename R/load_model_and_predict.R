#' Load BPVS model (FE or RE)
#'
#' Load the saved BPVS model for Fixed Effects (FE) or Random Effects (RE).
#'
#' @param type Character, one of "FE" or "RE". Which model to load.
#' @param model_path Optional explicit path to a .rds file. If NULL the function will try
#'   the installed package location (system.file) and then inst/extdata/ in the repo.
#' @return The R object saved in the .rds file (your fitted model).
#' @examples
#' # model <- load_bpvs_model("FE")
#' @export
load_bpvs_model <- function(type = c("FE", "RE"), model_path = NULL) {
  type <- match.arg(type)
  fname <- switch(type, FE = "bpvs_fe_model.rds", RE = "bpvs_re_model.rds")

  if (!is.null(model_path)) {
    if (!file.exists(model_path)) stop("Provided model_path does not exist: ", model_path)
    return(readRDS(model_path))
  }

  # Try installed package location first
  candidate <- system.file("extdata", fname, package = "BPVS")
  if (nzchar(candidate) && file.exists(candidate)) {
    return(readRDS(candidate))
  }

  # Fallback to repository path when running from source
  candidate_local <- file.path("inst", "extdata", fname)
  if (file.exists(candidate_local)) {
    return(readRDS(candidate_local))
  }

  stop("Model file not found. Save your model as inst/extdata/", fname, " or provide model_path.")
}

#' Generic compatibility loader (defaults to FE)
#'
#' Backwards-compatible loader used by older code expecting load_model().
#' @param model_path Optional path; if provided it will be used to load the model file.
#' @param default_type Default model type to load when not specified ("FE" or "RE").
#' @export
load_model <- function(model_path = NULL, default_type = c("FE", "RE")) {
  default_type <- match.arg(default_type)
  load_bpvs_model(type = default_type, model_path = model_path)
}

#' Predict using BPVS model (FE or RE)
#'
#' A wrapper to call the appropriate predict method for your saved BPVS model. By default
#' it returns a list with element `predictions` (the raw prediction output) and `selected`
#' (if your model object stores selected covariates in $selected or a named element).
#'
#' @param newdata A data.frame with the features required by the model.
#' @param type Character, one of "FE" or "RE". Which model to use for prediction.
#' @param model Optional pre-loaded model object. If provided, model_path is ignored.
#' @param model_path Optional path to a .rds file used only when model is NULL.
#' @return A list with at least element `predictions`. If your model records selected covariates
#'   (e.g., model$selected) that will be included as `selected`.
#' @export
predict_bpvs <- function(newdata, type = c("FE", "RE"), model = NULL, model_path = NULL) {
  type <- match.arg(type)
  if (is.null(model)) model <- load_bpvs_model(type = type, model_path = model_path)

  # Default predict call - adapt to your model class if necessary.
  preds <- tryCatch({
    stats::predict(model, newdata)
  }, error = function(e) {
    stop("Prediction failed; adapt predict_bpvs() for your model class. Original error: ", e$message)
  })

  selected <- NULL
  if (!is.null(model$selected)) selected <- model$selected
  if (!is.null(model$selected_covariates)) selected <- model$selected_covariates

  list(predictions = preds, selected = selected)
}

#' Backwards-compatible predict_model wrapper
#'
#' This wrapper keeps the original simple interface while allowing type selection via an argument.
#' @param newdata A data.frame with the features required by the model.
#' @param model Optional model object.
#' @param model_path Optional path to model file.
#' @param type Optional type, one of "FE" or "RE". Defaults to "FE" for compatibility.
#' @export
predict_model <- function(newdata, model = NULL, model_path = NULL, type = c("FE", "RE")) {
  # This snippet is required to be shown exactly as provided in the thread-scoped file:
  if (is.null(model)) model <- load_model(model_path)

  type <- match.arg(type)
  # If model was provided but doesn't correspond to the requested type, we still use it.
  if (!is.null(model)) {
    preds <- tryCatch({ stats::predict(model, newdata) }, error = function(e) stop(e))
    selected <- if (!is.null(model$selected)) model$selected else NULL
    return(list(predictions = preds, selected = selected))
  }

  # Otherwise delegate to predict_bpvs which will load the appropriate model file
  predict_bpvs(newdata = newdata, type = type, model = NULL, model_path = model_path)
}
