test_that("bpvs_backend returns a valid backend", {
  expect_true(bpvs_backend() %in% c("cmdstanr", "rstan"))
})

test_that("extract_slope_terms pulls RHS terms", {
  f <- y ~ x1 + x2 + x3
  expect_identical(BPVS:::extract_slope_terms(f), c("x1", "x2", "x3"))
})

test_that("null-coalescing operator works", {
  expect_null(NULL %||% NULL)
  expect_identical(1 %||% 2, 1)
  expect_identical(NULL %||% "a", "a")
})

test_that("compute_pip computes plausible PIPs", {
  set.seed(1)
  fit <- posterior::as_draws_df(data.frame(
    b_x1 = rnorm(500, mean = 0.8, sd = 0.1),
    b_x2 = rnorm(500, mean = 0.0, sd = 0.01)
  ))
  attr(fit, "bpvs_effects") <- "re"
  attr(fit, "bpvs_slopes")  <- c("x1", "x2")
  data <- data.frame(Income = rnorm(100, sd = 5))
  pip <- compute_pip(fit, data, outcome_var = "Income", effect_threshold = 0.01)
  expect_equal(nrow(pip), 2)
  expect_gt(pip$PIP[pip$Variable == "x1"], 0.9)
  expect_lt(pip$PIP[pip$Variable == "x2"], 0.1)
})

test_that("rank_variables assigns selection labels", {
  df <- data.frame(
    Variable = c("a", "b", "c"),
    Estimate = c(0.1, 0.2, 0.3),
    CI_lower = c(-0.1, 0.0, 0.2),
    CI_upper = c(0.3, 0.4, 0.4),
    PIP      = c(0.95, 0.3, 0.05),
    stringsAsFactors = FALSE
  )
  r <- rank_variables(df)
  expect_equal(r$Selection, c("Selected", "Weak", "Excluded"))
  expect_equal(r$Rank, 1:3)
  expect_equal(r$Variable, c("a", "b", "c"))
})

test_that("simulate_panel_data produces correct structure", {
  set.seed(456)
  dat <- simulate_panel_data(N_units = 10, T_periods = 5,
                             beta_true = c(1.5, -0.8, 0),
                             sigma = 0.5, tau_unit = 0.4)
  expect_equal(nrow(dat), 50)
  expect_equal(ncol(dat), 3 + 3)
  expect_true(all(c("id", "year", "Income", "x1", "x2", "x3") %in% names(dat)))
  expect_s3_class(dat$id, "factor")
})

test_that("developing dataset is bundled", {
  data(developing, package = "BPVS")
  expect_equal(dim(developing), c(600, 13))
  expect_true(all(c("Growth", "TOP", "FDI") %in% names(developing)))
})
