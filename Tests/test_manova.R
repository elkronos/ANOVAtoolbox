# =========================
# UAT: Full Test Suite Using testthat
# =========================
if (requireNamespace("testthat", quietly = TRUE)) {
  
  message("\n==================== Running UAT Tests ====================")
  library(testthat)
  
  # Test 1: Basic functionality with iris dataset (2 DVs, single grouping variable)
  test_that("manova_analysis works on iris dataset", {
    data(iris)
    res <- manova_analysis(cbind(Sepal.Length, Sepal.Width) ~ Species, data = iris)
    
    expect_s3_class(res$manova_fit, "manova")
    expect_true(!is.null(res$manova_summary))
    expect_true(is.list(res$aov_results))
    expect_true(is.list(res$effect_sizes))
    # If emmeans is installed, ci_results should be a non-empty list; else it should be empty.
    if (requireNamespace("emmeans", quietly = TRUE)) {
      expect_true(length(res$ci_results) > 0)
    }
  })
  
  # Test 2: Provide group_var explicitly when multiple predictors exist.
  test_that("manova_analysis handles explicit group_var when multiple predictors exist", {
    # Create synthetic data with two predictors and two DVs.
    set.seed(123)
    dat <- data.frame(
      DV1 = rnorm(100),
      DV2 = rnorm(100),
      Group = sample(letters[1:3], 100, replace = TRUE),
      Covar = rnorm(100)
    )
    # Formula has two predictors; explicitly provide group_var.
    res <- manova_analysis(cbind(DV1, DV2) ~ Group + Covar, data = dat, group_var = "Group")
    expect_s3_class(res$manova_fit, "manova")
  })
  
  # Test 3: When group_var is missing from data, function should error.
  test_that("manova_analysis errors when group_var is not in data", {
    dat <- iris
    expect_error(manova_analysis(cbind(Sepal.Length, Sepal.Width) ~ Species, data = dat, group_var = "Nonexistent"))
  })
  
  # Test 4: Check that diagnostic plots are generated (if plotDiagnostics = TRUE)
  test_that("Diagnostic plots are generated when requested", {
    data(iris)
    res <- manova_analysis(cbind(Sepal.Length, Sepal.Width) ~ Species, data = iris, plotDiagnostics = TRUE)
    expect_true(is.list(res$diagnostic_plots))
    # For each DV, expect a list with two plots.
    for (dv in names(res$diagnostic_plots)) {
      expect_true(all(c("residual_plot", "qq_plot") %in% names(res$diagnostic_plots[[dv]])))
    }
  })
  
  # Test 5: Check that when plotDiagnostics is FALSE, diagnostic_plots is empty.
  test_that("Diagnostic plots are skipped when plotDiagnostics is FALSE", {
    data(iris)
    res <- manova_analysis(cbind(Sepal.Length, Sepal.Width) ~ Species, data = iris, plotDiagnostics = FALSE)
    expect_equal(length(res$diagnostic_plots), 0)
  })
  
  message("All UAT tests passed.")
  
} else {
  message("Package 'testthat' not installed; skipping UAT tests.")
}
