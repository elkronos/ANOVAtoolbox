# test_anova_glm.R
#
# UAT for the anova_glm() function.
# This test file uses the testthat framework to validate the behavior of each parameter.
# Some warnings from data.table and ggplot2 are suppressed in the tests.

library(testthat)
library(data.table)
library(ggplot2)
library(multcomp)
library(car)
library(broom)

# Source the script containing anova_glm()
# (Change the path below if needed)
source("C:/repos/ANOVAtoolbox/R/anova_glm.R")

context("UAT for anova_glm function")

# ------------------------------------------------------------------------------
# Prepare sample datasets for testing
# ------------------------------------------------------------------------------

# Data for a single grouping variable (3 groups)
set.seed(123)
data_single <- data.frame(
  group = rep(c("Group A", "Group B", "Group C"), each = 50),
  value = c(rnorm(50, mean = 5, sd = 1),
            rnorm(50, mean = 7, sd = 1.5),
            rnorm(50, mean = 4, sd = 0.8))
)

# Data for two grouping variables
data_double <- data.frame(
  group1 = rep(c("Group A", "Group B"), each = 50),
  group2 = rep(c("Group X", "Group Y"), times = 50),
  value = c(rnorm(50, mean = 5, sd = 1),
            rnorm(50, mean = 7, sd = 1.5),
            rnorm(50, mean = 4, sd = 0.8),
            rnorm(50, mean = 6, sd = 1.2))
)

# Data for binary outcome (to test non-gaussian family)
data_binomial <- data.frame(
  group = rep(c("Group A", "Group B"), 100),
  outcome = rbinom(200, size = 1, prob = 0.5)
)

# ------------------------------------------------------------------------------
# Test 1: Error when response_var is not found in the data
# ------------------------------------------------------------------------------
test_that("Error when response_var is not in the data", {
  expect_error(
    anova_glm(data_single, "nonexistent", c("group")),
    regexp = "The response variable"
  )
})

# ------------------------------------------------------------------------------
# Test 2: Error when one or more grouping variables are not found
# ------------------------------------------------------------------------------
test_that("Error when a grouping variable is missing", {
  expect_error(
    anova_glm(data_single, "value", c("nonexistent")),
    regexp = "grouping variable"
  )
})

# ------------------------------------------------------------------------------
# Test 3: Grouping variable is converted to factor if not already
# ------------------------------------------------------------------------------
test_that("Grouping variable conversion to factor", {
  # Create a dataset with grouping variable as character
  data_char <- data_single
  data_char$group <- as.character(data_char$group)
  
  # Suppress warnings that come from data.table copying and aes_string deprecation.
  expect_message(
    res <- suppressWarnings(anova_glm(data_char, "value", c("group"))),
    regexp = "Variable 'group' was converted to a factor"
  )
  
  # Retrieve the x-axis labels using the modern ggplot2 API.
  x_labels <- suppressWarnings(
    ggplot_build(res$boxplot)$layout$panel_scales_x[[1]]$get_labels()
  )
  
  # Check that we successfully retrieved a non-null set of labels.
  expect_false(is.null(x_labels))
  
  # Compare the labels with the expected factor levels.
  expected_labels <- levels(as.factor(data_char$group))
  expect_equal(sort(as.character(x_labels)), sort(expected_labels))
})

# ------------------------------------------------------------------------------
# Test 4: Valid output with a single grouping variable and plot_residuals = TRUE
# ------------------------------------------------------------------------------
test_that("Output structure with single grouping variable (plot_residuals = TRUE)", {
  res <- suppressWarnings(anova_glm(data_single, "value", c("group"), plot_residuals = TRUE))
  
  expect_type(res, "list")
  expect_true("model" %in% names(res))
  expect_true("model_stats" %in% names(res))
  expect_true("anova_test" %in% names(res))
  expect_true("posthoc" %in% names(res))
  expect_true("effect_size" %in% names(res))
  expect_true("confidence_intervals" %in% names(res))
  expect_true("residuals" %in% names(res))
  expect_true("residuals_plot" %in% names(res))
  expect_true("boxplot" %in% names(res))
  
  # Check that residuals_plot is a ggplot object.
  expect_true(inherits(res$residuals_plot, "ggplot"))
  # Check that boxplot is a ggplot object.
  expect_true(inherits(res$boxplot, "ggplot"))
})

# ------------------------------------------------------------------------------
# Test 5: Valid output with multiple grouping variables and sum_squares_type = "III"
# ------------------------------------------------------------------------------
test_that("Output structure with multiple grouping variables", {
  res <- suppressWarnings(anova_glm(data_double, "value", c("group1", "group2"),
                                    sum_squares_type = "III", plot_residuals = FALSE))
  
  expect_type(res, "list")
  expect_true("model" %in% names(res))
  expect_true("model_stats" %in% names(res))
  expect_true("anova_test" %in% names(res))
  expect_true("posthoc" %in% names(res))
  expect_true("effect_size" %in% names(res))
  expect_true("confidence_intervals" %in% names(res))
  expect_true("residuals" %in% names(res))
  expect_true("residuals_plot" %in% names(res))
  expect_true("boxplot" %in% names(res))
  
  # plot_residuals = FALSE should yield a NULL residuals_plot.
  expect_null(res$residuals_plot)
  
  # The boxplot should be a ggplot object and include facetting.
  expect_true(inherits(res$boxplot, "ggplot"))
  expect_true(!is.null(res$boxplot$facet))
})

# ------------------------------------------------------------------------------
# Test 6: Error for unsupported sum_squares_type parameter
# ------------------------------------------------------------------------------
test_that("Error for unsupported sum_squares_type", {
  expect_error(
    suppressWarnings(anova_glm(data_single, "value", c("group"), sum_squares_type = "IV")),
    regexp = "Unsupported sum_squares_type"
  )
})

# ------------------------------------------------------------------------------
# Test 7: Post-hoc analysis output when there are >2 levels
# ------------------------------------------------------------------------------
test_that("Post-hoc analysis for >2 groups", {
  res <- suppressWarnings(anova_glm(data_single, "value", c("group")))
  
  # data_single has 3 groups so posthoc should be computed (i.e. not NULL).
  expect_false(is.null(res$posthoc))
  # Check that the posthoc object contains an element "test" (from summary.glht).
  expect_true("test" %in% names(res$posthoc))
})

# ------------------------------------------------------------------------------
# Test 8: Post-hoc analysis is NULL when there are only 2 groups
# ------------------------------------------------------------------------------
test_that("No post-hoc analysis for only 2 groups", {
  data_two <- data.frame(
    group = rep(c("Group A", "Group B"), each = 50),
    value = c(rnorm(50, mean = 5, sd = 1),
              rnorm(50, mean = 7, sd = 1.5))
  )
  res <- suppressWarnings(anova_glm(data_two, "value", c("group")))
  
  # With only two groups, posthoc output should be NULL.
  expect_null(res$posthoc)
})

# ------------------------------------------------------------------------------
# Test 9: Confidence intervals output is a data frame with proper dimensions
# ------------------------------------------------------------------------------
test_that("Confidence intervals output", {
  res <- suppressWarnings(anova_glm(data_single, "value", c("group")))
  
  expect_true(is.data.frame(res$confidence_intervals))
  # The number of rows should equal the number of model coefficients (including intercept).
  expect_equal(nrow(res$confidence_intervals), length(coef(res$model)))
})

# ------------------------------------------------------------------------------
# Test 10: Effect size is computed correctly
# ------------------------------------------------------------------------------
test_that("Effect size computation", {
  res <- suppressWarnings(anova_glm(data_single, "value", c("group")))
  
  # Compute pseudo R-squared manually.
  pseudo_r2_manual <- 1 - (res$model$deviance / res$model$null.deviance)
  expect_equal(res$effect_size, pseudo_r2_manual)
})

# ------------------------------------------------------------------------------
# Test 11: GLM family parameter works (using binomial family)
# ------------------------------------------------------------------------------
test_that("GLM family parameter works with binomial family", {
  res <- suppressWarnings(anova_glm(data_binomial, "outcome", c("group"),
                                    family = binomial(), plot_residuals = TRUE))
  
  expect_true(inherits(res$model, "glm"))
  expect_equal(res$model$family$family, "binomial")
  expect_true(inherits(res$residuals_plot, "ggplot"))
})
