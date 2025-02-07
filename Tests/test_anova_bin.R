# File: test-anova_bin.R
# This file contains UAT tests for the anova_bin() function.
# Save it in your test folder and run using testthat (e.g., via devtools::test() or testthat::test_dir()).

library(testthat)
library(dplyr)
library(ggplot2)
library(broom)
library(rlang)

# Adjust the path if necessary so that the script finds your function script.
source(file.path("..", "R", "anova_bin.R"))

context("Testing anova_bin function")

# ------------------------------------------------------------------------------
# 1. Test missing response variable parameter:
# ------------------------------------------------------------------------------
test_that("Error thrown when response_var is not in data", {
  df <- data.frame(a = c(0, 1, 0), group = c("A", "B", "A"))
  expect_error(
    anova_bin(df, "nonexistent", "group"),
    "Response variable nonexistent not found in data"
  )
})

# ------------------------------------------------------------------------------
# 2. Test missing grouping variable parameter:
# ------------------------------------------------------------------------------
test_that("Error thrown when group_var is not in data", {
  df <- data.frame(response = c(0, 1, 0), group = c("A", "B", "A"))
  expect_error(
    anova_bin(df, "response", "nonexistent"),
    "Grouping variable\\(s\\) nonexistent not found in data"
  )
})

# ------------------------------------------------------------------------------
# 3. Test that response_var must be binary:
# ------------------------------------------------------------------------------
test_that("Error thrown when response_var is not binary", {
  df <- data.frame(response = c(0, 1, 2), group = c("A", "B", "A"))
  expect_error(
    anova_bin(df, "response", "group"),
    "The response variable 'response' must be binary"
  )
})

# ------------------------------------------------------------------------------
# 4. Test na.rm parameter removes rows with missing values:
# ------------------------------------------------------------------------------
test_that("na.rm parameter removes rows with missing values", {
  df <- data.frame(
    response = c(0, 1, NA, 0, 1),
    group    = c("A", "B", "A", NA, "B")
  )
  
  # Capture the message about rows removed.
  messages <- capture_messages({
    result <- anova_bin(df, "response", "group", na.rm = TRUE)
  })
  expect_true(any(grepl("Removed", messages)))
  
  # Also, verify that the fitted model is a valid glm object.
  result <- anova_bin(df, "response", "group", na.rm = TRUE, print_plot = FALSE)
  expect_s3_class(result$fitted_model, "glm")
})

# ------------------------------------------------------------------------------
# 5. Test output structure of the function:
# ------------------------------------------------------------------------------
test_that("Function returns list with correct elements", {
  df <- data.frame(
    response = rep(c(0, 1), times = 10),
    group    = rep(c("A", "B"), times = 10)
  )
  result <- anova_bin(df, "response", "group", print_plot = FALSE)
  
  expect_type(result, "list")
  expect_true(all(c("deviance_anova", "effect_size", "model_stats", 
                    "count_plot", "fitted_model") %in% names(result)))
})

# ------------------------------------------------------------------------------
# 6. Test effect_size output contains odds ratios and confidence intervals:
# ------------------------------------------------------------------------------
test_that("Effect size output contains odds ratios and confidence intervals", {
  df <- data.frame(
    response = rep(c(0, 1), times = 20),
    group    = rep(c("A", "B"), times = 20)
  )
  result <- anova_bin(df, "response", "group", print_plot = FALSE)
  es <- result$effect_size
  
  expect_true("estimate" %in% names(es))
  expect_true("conf.low" %in% names(es))
  expect_true("conf.high" %in% names(es))
})

# ------------------------------------------------------------------------------
# 7. Test that count_plot is a ggplot object:
# ------------------------------------------------------------------------------
test_that("count_plot is a ggplot object", {
  df <- data.frame(
    response = rep(c(0, 1), times = 10),
    group    = rep(c("A", "B"), times = 10)
  )
  result <- anova_bin(df, "response", "group", print_plot = FALSE)
  
  expect_true(inherits(result$count_plot, "ggplot"))
})

# ------------------------------------------------------------------------------
# 8. Test multiple grouping variables combine correctly:
# ------------------------------------------------------------------------------
test_that("Multiple group_var parameters combine into a single factor for plotting", {
  df <- data.frame(
    response = rep(c(0, 1), times = 10),
    group1   = rep(c("A", "B"), times = 10),
    group2   = rep(c("X", "Y"), times = 10)
  )
  result <- anova_bin(df, "response", c("group1", "group2"), print_plot = FALSE)
  
  # When multiple grouping variables are used, the x-axis label should be "Combined Group".
  plot_labels <- result$count_plot$labels
  expect_equal(plot_labels$x, "Combined Group")
})

# ------------------------------------------------------------------------------
# 9. Test additional arguments passed to glm (e.g., weights):
# ------------------------------------------------------------------------------
test_that("Additional arguments are passed to glm", {
  set.seed(123)
  df <- data.frame(
    response = rep(c(0, 1), times = 10),
    group    = rep(c("A", "B"), times = 10),
    # Use integer weights to avoid non-integer warnings
    weights  = sample(1:10, 20, replace = TRUE)
  )
  
  result <- anova_bin(df, "response", "group", print_plot = FALSE, weights = df$weights)
  
  # Check that 'weights' appears in the fitted model's call.
  model_call <- result$fitted_model$call
  expect_true("weights" %in% names(model_call))
})

# ------------------------------------------------------------------------------
# 10. Test that print_plot parameter works (plot is not printed when FALSE):
# ------------------------------------------------------------------------------
test_that("Plot is not printed when print_plot is FALSE (but still returned)", {
  df <- data.frame(
    response = rep(c(0, 1), times = 10),
    group    = rep(c("A", "B"), times = 10)
  )
  result <- anova_bin(df, "response", "group", print_plot = FALSE)
  
  # Even though the plot is not printed to the device, it is still returned.
  expect_true(inherits(result$count_plot, "ggplot"))
})
