# File: tests/test-anova_welch.R

# It is common to load libraries at the top level.
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(nortest))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(effsize))

# Source the script containing anova_welch.
# Adjust the path below if needed (e.g., "../R/anova_welch.R")
source(file.path("..", "R", "anova_welch.R"))

# ------------------------------------------------------------------------------
# Test 1: Validate that an invalid 'data' argument throws an error
# ------------------------------------------------------------------------------
test_that("Error when 'data' is not a data frame", {
  expect_error(
    anova_welch("not_a_data_frame", "value", c("group")),
    "The 'data' argument must be a data frame."
  )
})

# ------------------------------------------------------------------------------
# Test 2: Validate that a non-existent response variable throws an error
# ------------------------------------------------------------------------------
test_that("Error when response variable is not found in data", {
  df <- data.frame(a = 1:10, b = rnorm(10))
  expect_error(
    anova_welch(df, "nonexistent", c("a")),
    "Response variable nonexistent not found in data."
  )
})

# ------------------------------------------------------------------------------
# Test 3: Validate that a missing grouping variable throws an error
# ------------------------------------------------------------------------------
test_that("Error when a grouping variable is not found in data", {
  df <- data.frame(a = 1:10, value = rnorm(10))
  expect_error(
    anova_welch(df, "value", c("nonexistent")),
    "The following grouping variables were not found in data: nonexistent"
  )
})

# ------------------------------------------------------------------------------
# Test 4: Validate that a non-numeric response variable throws an error
# ------------------------------------------------------------------------------
test_that("Error when the response variable is not numeric", {
  df <- data.frame(value = letters[1:10], group = rep("A", 10))
  expect_error(
    anova_welch(df, "value", c("group")),
    "The response variable must be numeric."
  )
})

# ------------------------------------------------------------------------------
# Test 5: Validate that rows with missing values are removed properly
# ------------------------------------------------------------------------------
test_that("Missing values in response or grouping variables are removed", {
  df <- data.frame(
    value = c(1, 2, NA, 4, 5),
    group = c("A", "B", "A", NA, "B"),
    stringsAsFactors = FALSE
  )
  # Although the cleaned data will be small (and the AD test skipped with a warning),
  # we still check that the summary_stats reflect the removed rows.
  result <- suppressWarnings(anova_welch(df, "value", c("group")))
  
  # In this case, rows 1, 2, and 5 are complete.
  # We expect two groups: one group (A) with 1 observation and one group (B) with 2 observations.
  expect_equal(nrow(result$summary_stats), 2)
  expect_true(all(result$summary_stats$n %in% c(1, 2)))
})

# ------------------------------------------------------------------------------
# Test 6: Validate that the function returns all expected elements for a valid input
# ------------------------------------------------------------------------------
test_that("Valid input returns a result list with expected keys", {
  set.seed(1)
  group <- rep(c("Group A", "Group B", "Group C"), each = 50)
  value <- c(rnorm(50, mean = 5, sd = 1),
             rnorm(50, mean = 7, sd = 1.5),
             rnorm(50, mean = 4, sd = 0.8))
  df <- data.frame(group = factor(group), value = value)
  
  result <- anova_welch(df, "value", c("group"))
  
  # Check for expected list elements
  expect_true(is.list(result))
  expected_keys <- c("assumptions", "welch_anova", "posthoc", 
                     "summary_stats", "effect_sizes", "means_plot")
  expect_true(all(expected_keys %in% names(result)))
  
  # Since there are more than 2 groups, posthoc should not be NULL.
  expect_false(is.null(result$posthoc))
})

# ------------------------------------------------------------------------------
# Test 7: Validate that post-hoc tests are not performed when there are exactly 2 groups
# ------------------------------------------------------------------------------
test_that("Post-hoc tests are not performed for two-group comparisons", {
  set.seed(1)
  group <- rep(c("Group A", "Group B"), each = 30)
  value <- c(rnorm(30, mean = 5, sd = 1),
             rnorm(30, mean = 7, sd = 1.5))
  df <- data.frame(group = factor(group), value = value)
  
  result <- anova_welch(df, "value", c("group"))
  
  # With two groups, the posthoc element should be NULL.
  expect_null(result$posthoc)
})

# ------------------------------------------------------------------------------
# Test 8: Validate that effect sizes are computed for all unique pairwise comparisons
# ------------------------------------------------------------------------------
test_that("Effect sizes computed for all pairwise group comparisons", {
  set.seed(1)
  group <- rep(c("Group A", "Group B", "Group C"), each = 40)
  value <- c(rnorm(40, mean = 5, sd = 1),
             rnorm(40, mean = 7, sd = 1.5),
             rnorm(40, mean = 4, sd = 0.8))
  df <- data.frame(group = factor(group), value = value)
  
  result <- anova_welch(df, "value", c("group"))
  
  # With 3 groups, there should be choose(3,2) = 3 pairwise comparisons.
  expect_equal(length(result$effect_sizes), 3)
  
  # Check that the expected comparison names appear.
  expected_comparisons <- c("Group A vs Group B", 
                            "Group A vs Group C", 
                            "Group B vs Group C")
  expect_true(all(expected_comparisons %in% names(result$effect_sizes)))
})

# ------------------------------------------------------------------------------
# Test 9: Validate that the means plot is a ggplot object
# ------------------------------------------------------------------------------
test_that("Means plot is a valid ggplot object", {
  set.seed(1)
  group <- rep(c("Group A", "Group B", "Group C"), each = 50)
  value <- c(rnorm(50, mean = 5, sd = 1),
             rnorm(50, mean = 7, sd = 1.5),
             rnorm(50, mean = 4, sd = 0.8))
  df <- data.frame(group = factor(group), value = value)
  
  result <- anova_welch(df, "value", c("group"))
  
  expect_true(inherits(result$means_plot, "ggplot"))
})

# ------------------------------------------------------------------------------
# Test 10: Validate that the assumptions output contains a QQ plot and an AD test result
# ------------------------------------------------------------------------------
test_that("Assumptions output is valid", {
  set.seed(1)
  group <- rep(c("Group A", "Group B"), each = 30)
  value <- c(rnorm(30, mean = 5, sd = 1),
             rnorm(30, mean = 7, sd = 1.5))
  df <- data.frame(group = factor(group), value = value)
  
  result <- anova_welch(df, "value", c("group"))
  
  # Check that the QQ plot is a ggplot object.
  expect_true(inherits(result$assumptions$qq_plot, "ggplot"))
  
  # With this dataset the sample size is large enough, so AD test result should be of class "htest"
  expect_true("htest" %in% class(result$assumptions$ad_test))
})
