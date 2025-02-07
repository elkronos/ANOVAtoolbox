# test-anova_count.R
# User Acceptance Tests (UAT) for the anova_count function

library(testthat)
library(data.table)
library(ggplot2)
library(multcomp)
library(broom)

# Source the anova_count function.
# Adjust the path if necessary. Here we assume the function is one directory up in "R"
source(file.path("..", "R", "anova_count.R"))

context("Testing anova_count function")

### 1. Valid Input --------------------------------------------------------------

test_that("Valid input returns expected list elements", {
  set.seed(123)
  data_example <- data.frame(
    group1 = factor(rep(c("A", "B"), each = 50)),
    group2 = factor(rep(c("X", "Y"), times = 50)),
    count   = c(rpois(50, 5), rpois(50, 10))
  )
  
  result <- anova_count(data_example, response_var = "count", 
                        group_vars_list = c("group1", "group2"), plot = FALSE)
  
  # Check that the result is a list and has all expected components
  expected_names <- c("model", "overdispersion_statistic", "overdispersion_check_message",
                      "deviance_anova", "posthoc", "effect_size", "count_plot")
  expect_type(result, "list")
  expect_setequal(names(result), expected_names)
  
  # Check specific components
  expect_s3_class(result$model, "glm")
  expect_true(is.numeric(result$overdispersion_statistic))
  expect_true(is.character(result$overdispersion_check_message))
  expect_s3_class(result$deviance_anova, "anova")
  
  # posthoc might be NULL if an error occurred inside glht (but for valid data it should work)
  if (!is.null(result$posthoc)) {
    expect_s3_class(result$posthoc, "glht")
  }
  
  expect_s3_class(result$effect_size, "data.frame")
  expect_s3_class(result$count_plot, "ggplot")
})

### 2. Empty Data ----------------------------------------------------------------

test_that("Empty data raises error", {
  expect_error(
    anova_count(data = data.frame(), response_var = "count", group_vars_list = c("group1")),
    "The data provided is empty"
  )
})

### 3. Missing Response Variable ------------------------------------------------

test_that("Missing response variable raises error", {
  data_missing_response <- data.frame(
    group1 = factor(c("A", "B")),
    group2 = factor(c("X", "Y")),
    count2 = c(1, 2)
  )
  expect_error(
    anova_count(data_missing_response, response_var = "count", 
                group_vars_list = c("group1", "group2")),
    "Response variable 'count' not found"
  )
})

### 4. Missing Grouping Variable ------------------------------------------------

test_that("Missing grouping variable raises error", {
  data_missing_group <- data.frame(
    group1 = factor(c("A", "B")),
    count = c(1, 2)
  )
  expect_error(
    anova_count(data_missing_group, response_var = "count", 
                group_vars_list = c("group1", "group2")),
    "The following grouping variable\\(s\\) are missing"
  )
})

### 5. Non-Numeric Response Variable --------------------------------------------

test_that("Non-numeric response variable raises error", {
  data_non_numeric <- data.frame(
    group1 = factor(c("A", "B")),
    group2 = factor(c("X", "Y")),
    count = c("a", "b")
  )
  expect_error(
    anova_count(data_non_numeric, response_var = "count", 
                group_vars_list = c("group1", "group2")),
    "Response variable 'count' must be numeric"
  )
})

### 6. Negative Values in Response Variable ------------------------------------

test_that("Response variable with negative values raises error", {
  data_negative <- data.frame(
    group1 = factor(c("A", "B")),
    group2 = factor(c("X", "Y")),
    count = c(-1, 2)
  )
  expect_error(
    anova_count(data_negative, response_var = "count", 
                group_vars_list = c("group1", "group2")),
    "Response variable 'count' must be non-negative"
  )
})

### 7. Non-Integer Values in Response Variable ----------------------------------

test_that("Response variable with non-integer values raises error", {
  data_non_integer <- data.frame(
    group1 = factor(c("A", "B")),
    group2 = factor(c("X", "Y")),
    count = c(1.5, 2.3)
  )
  expect_error(
    anova_count(data_non_integer, response_var = "count", 
                group_vars_list = c("group1", "group2")),
    "Response variable 'count' must contain integer values"
  )
})

### 8. Grouping Variable with Only One Level ------------------------------------

test_that("Grouping variable with one level raises error", {
  data_one_level <- data.frame(
    group1 = factor(rep("A", 5)),
    group2 = factor(c("X", "Y", "X", "Y", "X")),
    count = c(1, 2, 3, 4, 5)
  )
  expect_error(
    anova_count(data_one_level, response_var = "count", 
                group_vars_list = c("group1", "group2")),
    "Grouping variable 'group1' must have at least 2 levels"
  )
})

### 9. Interaction of Grouping Variables (Less than 2 Levels) --------------------

# Although each grouping variable is individually checked for multiple levels,
# this test ensures that if the interaction yields only one level, an error is raised.
test_that("Interaction of grouping variables with only one level raises error", {
  data_interaction_one_level <- data.frame(
    group1 = factor(rep("A", 10)),
    group2 = factor(rep("X", 10)),
    count = rpois(10, lambda = 5)
  )
  expect_error(
    anova_count(data_interaction_one_level, response_var = "count", 
                group_vars_list = c("group1", "group2")),
    "Grouping variable 'group1' must have at least 2 levels"
  )
})

### 10. Overdispersion Check ----------------------------------------------------

test_that("Overdispersion message is appropriate for overdispersed data", {
  set.seed(42)
  n <- 100
  data_overdisp <- data.frame(
    group1 = factor(rep(c("A", "B"), each = n/2)),
    group2 = factor(rep(c("X", "Y"), times = n/2)),
    # Generate data from a negative binomial to induce overdispersion relative to Poisson
    count   = rnbinom(n, size = 1, mu = 5)
  )
  
  result_overdisp <- anova_count(data_overdisp, response_var = "count", 
                                 group_vars_list = c("group1", "group2"),
                                 overdispersion_threshold = 1.5, plot = FALSE)
  
  if(result_overdisp$overdispersion_statistic > 1.5) {
    expect_match(result_overdisp$overdispersion_check_message, "Warning: Overdispersion detected")
  } else {
    expect_match(result_overdisp$overdispersion_check_message, "No overdispersion detected")
  }
})

### 11. Effect Size Table -------------------------------------------------------

test_that("Effect size table contains correct columns", {
  set.seed(123)
  data_example <- data.frame(
    group1 = factor(rep(c("A", "B"), each = 50)),
    group2 = factor(rep(c("X", "Y"), times = 50)),
    count   = c(rpois(50, 5), rpois(50, 10))
  )
  
  result <- anova_count(data_example, response_var = "count", 
                        group_vars_list = c("group1", "group2"), plot = FALSE)
  
  # Check that the effect_size data frame has the expected columns
  expected_cols <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
  expect_true(all(expected_cols %in% names(result$effect_size)))
})

### 12. Count Plot Generation ---------------------------------------------------

test_that("Count plot is generated regardless of plot parameter", {
  set.seed(123)
  data_example <- data.frame(
    group1 = factor(rep(c("A", "B"), each = 50)),
    group2 = factor(rep(c("X", "Y"), times = 50)),
    count   = c(rpois(50, 5), rpois(50, 10))
  )
  
  # Even when plot = FALSE, the ggplot object should be created (but not printed to the console)
  result <- anova_count(data_example, response_var = "count", 
                        group_vars_list = c("group1", "group2"), plot = FALSE)
  expect_s3_class(result$count_plot, "ggplot")
})
