# test_anova_kw.R

# Load required libraries for testing
library(testthat)
library(data.table)
library(ggplot2)
library(nortest)
library(dunn.test)

# Source the anova_kw function
# Adjust the path to where your script is located.
source(file.path("..", "R", "anova_kw.R"))

#####
# Test 1: Data Validation and Parameter Checks
#####

test_that("Error when response_var is not found in data", {
  df <- data.frame(x = 1:10, group = rep(c("A", "B"), each = 5))
  expect_error(anova_kw(df, "y", "group"), "response_var not found in data")
})

test_that("Error when response_var is non-numeric", {
  df <- data.frame(a = letters[1:10], group = rep(c("A", "B"), each = 5))
  expect_error(anova_kw(df, "a", "group"), "response_var must be numeric")
})

test_that("Error when a grouping variable is missing from data", {
  df <- data.frame(a = 1:10, group = rep(c("A", "B"), each = 5))
  expect_error(anova_kw(df, "a", "nonexistent"),
               "The following group_vars_vec elements are not found in data")
})

test_that("Grouping variables are converted to factors if not already", {
  # Provide a data frame where group is character (not factor)
  df <- data.frame(a = rnorm(30), group = rep(c("A", "B", "C"), each = 10),
                   stringsAsFactors = FALSE)
  # Should not error and should produce a ggplot boxplot
  result <- anova_kw(df, "a", "group", plot_qq = FALSE)
  expect_true(inherits(result$boxplot, "ggplot"))
})

test_that("Error when there is only one unique group", {
  df <- data.frame(a = rnorm(10), group = rep("A", 10))
  expect_error(anova_kw(df, "a", "group"),
               "Not enough groups for analysis. At least 2 unique groups are required.")
})

#####
# Test 2: Function Output for a Single Grouping Variable
#####

test_that("Valid output for a single grouping variable", {
  set.seed(123)
  df <- data.frame(
    group = rep(c("A", "B", "C"), each = 30),
    a = c(rnorm(30, mean = 5), rnorm(30, mean = 7), rnorm(30, mean = 4))
  )
  # Ensure group is a factor
  df$group <- as.factor(df$group)
  
  result <- anova_kw(df, "a", "group", plot_qq = TRUE, posthoc_method = "bh")
  
  # Check that result is a list with the required components
  expect_true(is.list(result))
  
  # Check assumptions list
  expect_true(is.list(result$assumptions))
  expect_true(is.numeric(result$assumptions$residuals))
  expect_true(inherits(result$assumptions$ad_test, "htest"))
  expect_true(!is.null(result$assumptions$qqplot))
  
  # Check the Kruskal-Wallis test result
  expect_true(inherits(result$kruskal_test, "htest"))
  
  # Check the post-hoc (Dunn's) test result
  expect_true(!is.null(result$posthoc))
  
  # Check that the boxplot is a ggplot object
  expect_true(inherits(result$boxplot, "ggplot"))
})

#####
# Test 3: Function Output for Multiple Grouping Variables
#####

test_that("Valid output for multiple grouping variables", {
  set.seed(123)
  df <- data.frame(
    group1 = rep(c("A", "B"), each = 30),
    group2 = rep(c("X", "Y"), times = 30),
    a = c(rnorm(30, mean = 5), rnorm(30, mean = 7))
  )
  # Ensure both grouping variables are factors
  df$group1 <- as.factor(df$group1)
  df$group2 <- as.factor(df$group2)
  
  result <- anova_kw(df, "a", c("group1", "group2"),
                     plot_qq = FALSE, posthoc_method = "bonferroni")
  
  # Basic structure checks
  expect_true(is.list(result))
  expect_true(is.numeric(result$assumptions$residuals))
  expect_true(inherits(result$assumptions$ad_test, "htest"))
  
  # QQ plot should be NULL because plot_qq = FALSE
  expect_null(result$assumptions$qqplot)
  
  expect_true(inherits(result$kruskal_test, "htest"))
  expect_true(!is.null(result$posthoc))
  expect_true(inherits(result$boxplot, "ggplot"))
})

#####
# Test 4: Testing the plot_qq Parameter Behavior
#####

test_that("plot_qq parameter toggles the QQ plot output", {
  set.seed(123)
  df <- data.frame(
    group = rep(c("A", "B", "C"), each = 30),
    a = c(rnorm(30, mean = 5), rnorm(30, mean = 7), rnorm(30, mean = 4))
  )
  df$group <- as.factor(df$group)
  
  # With plot_qq = TRUE, the qqplot should be generated (not NULL)
  result1 <- anova_kw(df, "a", "group", plot_qq = TRUE)
  expect_true(!is.null(result1$assumptions$qqplot))
  
  # With plot_qq = FALSE, the qqplot should be NULL
  result2 <- anova_kw(df, "a", "group", plot_qq = FALSE)
  expect_null(result2$assumptions$qqplot)
})

#####
# Test 5: Testing the posthoc_method Parameter
#####

test_that("Valid posthoc_method produces a Dunn test result", {
  set.seed(123)
  df <- data.frame(
    group = rep(c("A", "B", "C"), each = 30),
    a = c(rnorm(30, mean = 5), rnorm(30, mean = 7), rnorm(30, mean = 4))
  )
  df$group <- as.factor(df$group)
  
  # Use "bonferroni" as an alternative adjustment method
  result <- anova_kw(df, "a", "group", posthoc_method = "bonferroni")
  expect_true(!is.null(result$posthoc))
})

test_that("Invalid posthoc_method produces an error", {
  set.seed(123)
  df <- data.frame(
    group = rep(c("A", "B", "C"), each = 30),
    a = c(rnorm(30, mean = 5), rnorm(30, mean = 7), rnorm(30, mean = 4))
  )
  df$group <- as.factor(df$group)
  
  # Passing an invalid method should trigger an error from dunn.test
  expect_error(anova_kw(df, "a", "group", posthoc_method = "invalid_method"))
})
