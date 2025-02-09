## ----------------------------------------
## UAT - Repeated measures ANOVA
## ----------------------------------------

# --- Load Required Packages for Testing ---
suppressWarnings(suppressPackageStartupMessages(library(afex)))
suppressWarnings(suppressPackageStartupMessages(library(emmeans)))
suppressWarnings(suppressPackageStartupMessages(library(car)))
suppressWarnings(suppressPackageStartupMessages(library(effectsize)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
library(testthat)

## ----------------------------
## Create Sample Test Data
## ----------------------------
set.seed(123)
subject_ids <- rep(1:20, each = 3)
Time <- factor(rep(c("Pre", "Post", "FollowUp"), times = 20))
Group <- factor(rep(sample(c("Control", "Treatment"), 20, replace = TRUE), each = 3))
Score <- rnorm(60, mean = 100, sd = 15)
test_data <- data.frame(ID = subject_ids, Score = Score, Time = Time, Group = Group)

# Create a faulty dataset missing the subject column (for negative testing)
faulty_data <- test_data[, -1]

## ----------------------------
## 1. Test: Input Checks and Factor Conversion
## ----------------------------
test_that("anova_rm: Error when required columns are missing", {
  expect_error(
    anova_rm(data = faulty_data,
             subject = "ID",
             dv = "Score",
             within = "Time"),
    "missing"
  )
})

test_that("anova_rm: Factorize parameter works correctly", {
  res1 <- anova_rm(data = test_data,
                   subject = "ID",
                   dv = "Score",
                   within = "Time",
                   between = "Group",
                   factorize = TRUE,
                   diagnostics = FALSE,
                   verbose = FALSE)
  expect_true(is.factor(res1$anova_obj$data[["ID"]]))
  
  res2 <- anova_rm(data = test_data,
                   subject = "ID",
                   dv = "Score",
                   within = "Time",
                   between = "Group",
                   factorize = FALSE,
                   diagnostics = FALSE,
                   verbose = FALSE)
  expect_false(is.factor(res2$anova_obj$data[["ID"]]))
})

## ----------------------------
## 2. Test: Diagnostics Parameter
## ----------------------------
test_that("anova_rm: Diagnostics parameter controls diagnostic plot generation", {
  res_diag <- anova_rm(data = test_data,
                       subject = "ID",
                       dv = "Score",
                       within = "Time",
                       between = "Group",
                       diagnostics = TRUE,
                       verbose = FALSE)
  # Only expect diagnostic_plots if the underlying model is univariate (i.e., not an mlm)
  if (inherits(res_diag$anova_obj$lm, "lm") && !inherits(res_diag$anova_obj$lm, "mlm")) {
    expect_true(!is.null(res_diag$diagnostic_plots))
  } else {
    expect_null(res_diag$diagnostic_plots)
  }
  
  res_no_diag <- anova_rm(data = test_data,
                          subject = "ID",
                          dv = "Score",
                          within = "Time",
                          between = "Group",
                          diagnostics = FALSE,
                          verbose = FALSE)
  expect_null(res_no_diag$diagnostic_plots)
})
