###############################################################################
# Testing Infrastructure - anova_glm UAT (Generalized Linear Model "ANOVA")
###############################################################################

# Package management for tests -------------------------------------------------
required_packages <- c(
  "car","emmeans","multcomp","ggplot2","dplyr","broom","withr","rlang",
  "scales","magrittr","data.table","MASS"
)
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
suppressPackageStartupMessages({
  library(car)
  library(emmeans)
  library(multcomp)
  library(ggplot2)
  library(dplyr)
  library(broom)
  library(withr)
  library(rlang)
  library(scales)
  library(magrittr)
  library(data.table)
  library(MASS)  # for NB-like generation helpers if needed
})

# NOTE: This UAT assumes `anova_glm()` is already defined in the session.

# Global container for test results -------------------------------------------
test_results <- data.frame(Test = character(), Result = character(), stringsAsFactors = FALSE)

#' Helper: Print and Store Test Result
#'
#' @param test_name Character. Name of the test.
#' @param passed Logical. Whether the test passed.
#' @param note Optional character. Additional notes.
print_and_store_result <- function(test_name, passed, note = NULL) {
  result <- if (isTRUE(passed)) "PASS" else "FAIL"
  cat(sprintf("%-90s [%s]\n", test_name, result))
  if (!is.null(note)) cat("  Note: ", note, "\n", sep = "")
  test_results <<- rbind(
    test_results,
    data.frame(Test = test_name, Result = result, stringsAsFactors = FALSE)
  )
}

###############################################################################
# Synthetic Data Generators (GLM families)  -- FIXED to handle any n
###############################################################################

# Gaussian data: supports 1–2 factors + optional interaction on mean
gen_gaussian_data <- function(n = 400, include_g2 = TRUE, sigma = 1.0) {
  if (include_g2) {
    lv1 <- c("A","B"); lv2 <- c("X","Y")
    cells <- expand.grid(g1 = lv1, g2 = lv2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    idx <- rep(seq_len(nrow(cells)), length.out = n)
    g1 <- factor(cells$g1[idx], levels = lv1)
    g2 <- factor(cells$g2[idx], levels = lv2)
    
    mu <- ifelse(g1=="A" & g2=="X", 0.0,
                 ifelse(g1=="A" & g2=="Y", 0.5,
                        ifelse(g1=="B" & g2=="X", 1.0, 1.5)))
    y <- rnorm(n, mean = mu, sd = sigma)
    data.frame(y = y, g1 = g1, g2 = g2)
  } else {
    g1 <- factor(rep(c("A","B","C"), length.out = n))
    mu <- ifelse(g1=="A", 0.0, ifelse(g1=="B", 0.75, 1.5))
    y <- rnorm(n, mean = mu, sd = sigma)
    data.frame(y = y, g1 = g1)
  }
}

# Binomial (binary) data: 1–2 factors; optional interaction on log-odds
gen_binomial_data <- function(n = 500, include_g2 = TRUE) {
  invlogit <- function(x) 1/(1+exp(-x))
  if (include_g2) {
    lv1 <- c("A","B"); lv2 <- c("X","Y")
    cells <- expand.grid(g1 = lv1, g2 = lv2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    idx <- rep(seq_len(nrow(cells)), length.out = n)
    g1 <- factor(cells$g1[idx], levels = lv1)
    g2 <- factor(cells$g2[idx], levels = lv2)
    
    eta <- ifelse(g1=="A" & g2=="X", -1.0,
                  ifelse(g1=="A" & g2=="Y", -0.3,
                         ifelse(g1=="B" & g2=="X",  0.3,  1.0)))
    p <- invlogit(eta)
    y <- rbinom(n, size = 1, prob = p)
    data.frame(y = y, g1 = g1, g2 = g2)
  } else {
    g1 <- factor(rep(c("A","B","C"), length.out = n))
    eta <- ifelse(g1=="A", -1.2, ifelse(g1=="B", -0.1, 0.9))
    p <- invlogit(eta)
    y <- rbinom(n, size = 1, prob = p)
    data.frame(y = y, g1 = g1)
  }
}

# Poisson data: 1–2 factors; optional interaction effect on log-rate
gen_poisson_data <- function(n = 400, include_g2 = TRUE,
                             rates = matrix(c(5,10,15,20), nrow = 2, byrow = TRUE)) {
  if (include_g2) {
    lv1 <- c("A","B"); lv2 <- c("X","Y")
    cells <- expand.grid(g1 = lv1, g2 = lv2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    idx <- rep(seq_len(nrow(cells)), length.out = n)
    g1 <- factor(cells$g1[idx], levels = lv1)
    g2 <- factor(cells$g2[idx], levels = lv2)
    
    # rates matrix mapping: rows = g1 (A,B), cols = g2 (X,Y)
    lam <- ifelse(g1=="A" & g2=="X", rates[1,1],
                  ifelse(g1=="A" & g2=="Y", rates[1,2],
                         ifelse(g1=="B" & g2=="X", rates[2,1], rates[2,2])))
    y <- rpois(n, lam)
    data.frame(y = y, g1 = g1, g2 = g2)
  } else {
    g1 <- factor(rep(c("A","B","C"), length.out = n))
    lam <- ifelse(g1=="A", 6, ifelse(g1=="B", 12, 18))
    y <- rpois(n, lam)
    data.frame(y = y, g1 = g1)
  }
}

# Overdispersed count data (NegBin generator; often fit as Poisson to trigger note)
gen_overdispersed_counts <- function(n = 500, theta = 0.6) {
  lv1 <- c("A","B"); lv2 <- c("X","Y")
  cells <- expand.grid(g1 = lv1, g2 = lv2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  idx <- rep(seq_len(nrow(cells)), length.out = n)
  g1 <- factor(cells$g1[idx], levels = lv1)
  g2 <- factor(cells$g2[idx], levels = lv2)
  
  mu <- ifelse(g1=="A" & g2=="X", 6,
               ifelse(g1=="A" & g2=="Y", 12,
                      ifelse(g1=="B" & g2=="X", 18, 24)))
  y <- rnbinom(n, size = theta, mu = mu)  # strong overdispersion
  data.frame(y = y, g1 = g1, g2 = g2)
}

# Unbalanced design helper (for Type II vs III checks)
gen_unbalanced_gaussian <- function(nA = 120, nB = 80) {
  g1 <- factor(c(rep("A", nA), rep("B", nB)))
  g2 <- factor(rep(c("X","Y"), length.out = nA + nB))
  set.seed(1)
  mu <- ifelse(g1=="A" & g2=="X", 0.0,
               ifelse(g1=="A" & g2=="Y", 0.6,
                      ifelse(g1=="B" & g2=="X", 1.2, 1.8)))
  y <- rnorm(length(g1), mean = mu, sd = 1.0)
  data.frame(y = y, g1 = g1, g2 = g2)
}

###############################################################################
# Unit Tests
###############################################################################

# Test: input validation & factor coercion ------------------------------------
test_input_validation_glm <- function() {
  df <- gen_gaussian_data(200, include_g2 = FALSE)
  names(df)[names(df) == "y"] <- "resp"
  e1 <- inherits(tryCatch(anova_glm(df, "missing_resp", "g1"), error = identity), "error")
  e2 <- inherits(tryCatch(anova_glm(df, "resp", "missing_g"), error = identity), "error")
  
  # Non-factor coerced to factor
  df2 <- df; df2$g1 <- as.character(df2$g1)
  res <- anova_glm(df2, "resp", "g1")
  coerced <- is.factor(as.data.frame(res$model$model)[["g1"]])
  
  print_and_store_result("input_validation: errors when response missing", e1)
  print_and_store_result("input_validation: errors when group var missing", e2)
  print_and_store_result("input_validation: character group coerced to factor", coerced)
}

# Test: NA handling ------------------------------------------------------------
test_na_handling_glm <- function() {
  df <- gen_gaussian_data(300, include_g2 = TRUE)
  idx_y <- sample.int(nrow(df), 15); idx_g <- sample.int(nrow(df), 11)
  df$y[idx_y] <- NA
  df$g1[idx_g] <- NA
  
  res <- anova_glm(df, "y", c("g1","g2"), plot_residuals = FALSE)
  # Residuals should correspond to complete cases under model frame
  used_n <- length(res$residuals)
  expect_n <- sum(stats::complete.cases(df[, c("y","g1","g2")]))
  t1 <- identical(used_n, expect_n)
  
  print_and_store_result("na_handling: model uses only complete cases (residual length)", t1)
}

# Test: basic Gaussian fit, Type III, residual plot, R2 -----------------------
test_basic_gaussian <- function() {
  df <- gen_gaussian_data(500, include_g2 = TRUE, sigma = 1.0)
  res <- anova_glm(df, "y", c("g1","g2"),
                   family = gaussian(),
                   sum_squares_type = "III",
                   plot_residuals = TRUE,
                   full_factorial = TRUE,
                   posthoc_method = "none")
  
  need <- c("model","model_stats","anova_test","effect_sizes","confidence_intervals",
            "residuals","residuals_plot","boxplot")
  t1 <- is.list(res) && all(need %in% names(res)) && inherits(res$model, "glm")
  t2 <- !is.null(res$anova_test) && is.data.frame(res$anova_test) && nrow(res$anova_test) >= 1
  t3 <- inherits(res$residuals_plot, "ggplot")
  r2_ok <- is.list(res$effect_sizes) && "r2_gaussian" %in% names(res$effect_sizes) &&
    (is.na(res$effect_sizes$r2_gaussian) || (res$effect_sizes$r2_gaussian >= 0 && res$effect_sizes$r2_gaussian <= 1))
  
  print_and_store_result("gaussian_fit: returns core objects", t1)
  print_and_store_result("gaussian_fit: anova_test present", t2)
  print_and_store_result("gaussian_fit: residuals_plot is ggplot", t3)
  print_and_store_result("gaussian_fit: R2 (gaussian) in [0,1] or NA", r2_ok)
}

# Test: binomial fit, LR tests, non-gaussian notes ----------------------------
test_basic_binomial <- function() {
  df <- gen_binomial_data(600, include_g2 = TRUE)
  res <- anova_glm(df, "y", c("g1","g2"),
                   family = binomial(),
                   sum_squares_type = "III",
                   plot_residuals = FALSE,
                   full_factorial = TRUE,
                   posthoc_method = "none")
  
  t1 <- is.data.frame(res$anova_test) && nrow(res$anova_test) >= 1
  t2 <- is.null(res$residuals_plot)  # we requested FALSE
  t3 <- isTRUE(is.list(res$effect_sizes)) && "pseudo_r2_mcfadden" %in% names(res$effect_sizes)
  t4 <- is.null(res$notes) || any(grepl("non-Gaussian", paste(res$notes, collapse = " ")))
  
  print_and_store_result("binomial_fit: anova_test present (LR)", t1)
  print_and_store_result("binomial_fit: residuals_plot absent when disabled", t2)
  print_and_store_result("binomial_fit: pseudo R2 returned", t3)
  print_and_store_result("binomial_fit: non-Gaussian note present", t4)
}

# Test: Poisson with overdispersion note --------------------------------------
test_poisson_overdispersion_note <- function() {
  set.seed(42)
  df <- gen_overdispersed_counts(600, theta = 0.6)
  res <- anova_glm(df, "y", c("g1","g2"),
                   family = poisson(),
                   sum_squares_type = "II",
                   plot_residuals = FALSE,
                   full_factorial = TRUE,
                   posthoc_method = "none")
  msg <- paste(res$notes, collapse = " | ")
  flagged <- any(grepl("overdispersion", tolower(msg)))
  print_and_store_result("poisson_overdispersion: note indicates possible overdispersion", flagged)
}

# Test: Type II vs Type III (unbalanced factorial) ----------------------------
test_type_ii_vs_iii <- function() {
  df <- gen_unbalanced_gaussian(150, 90)
  res2 <- anova_glm(df, "y", c("g1","g2"),
                    family = gaussian(), sum_squares_type = "II",
                    plot_residuals = FALSE, full_factorial = TRUE, posthoc_method = "none")
  res3 <- anova_glm(df, "y", c("g1","g2"),
                    family = gaussian(), sum_squares_type = "III",
                    plot_residuals = FALSE, full_factorial = TRUE, posthoc_method = "none")
  
  t1 <- is.data.frame(res2$anova_test) && is.data.frame(res3$anova_test)
  rn2 <- rownames(res2$anova_test); rn3 <- rownames(res3$anova_test)
  has_terms_ii  <- any(grepl("g1", rn2)) && any(grepl("g2", rn2))
  has_terms_iii <- any(grepl("g1", rn3)) && any(grepl("g2", rn3))
  print_and_store_result("type_ii_vs_iii: both ANOVA tables present", t1)
  print_and_store_result("type_ii_vs_iii: terms detected in both tables", has_terms_ii && has_terms_iii)
}

# Test: post-hoc via multcomp (single factor, >2 levels) ----------------------
test_posthoc_multcomp_single_factor <- function() {
  df <- gen_gaussian_data(450, include_g2 = FALSE, sigma = 1.2)  # 3 levels in g1
  res <- anova_glm(df, "y", "g1",
                   family = gaussian(),
                   sum_squares_type = "II",
                   plot_residuals = FALSE,
                   full_factorial = FALSE,
                   posthoc_method = "multcomp")
  t1 <- is.null(res$posthoc) || inherits(res$posthoc, "summary.glht")
  print_and_store_result("posthoc_multcomp: returns summary.glht or NULL (if degenerate)", t1)
}

# Test: post-hoc via emmeans (full factorial) ---------------------------------
test_posthoc_emmeans_full_factorial <- function() {
  df <- gen_gaussian_data(480, include_g2 = TRUE, sigma = 0.9)
  res <- anova_glm(df, "y", c("g1","g2"),
                   family = gaussian(),
                   sum_squares_type = "III",
                   plot_residuals = FALSE,
                   full_factorial = TRUE,
                   posthoc_method = "emmeans",
                   emmeans_simple_for = "g1")
  
  ok <- is.null(res$posthoc) ||
    (is.list(res$posthoc) && all(c("table","emm_obj") %in% names(res$posthoc)))
  print_and_store_result("posthoc_emmeans: returns table & emm_obj (or NULL if unavailable)", ok)
}

# Test: boxplot object structure ----------------------------------------------
test_boxplot_object <- function() {
  df <- gen_poisson_data(320, include_g2 = TRUE)
  res <- anova_glm(df, "y", c("g1","g2"),
                   family = poisson(),
                   sum_squares_type = "II",
                   plot_residuals = FALSE,
                   full_factorial = FALSE)  # combined groups
  t1 <- inherits(res$boxplot, "ggplot") || inherits(res$boxplot, "gg")
  gb <- ggplot_build(res$boxplot)
  layer_idx <- which(lengths(gb$data) > 0)[1]
  n_found <- length(unique(gb$data[[layer_idx]]$x))
  n_expected <- nlevels(interaction(df$g1, df$g2))
  t2 <- n_found == n_expected
  print_and_store_result("boxplot: returns ggplot object", t1)
  print_and_store_result("boxplot: x categories match combined groups", t2)
}

# Test: interaction (combined vs full_factorial) ------------------------------
test_interaction_modes <- function() {
  df <- gen_gaussian_data(400, include_g2 = TRUE)
  res_combined <- anova_glm(df, "y", c("g1","g2"),
                            family = gaussian(),
                            sum_squares_type = "II",
                            full_factorial = FALSE,
                            posthoc_method = "none")
  res_full <- anova_glm(df, "y", c("g1","g2"),
                        family = gaussian(),
                        sum_squares_type = "II",
                        full_factorial = TRUE,
                        posthoc_method = "none")
  rn_combined <- rownames(res_combined$anova_test)
  rn_full <- rownames(res_full$anova_test)
  t1 <- any(grepl("interaction_term", rn_combined))
  t2 <- any(grepl("g1:g2", rn_full)) || any(grepl("g1\\:g2", rn_full))
  print_and_store_result("interaction_modes: combined mode includes 'interaction_term'", t1)
  print_and_store_result("interaction_modes: full factorial includes 'g1:g2'", t2)
}

# Test: confidence intervals table present ------------------------------------
test_confidence_intervals_present <- function() {
  df <- gen_gaussian_data(300, include_g2 = TRUE)
  res <- anova_glm(df, "y", c("g1","g2"),
                   family = gaussian(),
                   sum_squares_type = "III",
                   full_factorial = TRUE,
                   posthoc_method = "none")
  t1 <- is.data.frame(res$confidence_intervals) && ncol(res$confidence_intervals) == 2
  print_and_store_result("confint: confidence intervals data frame present", t1)
}

# Test: residuals are finite and length matches df.residual -------------------
test_residuals_integrity <- function() {
  df <- gen_gaussian_data(350, include_g2 = TRUE)  # n not multiple of 4; generators handle it
  res <- anova_glm(df, "y", c("g1","g2"),
                   family = gaussian(),
                   sum_squares_type = "II",
                   full_factorial = TRUE)
  t1 <- all(is.finite(res$residuals))
  t2 <- length(res$residuals) == res$model$df.residual + length(coef(res$model))
  print_and_store_result("residuals: all finite", t1)
  print_and_store_result("residuals: length equals observations used", t2)
}

# End-to-end scenarios ---------------------------------------------------------
test_end_to_end_glm <- function() {
  # Scenario 1: One factor Gaussian, Type II, multcomp post-hoc
  df1 <- gen_gaussian_data(500, include_g2 = FALSE)
  res1 <- anova_glm(df1, "y", "g1",
                    family = gaussian(),
                    sum_squares_type = "II",
                    full_factorial = FALSE,
                    posthoc_method = "multcomp")
  need1 <- c("model","anova_test","confidence_intervals","boxplot","effect_sizes")
  t1 <- is.list(res1) && all(need1 %in% names(res1)) && inherits(res1$model, "glm")
  
  # Scenario 2: Two factors binomial, Type III, emmeans simple-effects
  df2 <- gen_binomial_data(700, include_g2 = TRUE)
  res2 <- anova_glm(df2, "y", c("g1","g2"),
                    family = binomial(),
                    sum_squares_type = "III",
                    full_factorial = TRUE,
                    posthoc_method = "emmeans",
                    emmeans_simple_for = "g1")
  t2 <- is.data.frame(res2$anova_test) &&
    (is.null(res2$posthoc) || (is.list(res2$posthoc) && "table" %in% names(res2$posthoc)))
  
  # Scenario 3: Two factors Poisson, Type II, no post-hoc, residuals plot
  df3 <- gen_poisson_data(600, include_g2 = TRUE)
  res3 <- anova_glm(df3, "y", c("g1","g2"),
                    family = poisson(),
                    sum_squares_type = "II",
                    full_factorial = TRUE,
                    plot_residuals = TRUE,
                    posthoc_method = "none")
  t3 <- inherits(res3$residuals_plot, "ggplot")
  
  print_and_store_result("end_to_end: Gaussian one-factor + multcomp", t1)
  print_and_store_result("end_to_end: Binomial two-factor + emmeans", t2)
  print_and_store_result("end_to_end: Poisson two-factor + residuals plot", t3)
}

###############################################################################
# Test Runner
###############################################################################

#' Run All anova_glm Tests
#'
#' Execute a comprehensive suite of tests for the `anova_glm()` function.
#'
#' @details
#' All tests write a PASS/FAIL line and populate a global \code{test_results} table.
run_all_tests_anova_glm <- function() {
  cat("========== Running Comprehensive anova_glm UAT ==========\n")
  test_input_validation_glm()
  test_na_handling_glm()
  test_basic_gaussian()
  test_basic_binomial()
  test_poisson_overdispersion_note()
  test_type_ii_vs_iii()
  test_posthoc_multcomp_single_factor()
  test_posthoc_emmeans_full_factorial()
  test_boxplot_object()
  test_interaction_modes()
  test_confidence_intervals_present()
  test_residuals_integrity()
  test_end_to_end_glm()
  cat("========== Tests Completed ==========\n\n")
  cat("Test Summary:\n")
  print(table(test_results$Result))
  cat("\nDetailed Results:\n")
  print(test_results)
}

# Uncomment to run all tests:
# run_all_tests_anova_glm()
