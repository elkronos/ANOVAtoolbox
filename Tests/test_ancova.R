###############################################################################
# Testing Infrastructure - ANCOVA Toolkit
###############################################################################

# Package management for tests -------------------------------------------------
required_packages <- c("car", "ggplot2", "emmeans", "withr", "sandwich", "lmtest")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
suppressPackageStartupMessages({
  library(car)
  library(ggplot2)
  library(emmeans)
  library(withr)
  library(sandwich)
  library(lmtest)
})

# Global container for test results -------------------------------------------
test_results <- data.frame(Test = character(), Result = character(), stringsAsFactors = FALSE)

#' Helper: Print and Store Test Result
#'
#' @param test_name Character. Name of the test.
#' @param passed Logical. Whether the test passed.
#' @param note Optional character. Additional notes.
print_and_store_result <- function(test_name, passed, note = NULL) {
  result <- if (isTRUE(passed)) "PASS" else "FAIL"
  cat(sprintf("%-65s [%s]\n", test_name, result))
  if (!is.null(note)) cat("  Note: ", note, "\n", sep = "")
  test_results <<- rbind(
    test_results,
    data.frame(Test = test_name, Result = result, stringsAsFactors = FALSE)
  )
}

###############################################################################
# Unit Tests
###############################################################################

# Test: validate_vars ----------------------------------------------------------
test_validate_vars <- function() {
  df <- data.frame(a = 1:3, b = 4:6)
  ok <- tryCatch({ validate_vars(df, c("a","b")); TRUE }, error = function(e) FALSE)
  missing <- tryCatch({ validate_vars(df, c("a","c")); FALSE }, error = function(e) TRUE)
  print_and_store_result("validate_vars: all variables present", ok)
  print_and_store_result("validate_vars: missing variable detected", missing)
}

# Test: ensure_factor_iv -------------------------------------------------------
test_ensure_factor_iv <- function() {
  df <- data.frame(y = rnorm(6), g = c(1,1,1,2,2,2))
  df2 <- tryCatch(ensure_factor_iv(df, "g"), error = function(e) e)
  coerced <- is.data.frame(df2) && is.factor(df2$g) && nlevels(df2$g) == 2
  single_level_error <- inherits(tryCatch({
    ensure_factor_iv(data.frame(y = 1:3, g = c(1,1,1)), "g")
  }, error = function(e) e), "error")
  print_and_store_result("ensure_factor_iv: coerces to factor and checks >= 2 levels", coerced)
  print_and_store_result("ensure_factor_iv: errors for single-level factor", single_level_error)
}

# Test: check_homogeneity_slopes ----------------------------------------------
test_check_homogeneity_slopes <- function() {
  set.seed(42)
  n <- 120
  g <- factor(rep(c("A","B"), each = n/2))
  x <- rnorm(n)
  # Equal slopes case: y = 2 + 1.5*x + group effect
  y_equal <- 2 + 1.5*x + ifelse(g=="B", 0.5, 0) + rnorm(n, sd = 0.5)
  d_equal <- data.frame(y = y_equal, g = g, x = x)
  
  # Different slopes case: y = 2 + 1.5*x + 1.0*(x*I[g==\"B\"]) + group effect
  y_diff <- 2 + 1.5*x + ifelse(g=="B", 1.0*x, 0) + ifelse(g=="B", 0.5, 0) + rnorm(n, sd = 0.5)
  d_diff <- data.frame(y = y_diff, g = g, x = x)
  
  res_equal <- check_homogeneity_slopes(d_equal, "y", "g", "x", alpha = 0.05)
  res_diff  <- check_homogeneity_slopes(d_diff,  "y", "g", "x", alpha = 0.05)
  
  t1 <- is.list(res_equal) && isTRUE(res_equal$homogeneous)
  t2 <- is.list(res_diff)  && identical(res_diff$homogeneous, FALSE)
  t3 <- is.numeric(res_equal$p_value) && is.numeric(res_diff$p_value)
  
  print_and_store_result("check_homogeneity_slopes: detects equal slopes", t1)
  print_and_store_result("check_homogeneity_slopes: detects different slopes", t2)
  print_and_store_result("check_homogeneity_slopes: returns numeric p-values", t3)
}

# Test: check_normality --------------------------------------------------------
test_check_normality <- function() {
  # Case 1: n within [3, 5000]
  m_ok <- lm(mpg ~ wt, data = mtcars)
  res_ok <- check_normality(m_ok)
  case_ok <- inherits(res_ok$shapiro, "htest") && inherits(res_ok$qq_plot, "gg")
  
  # Case 2: n < 3
  m_small <- lm(y ~ 1, data = data.frame(y = c(1, 2)))  # residuals length = 2
  res_small <- check_normality(m_small)
  case_small <- is.null(res_small$shapiro) && inherits(res_small$qq_plot, "gg")
  
  # Case 3: n > 5000
  set.seed(7)
  big_n <- 5001
  df_big <- data.frame(y = rnorm(big_n), x = rnorm(big_n))
  m_big <- lm(y ~ x, data = df_big)
  res_big <- check_normality(m_big)
  case_big <- is.null(res_big$shapiro) && inherits(res_big$qq_plot, "gg")
  
  print_and_store_result("check_normality: shapiro + QQ for n in [3,5000]", case_ok)
  print_and_store_result("check_normality: skips shapiro for n < 3", case_small)
  print_and_store_result("check_normality: skips shapiro for n > 5000", case_big)
}

# Test: check_homogeneity_variance --------------------------------------------
test_check_homogeneity_variance <- function() {
  m <- lm(mpg ~ factor(cyl), data = mtcars)
  lev <- check_homogeneity_variance(mtcars, m, "cyl")
  ok <- "Pr(>F)" %in% colnames(lev)
  print_and_store_result("check_homogeneity_variance: returns Levene table", ok)
}

# Test: anova_type3 ------------------------------------------------------------
test_anova_type3 <- function() {
  m <- withr::with_options(
    list(contrasts = c("contr.sum", "contr.poly")),
    lm(mpg ~ wt + factor(cyl), data = mtcars)
  )
  a_classical <- anova_type3(m, use_hc = FALSE)
  a_hc        <- anova_type3(m, use_hc = TRUE, hc_type = "HC3")
  
  t1 <- is.list(a_classical) && "classical" %in% names(a_classical)
  t2 <- is.null(a_classical$hc_tests)
  t3 <- is.list(a_hc) && is.data.frame(a_hc$hc_tests)
  t4 <- all(c("Pr(>F)") %in% colnames(a_hc$hc_tests))
  
  print_and_store_result("anova_type3: classical Type-III table present", t1)
  print_and_store_result("anova_type3: hc_tests NULL when use_hc=FALSE", t2)
  print_and_store_result("anova_type3: hc_tests returned when use_hc=TRUE", t3)
  print_and_store_result("anova_type3: hc_tests contains p-values", t4)
}

# Test: calculate_effect_sizes -------------------------------------------------
test_calculate_effect_sizes <- function() {
  m <- withr::with_options(
    list(contrasts = c("contr.sum", "contr.poly")),
    lm(mpg ~ wt + factor(cyl), data = mtcars)
  )
  es <- calculate_effect_sizes(m)
  ok <- is.data.frame(es) &&
    all(c("Effect","Partial_Eta_Squared") %in% names(es)) &&
    all(es$Partial_Eta_Squared >= 0 & es$Partial_Eta_Squared <= 1)
  print_and_store_result("calculate_effect_sizes: returns valid partial eta squared", ok)
}

# Test: plot_diagnostics -------------------------------------------------------
test_plot_diagnostics <- function() {
  m <- lm(mpg ~ wt + factor(cyl), data = mtcars)
  pl <- plot_diagnostics(m, mtcars, "mpg", "cyl", "wt")
  ok <- is.list(pl) &&
    all(c("residuals_vs_fitted","qq_plot","ancova_plot") %in% names(pl)) &&
    all(vapply(pl, function(p) inherits(p, "gg"), logical(1)))
  print_and_store_result("plot_diagnostics: returns three ggplot objects", ok)
}

# Test: compute_emm ------------------------------------------------------------
test_compute_emm <- function() {
  m <- withr::with_options(
    list(contrasts = c("contr.sum", "contr.poly")),
    lm(mpg ~ wt + factor(cyl), data = mtcars)
  )
  e1 <- tryCatch(compute_emm(m, "cyl", use_hc = FALSE), error = function(e) e)
  e2 <- tryCatch(compute_emm(m, "cyl", use_hc = TRUE, hc_type = "HC3"), error = function(e) e)
  
  ok1 <- inherits(e1, "emmGrid")
  ok2 <- inherits(e2, "emmGrid")
  
  print_and_store_result("compute_emm: EMMs without HC covariance", ok1)
  print_and_store_result("compute_emm: EMMs with HC covariance", ok2)
}

# Test: ancova_analysis (end-to-end) ------------------------------------------
test_ancova_analysis <- function() {
  set.seed(123)
  n <- 150
  g <- factor(rep(c("G1","G2","G3"), length.out = n))
  x <- rnorm(n, 0, 1)
  # Construct data with additive structure (no interaction)
  y <- 50 + 5*x + ifelse(g=="G2", 2, ifelse(g=="G3", -1, 0)) + rnorm(n, sd = 2)
  d <- data.frame(dv = y, iv = g, covariate = x)
  
  # Add a few NAs to exercise drop_na branch
  d$dv[sample.int(n, 4)] <- NA
  d$covariate[sample.int(n, 3)] <- NA
  
  res <- ancova_analysis(
    d, dv = "dv", iv = "iv", covariate = "covariate",
    alpha = 0.05, use_hc = TRUE, hc_type = "HC3",
    drop_na = TRUE, plots = TRUE, verbose = FALSE
  )
  
  # Basic structure checks
  needed <- c("final_model","homogeneity_slopes","normality","levene_test",
              "anova_type3","effect_sizes","emm","plots")
  has_all <- is.list(res) && all(needed %in% names(res))
  
  # Model object check
  is_lm <- inherits(res$final_model, "lm")
  
  # ANOVA content check
  has_classical <- is.list(res$anova_type3) && "classical" %in% names(res$anova_type3)
  
  # Effect sizes are in range
  es <- res$effect_sizes
  es_ok <- is.data.frame(es) &&
    all(es$Partial_Eta_Squared >= 0 & es$Partial_Eta_Squared <= 1)
  
  # EMMs present
  emm_ok <- inherits(res$emm, "emmGrid") || is.null(res$emm)
  
  # Plots present
  plots_ok <- is.list(res$plots) &&
    all(c("residuals_vs_fitted","qq_plot","ancova_plot") %in% names(res$plots))
  
  print_and_store_result("ancova_analysis: returns named elements", has_all)
  print_and_store_result("ancova_analysis: final model is lm", is_lm)
  print_and_store_result("ancova_analysis: includes Type-III table", has_classical)
  print_and_store_result("ancova_analysis: effect sizes in [0,1]", es_ok)
  print_and_store_result("ancova_analysis: EMMs object present or NULL", emm_ok)
  print_and_store_result("ancova_analysis: plots list contains three ggplots", plots_ok)
  
  # Additional scenario: interaction selected when slopes differ
  set.seed(456)
  n2 <- 120
  g2 <- factor(rep(c("A","B"), each = n2/2))
  x2 <- rnorm(n2)
  y2 <- 10 + 2*x2 + ifelse(g2=="B", 1.5*x2, 0) + rnorm(n2, sd = 1)
  d2 <- data.frame(dv = y2, iv = g2, covariate = x2)
  res2 <- ancova_analysis(
    d2, "dv", "iv", "covariate",
    alpha = 0.05, use_hc = FALSE, drop_na = TRUE, plots = FALSE, verbose = FALSE
  )
  # Expect interaction in formula
  ftxt <- deparse(formula(res2$final_model))
  has_interaction <- grepl("covariate \\* iv", ftxt, fixed = FALSE)
  print_and_store_result("ancova_analysis: selects interaction when slopes differ", has_interaction)
}

###############################################################################
# Test Runner
###############################################################################

#' Run All Tests
#'
#' Execute a suite of tests to assess the ANCOVA Toolkit.
#'
#' @param include_large Logical. If TRUE, include the n > 5000 normality case.
#' @details
#' All tests write a PASS/FAIL line and populate a global \code{test_results} table.
run_all_tests <- function(include_large = TRUE) {
  cat("========== Running Comprehensive ANCOVA Toolkit Tests ==========\n")
  test_validate_vars()
  test_ensure_factor_iv()
  test_check_homogeneity_slopes()
  test_check_normality()
  test_check_homogeneity_variance()
  test_anova_type3()
  test_calculate_effect_sizes()
  test_plot_diagnostics()
  test_compute_emm()
  test_ancova_analysis()
  if (!include_large) {
    # Already handled inside test_check_normality; keep parameter for symmetry with other harnesses
  }
  cat("========== Tests Completed ==========\n\n")
  cat("Test Summary:\n")
  print(table(test_results$Result))
  cat("\nDetailed Results:\n")
  print(test_results)
}

# Uncomment to run all tests:
# run_all_tests(include_large = TRUE)
