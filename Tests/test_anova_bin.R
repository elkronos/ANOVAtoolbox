###############################################################################
# Testing Infrastructure - anova_bin UAT (Logistic "One-Way ANOVA")
###############################################################################

# Package management for tests -------------------------------------------------
required_packages <- c(
  "car","ggplot2","dplyr","broom","withr","sandwich","lmtest",
  "rlang","scales","magrittr"
)
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
suppressPackageStartupMessages({
  library(car)
  library(ggplot2)
  library(dplyr)
  library(broom)
  library(withr)
  library(sandwich)
  library(lmtest)
  library(rlang)
  library(scales)
  library(magrittr)
})

# NOTE: This UAT assumes `anova_bin()` is already defined in the session.

# Global container for test results -------------------------------------------
test_results <- data.frame(Test = character(), Result = character(), stringsAsFactors = FALSE)

#' Helper: Print and Store Test Result
#'
#' @param test_name Character. Name of the test.
#' @param passed Logical. Whether the test passed.
#' @param note Optional character. Additional notes.
print_and_store_result <- function(test_name, passed, note = NULL) {
  result <- if (isTRUE(passed)) "PASS" else "FAIL"
  cat(sprintf("%-80s [%s]\n", test_name, result))
  if (!is.null(note)) cat("  Note: ", note, "\n", sep = "")
  test_results <<- rbind(
    test_results,
    data.frame(Test = test_name, Result = result, stringsAsFactors = FALSE)
  )
}

###############################################################################
# Synthetic Data Generators
###############################################################################

# Balanced two-group binary outcome with signal in group
gen_data_two_group <- function(n = 200, p_success_A = 0.3, p_success_B = 0.7) {
  g <- factor(rep(c("A","B"), each = n/2))
  prob <- ifelse(g == "A", p_success_A, p_success_B)
  y <- rbinom(n, size = 1, prob = prob)
  data.frame(y = factor(y), g = g)
}

# Two factors with interaction on log-odds scale
gen_data_two_factors <- function(n = 400) {
  g1 <- factor(rep(c("A","B"), each = n/2))
  g2 <- factor(rep(rep(c("X","Y"), each = n/4), 2))
  # baseline logit = -0.5, main effects and interaction
  eta <- -0.5 + 0.8*(g1=="B") + 0.6*(g2=="Y") + 0.9*(g1=="B" & g2=="Y")
  p <- 1/(1+exp(-eta))
  y <- rbinom(n, 1, p)
  data.frame(y = factor(y), g1 = g1, g2 = g2)
}

###############################################################################
# Unit Tests
###############################################################################

# Test: input validation -------------------------------------------------------
test_input_validation <- function() {
  df <- gen_data_two_group(100)
  # Missing response
  e1 <- inherits(tryCatch(anova_bin(df, "y_missing", "g"), error = identity), "error")
  # Missing group
  e2 <- inherits(tryCatch(anova_bin(df, "y", "g_missing"), error = identity), "error")
  # Non-binary response (3+ levels)
  df_bad <- df; df_bad$y <- factor(c(0,1,2))[sample(1:3, nrow(df_bad), TRUE)]
  e3 <- inherits(tryCatch(anova_bin(df_bad, "y", "g"), error = identity), "error")
  # success_level not in levels
  e4 <- inherits(tryCatch(anova_bin(df, "y", "g", success_level = "notalevel"), error = identity), "error")
  # group_ref level not in factor
  e5 <- inherits(tryCatch(anova_bin(df, "y", "g", group_ref = list(g="Z")), error = identity), "error")
  
  print_and_store_result("input_validation: errors when response variable missing", e1)
  print_and_store_result("input_validation: errors when group variable missing", e2)
  print_and_store_result("input_validation: errors when response not binary", e3)
  print_and_store_result("input_validation: errors when success_level invalid", e4)
  print_and_store_result("input_validation: errors when group_ref invalid", e5)
}

# Test: NA handling ------------------------------------------------------------
test_na_handling <- function() {
  df <- gen_data_two_group(200)
  # Insert NAs
  idx_y <- sample.int(nrow(df), 10); idx_g <- sample.int(nrow(df), 7)
  df$y[idx_y] <- NA
  df$g[idx_g] <- NA
  
  # With na.rm = TRUE: rows dropped equals complete.cases count
  res_true <- anova_bin(df, "y", "g", na.rm = TRUE, print_plot = FALSE)
  used_n_true <- length(res_true$fitted_model$y)
  expect_n_true <- sum(stats::complete.cases(df[, c("y","g")]))
  t1 <- identical(used_n_true, expect_n_true)
  
  # With na.rm = FALSE: glm drops internally; prop_table includes NA groups
  res_false <- anova_bin(df, "y", "g", na.rm = FALSE, print_plot = FALSE)
  used_n_false <- length(res_false$fitted_model$y)
  expect_n_false <- sum(stats::complete.cases(df[, c("y","g")]))
  # proportions table should have some NA in g or y
  has_na_group <- any(is.na(res_false$prop_table$g)) || any(is.na(res_false$prop_table$y))
  t2 <- identical(used_n_false, expect_n_false)
  t3 <- isTRUE(has_na_group)
  
  print_and_store_result("na_handling: na.rm=TRUE uses only complete cases", t1)
  print_and_store_result("na_handling: na.rm=FALSE still fits on complete cases", t2)
  print_and_store_result("na_handling: na.rm=FALSE retains NA groups in prop_table", t3)
}

# Test: success_level handling -------------------------------------------------
test_success_level <- function() {
  set.seed(123)
  df <- gen_data_two_group(300, p_success_A = 0.2, p_success_B = 0.8)
  
  # success = "1" (default second level for factor(y))
  res1 <- anova_bin(df, "y", "g", print_plot = FALSE)
  # success = "0" flips the modeled event
  res0 <- anova_bin(df, "y", "g", success_level = "0", print_plot = FALSE)
  
  # Check model_stats success level
  t1 <- identical(res1$model_stats$success_level, "1")
  t2 <- identical(res0$model_stats$success_level, "0")
  
  # Effect of gB should invert between success levels
  or1 <- exp(coef(res1$fitted_model)[["gB"]])
  or0 <- exp(coef(res0$fitted_model)[["gB"]])
  t3 <- isTRUE(or1 > 1) && isTRUE(or0 < 1)
  
  print_and_store_result("success_level: default success is second level", t1)
  print_and_store_result("success_level: respects provided success_level", t2)
  print_and_store_result("success_level: odds ratio direction flips with success level", t3)
}

# Test: group_ref baseline setting --------------------------------------------
test_group_ref <- function() {
  df <- gen_data_two_group(200)
  # Default baseline is first level; force baseline = "B" and compare sign/inversion
  res_default <- anova_bin(df, "y", "g", print_plot = FALSE)
  res_refB    <- anova_bin(df, "y", "g", group_ref = list(g = "B"), print_plot = FALSE)
  
  # Term name changes when baseline changes; check that coefficients correspond
  coef_def <- coef(res_default$fitted_model)
  coef_ref <- coef(res_refB$fitted_model)
  
  # With baseline "B", the coefficient is for gA (vs B), which should be the negative of gB (vs A).
  t1 <- all(c("(Intercept)","gB") %in% names(coef_def))
  t2 <- all(c("(Intercept)","gA") %in% names(coef_ref))
  t3 <- isTRUE(all.equal(unname(coef_ref["gA"]), -unname(coef_def["gB"]), tolerance = 1e-6))
  
  print_and_store_result("group_ref: default coefficient term exists", t1)
  print_and_store_result("group_ref: ref change produces complementary term", t2)
  print_and_store_result("group_ref: coefficients invert when changing reference", t3)
}

# Test: ANOVA types (I / II / III) --------------------------------------------
test_anova_types <- function() {
  df <- gen_data_two_factors(400)
  
  # Type I (sequential)
  res_I <- anova_bin(df, "y", c("g1","g2"), anova_type = "I", print_plot = FALSE)
  tI <- !is.null(res_I$deviance_anova) && is.null(res_I$car_anova)
  
  # Type II via car::Anova
  res_II <- withr::with_options(
    list(contrasts = c("contr.sum","contr.poly")),
    anova_bin(df, "y", c("g1","g2"), anova_type = "II", anova_test = "LR", print_plot = FALSE)
  )
  tII <- is.null(res_II$deviance_anova) && inherits(res_II$car_anova, "anova")
  
  # Type III via car::Anova
  res_III <- withr::with_options(
    list(contrasts = c("contr.sum","contr.poly")),
    anova_bin(df, "y", c("g1","g2"), anova_type = "III", anova_test = "Wald", print_plot = FALSE)
  )
  tIII <- is.null(res_III$deviance_anova) && inherits(res_III$car_anova, "anova")
  
  print_and_store_result("anova_types: Type I returns deviance_anova", tI)
  print_and_store_result("anova_types: Type II returns car_anova table", tII)
  print_and_store_result("anova_types: Type III returns car_anova table", tIII)
}

# Test: robust SE path ---------------------------------------------------------
test_robust_se <- function() {
  set.seed(42)
  df <- gen_data_two_factors(300)
  
  # Conventional
  res_classic <- anova_bin(df, "y", c("g1","g2"), robust = FALSE, print_plot = FALSE)
  # Robust HC3
  res_robust  <- anova_bin(df, "y", c("g1","g2"), robust = TRUE, robust_type = "HC3", print_plot = FALSE)
  
  # Check fields and that robust flag is recorded
  has_coef <- is.data.frame(res_robust$coef_tests) &&
    all(c("term","estimate","std_error","z_value","p_value") %in% names(res_robust$coef_tests))
  flag_set <- isTRUE(res_robust$model_stats$robust) &&
    identical(res_robust$model_stats$robust_type, "HC3")
  
  # Effect size table numeric columns present
  es_ok <- is.data.frame(res_robust$effect_size) &&
    all(c("term","odds_ratio","conf.low","conf.high","p_value") %in% names(res_robust$effect_size))
  
  # Compare at least one std error differs in many cases (not guaranteed but often true)
  try_diff <- try({
    # Align by term names (may differ in ordering)
    ct1 <- res_classic$coef_tests; rownames(ct1) <- ct1$term
    ct2 <- res_robust$coef_tests;  rownames(ct2) <- ct2$term
    common <- intersect(rownames(ct1), rownames(ct2))
    any_diff <- any(abs(ct1[common, "std_error"] - ct2[common, "std_error"]) > 1e-8)
    any_diff
  }, silent = TRUE)
  diff_flag <- isTRUE(try_diff)
  
  print_and_store_result("robust_se: coef_tests present with required columns", has_coef)
  print_and_store_result("robust_se: model_stats records robust flags", flag_set)
  print_and_store_result("robust_se: effect_size table structure OK", es_ok)
  print_and_store_result("robust_se: robust vs classic std errors differ (often)", diff_flag)
}

# Test: LR null vs full comparison --------------------------------------------
test_lr_comparison <- function() {
  df <- gen_data_two_group(300, 0.25, 0.75)
  res <- anova_bin(df, "y", "g", print_plot = FALSE)
  
  ok_models <- is.data.frame(res$lr_comparison$models) &&
    all(c("Model","logLik","df") %in% names(res$lr_comparison$models))
  ok_test <- is.data.frame(res$lr_comparison$test) &&
    all(c("Test","Chisq","Df","Pr(>Chisq)") %in% names(res$lr_comparison$test)) &&
    res$lr_comparison$test$`Pr(>Chisq)` >= 0 && res$lr_comparison$test$`Pr(>Chisq)` <= 1
  
  print_and_store_result("lr_comparison: contains models table", ok_models)
  print_and_store_result("lr_comparison: contains valid LR test table", ok_test)
}

# Test: effect size & intercept toggle ----------------------------------------
test_effect_size_intercept_toggle <- function() {
  df <- gen_data_two_group(200)
  res_no_int <- anova_bin(df, "y", "g", include_intercept = FALSE, print_plot = FALSE)
  res_int    <- anova_bin(df, "y", "g", include_intercept = TRUE,  print_plot = FALSE)
  
  t1 <- all(res_no_int$effect_size$term != "(Intercept)")
  t2 <- any(res_int$effect_size$term == "(Intercept)")
  t3 <- all(res_int$effect_size$odds_ratio > 0)
  
  print_and_store_result("effect_size: excludes intercept when include_intercept=FALSE", t1)
  print_and_store_result("effect_size: includes intercept when include_intercept=TRUE", t2)
  print_and_store_result("effect_size: odds ratios are positive", t3)
}

# Test: proportions table & plot object ---------------------------------------
test_props_and_plot <- function() {
  df <- gen_data_two_factors(240)
  
  res1 <- anova_bin(df, "y", "g1", print_plot = FALSE)
  res2 <- anova_bin(df, "y", c("g1","g2"), print_plot = FALSE)
  
  # prop_table sums to 1 per group combination
  sum_check_1 <- res1$prop_table %>%
    group_by(g1) %>% summarise(s = sum(prop)) %>% pull(s)
  sum_check_2 <- res2$prop_table %>%
    group_by(g1, g2) %>% summarise(s = sum(prop)) %>% pull(s)
  
  t1 <- inherits(res1$count_plot, "gg")
  t2 <- all(abs(sum_check_1 - 1) < 1e-8)
  t3 <- all(abs(sum_check_2 - 1) < 1e-8)
  
  print_and_store_result("prop_table/plot: count_plot is a ggplot object", t1)
  print_and_store_result("prop_table/plot: proportions sum to 1 by g1", t2)
  print_and_store_result("prop_table/plot: proportions sum to 1 by g1:g2", t3)
}

# Test: multi-group plotting (combined label) ---------------------------------
test_multi_group_plot_label <- function() {
  df <- gen_data_two_factors(200)
  # The plot internally builds a "group_combined" aesthetic for 2+ groups; we cannot
  # directly access it, but we can ensure the number of x categories equals combinations.
  res <- anova_bin(df, "y", c("g1","g2"), print_plot = FALSE)
  # compute expected unique combos
  n_expected <- nlevels(df$g1) * nlevels(df$g2)
  # Extract data from ggplot build
  gb <- ggplot_build(res$count_plot)
  # number of x positions found in plot = expected (stacked bars per combo)
  n_found <- length(unique(gb$data[[1]]$x))
  print_and_store_result("multi_group_plot: x categories equal group combinations", n_found == n_expected)
}

# Test: print_plot flag --------------------------------------------------------
test_print_plot_flag <- function() {
  df <- gen_data_two_group(150)
  # Should not error when print_plot = FALSE
  ok <- !inherits(tryCatch(anova_bin(df, "y", "g", print_plot = FALSE), error = identity), "error")
  print_and_store_result("print_plot: no error when print_plot=FALSE", ok)
}

# End-to-end scenarios ---------------------------------------------------------
test_end_to_end <- function() {
  # Scenario 1: Single factor, Type I, robust HC0
  df1 <- gen_data_two_group(300, 0.3, 0.8)
  res1 <- anova_bin(df1, "y", "g", robust = TRUE, robust_type = "HC0", print_plot = FALSE)
  need1 <- c("deviance_anova","lr_comparison","coef_tests","effect_size",
             "model_stats","prop_table","count_plot","fitted_model")
  t1 <- is.list(res1) && all(need1 %in% names(res1)) && inherits(res1$fitted_model, "glm")
  
  # Scenario 2: Two factors, Type II (LR), success level explicitly set
  df2 <- gen_data_two_factors(500)
  res2 <- withr::with_options(
    list(contrasts = c("contr.sum","contr.poly")),
    anova_bin(df2, "y", c("g1","g2"), anova_type = "II", anova_test = "LR",
              success_level = "1", print_plot = FALSE)
  )
  t2 <- is.null(res2$deviance_anova) && inherits(res2$car_anova, "anova") &&
    identical(res2$model_stats$success_level, "1")
  
  # Scenario 3: Two factors, Type III (Wald), group references set, include_intercept
  res3 <- withr::with_options(
    list(contrasts = c("contr.sum","contr.poly")),
    anova_bin(df2, "y", c("g1","g2"),
              anova_type = "III", anova_test = "Wald",
              group_ref = list(g1 = "B", g2 = "Y"),
              include_intercept = TRUE, print_plot = FALSE)
  )
  t3 <- is.data.frame(res3$effect_size) && any(res3$effect_size$term == "(Intercept)")
  
  print_and_store_result("end_to_end: Scenario 1 returns all required objects", t1)
  print_and_store_result("end_to_end: Scenario 2 uses Type II (LR) with success level set", t2)
  print_and_store_result("end_to_end: Scenario 3 Type III with refs and intercept", t3)
}

###############################################################################
# Test Runner
###############################################################################

#' Run All anova_bin Tests
#'
#' Execute a comprehensive suite of tests for the `anova_bin()` function.
#'
#' @details
#' All tests write a PASS/FAIL line and populate a global \code{test_results} table.
run_all_tests_anova_bin <- function() {
  cat("========== Running Comprehensive anova_bin UAT ==========\n")
  test_input_validation()
  test_na_handling()
  test_success_level()
  test_group_ref()
  test_anova_types()
  test_robust_se()
  test_lr_comparison()
  test_effect_size_intercept_toggle()
  test_props_and_plot()
  test_multi_group_plot_label()
  test_print_plot_flag()
  test_end_to_end()
  cat("========== Tests Completed ==========\n\n")
  cat("Test Summary:\n")
  print(table(test_results$Result))
  cat("\nDetailed Results:\n")
  print(test_results)
}

# Uncomment to run all tests:
# run_all_tests_anova_bin()
