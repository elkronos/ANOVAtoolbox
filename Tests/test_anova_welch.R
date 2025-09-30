###############################################################################
# Testing Infrastructure - anova_welch UAT (Welch's One-Way ANOVA)
###############################################################################

# Package management for tests -------------------------------------------------
required_packages <- c("ggplot2","data.table","nortest","effsize")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

# NOTE: This UAT assumes `anova_welch()` is defined in the session.

# If a generator isn't present, define a simple one here so tests are runnable.
if (!exists("gen_between", mode = "function")) {
  #' Generate between-groups data (1 or 2 grouping variables), balanced cells
  #' @param n_per_cell Integer, observations per cell
  #' @param levels1 Character vector of levels for Group1 (default c("A","B","C"))
  #' @param levels2 Optional character vector of levels for Group2
  #' @param mu1 Numeric effect for each level in levels1 (same length)
  #' @param mu2 Optional numeric effect for each level in levels2 (same length)
  #' @param base Numeric grand mean
  #' @param sd Residual SD
  gen_between <- function(n_per_cell = 20,
                          levels1 = c("A","B","C"),
                          levels2 = NULL,
                          mu1 = c(0, 0.6, 1.2),
                          mu2 = NULL,
                          base = 10,
                          sd = 1) {
    stopifnot(length(levels1) == length(mu1))
    if (!is.null(levels2)) stopifnot(length(levels2) == length(mu2))
    if (is.null(levels2)) {
      Group1 <- factor(rep(levels1, each = n_per_cell), levels = levels1)
      eff1   <- mu1[as.integer(Group1)]
      DV     <- base + eff1 + rnorm(length(Group1), 0, sd)
      data.frame(DV = DV, Group1 = Group1, stringsAsFactors = FALSE)
    } else {
      grid <- expand.grid(Group1 = levels1, Group2 = levels2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      n_cells <- nrow(grid)
      Group1 <- factor(rep(grid$Group1, each = n_per_cell), levels = levels1)
      Group2 <- factor(rep(grid$Group2, each = n_per_cell), levels = levels2)
      eff1   <- mu1[as.integer(Group1)]
      eff2   <- mu2[as.integer(Group2)]
      DV     <- base + eff1 + eff2 + rnorm(n_cells * n_per_cell, 0, sd)
      data.frame(DV = DV, Group1 = Group1, Group2 = Group2, stringsAsFactors = FALSE)
    }
  }
}

# Global container for test results -------------------------------------------
test_results <- data.frame(Test = character(), Result = character(), stringsAsFactors = FALSE)

# Helper: Print and Store Test Result -----------------------------------------
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
# Unit Tests
###############################################################################

# Test: input validation & basic errors ---------------------------------------
test_input_validation_welch <- function() {
  # Non-data.frame -> error
  res1 <- tryCatch(anova_welch(list(a = 1), "y", "g"), error = identity)
  ok1  <- inherits(res1, "error")
  
  # Missing response var -> error
  df2  <- data.frame(g = factor(rep(letters[1:3], each = 5)))
  res2 <- tryCatch(anova_welch(df2, "y", "g"), error = identity)
  ok2  <- inherits(res2, "error")
  
  # Missing group var -> error
  df3  <- data.frame(y = rnorm(10))
  res3 <- tryCatch(anova_welch(df3, "y", "g"), error = identity)
  ok3  <- inherits(res3, "error")
  
  # Non-numeric response -> error (adjusted code requires numeric)
  df4 <- data.frame(y = as.character(rnorm(10)), g = factor(rep(c("A","B"), each = 5)))
  res4 <- tryCatch(anova_welch(df4, "y", "g"), error = identity)
  ok4  <- inherits(res4, "error")
  
  # Not enough complete rows -> error
  df5 <- data.frame(y = rep(NA_real_, 6), g = factor(rep(c("A","B","C"), each = 2)))
  res5 <- tryCatch(anova_welch(df5, "y", "g"), error = identity)
  ok5  <- inherits(res5, "error")
  
  print_and_store_result("input_validation: non-data.frame errors", ok1)
  print_and_store_result("input_validation: missing response variable errors", ok2)
  print_and_store_result("input_validation: missing group variable errors", ok3)
  print_and_store_result("input_validation: non-numeric response errors", ok4)
  print_and_store_result("input_validation: insufficient complete rows errors", ok5)
}

# Test: NA handling & model runs ----------------------------------------------
test_na_handling_welch <- function() {
  set.seed(101)
  df <- gen_between(n_per_cell = 25, levels1 = c("A","B","C"))
  df$DV[sample.int(nrow(df), 12)] <- NA
  res <- tryCatch(anova_welch(df, "DV", "Group1"), error = identity)
  ok  <- is.list(res) && inherits(res$welch_anova, "htest")
  print_and_store_result("na_handling: model runs and returns htest", ok)
}

# Test: two-group design -> no posthoc ----------------------------------------
test_two_group_no_posthoc <- function() {
  set.seed(202)
  df <- gen_between(n_per_cell = 30, levels1 = c("Ctl","Trt"), mu1 = c(0, 0.8))
  res <- anova_welch(df, "DV", "Group1")
  ok  <- is.null(res$posthoc)
  print_and_store_result("two_groups: posthoc is NULL when only two groups", ok)
}

# Test: multi-group posthoc present & class -----------------------------------
test_posthoc_presence_class <- function() {
  set.seed(303)
  df <- gen_between(n_per_cell = 20, levels1 = c("A","B","C"), mu1 = c(0, 0.6, 1.2))
  res <- anova_welch(df, "DV", "Group1", posthoc_adjust = "holm")
  t1  <- is.list(res) && inherits(res$welch_anova, "htest")
  t2  <- !is.null(res$posthoc) && inherits(res$posthoc, "pairwise.htest")
  t3  <- isTRUE(tryCatch(res$posthoc$p.adjust.method == "holm", error = function(e) FALSE))
  print_and_store_result("posthoc: pairwise Welch t-tests returned (pairwise.htest)", t2 && t1)
  print_and_store_result("posthoc: p.adjust method set to holm", t3)
}

# Test: multiple grouping variables supported ---------------------------------
test_multiple_group_vars <- function() {
  set.seed(404)
  df <- gen_between(n_per_cell = 15,
                    levels1 = c("A","B"),
                    levels2 = c("X","Y"),
                    mu1 = c(0, 1.0),
                    mu2 = c(0, 0.5))
  res <- anova_welch(df, "DV", c("Group1","Group2"))
  # Expect 4 cells => posthoc present
  t1 <- length(unique(interaction(df$Group1, df$Group2))) == 4
  t2 <- is.list(res) && inherits(res$welch_anova, "htest")
  t3 <- !is.null(res$posthoc) && inherits(res$posthoc, "pairwise.htest")
  print_and_store_result("multi_group_vars: runs with two grouping variables", t1 && t2 && t3)
}

# Test: effect sizes presence/shape (tolerate missing effsize pkg) ------------
test_effect_sizes_presence <- function() {
  set.seed(505)
  df <- gen_between(n_per_cell = 20, levels1 = c("A","B","C"))
  res <- anova_welch(df, "DV", "Group1")
  ok <- is.null(res$effect_sizes) ||
    (is.data.frame(res$effect_sizes) &&
       all(c("group1","group2","cohen_d","magnitude") %in% names(res$effect_sizes)))
  note <- if (is.null(res$effect_sizes)) "effsize not installed or no pairs computed" else NULL
  print_and_store_result("effect_sizes: NULL (allowed) or data.table with expected columns", ok, note)
}

# Test: normality diagnostics (AD test behavior) ------------------------------
test_ad_test_behavior <- function() {
  # Small sample: residuals < 8 => NA
  set.seed(606)
  df_small <- gen_between(n_per_cell = 2, levels1 = c("A","B","C"))  # n=6, k=3, residuals ~ 3
  res_small <- anova_welch(df_small, "DV", "Group1")
  t_small <- isTRUE(is.na(res_small$assumptions$ad_test))
  print_and_store_result("normality: AD test skipped (residuals < 8 => NA)", t_small)
  
  # Larger sample: if nortest available, expect htest; else NA
  set.seed(607)
  df_big <- gen_between(n_per_cell = 12, levels1 = c("A","B","C"))   # n=36
  res_big <- anova_welch(df_big, "DV", "Group1")
  ad_is_ok <- isTRUE(is.na(res_big$assumptions$ad_test)) ||
    inherits(res_big$assumptions$ad_test, "htest")
  note <- if (isTRUE(is.na(res_big$assumptions$ad_test))) "nortest not installed" else NULL
  print_and_store_result("normality: AD test returns htest when available (or NA if skipped)", ad_is_ok, note)
}

# Test: summary stats columns & conf_level effect ------------------------------
test_summary_stats_and_conflevel <- function() {
  set.seed(808)
  df <- gen_between(n_per_cell = 40, levels1 = c("A","B","C"))
  r95 <- anova_welch(df, "DV", "Group1", conf_level = 0.95)
  r90 <- anova_welch(df, "DV", "Group1", conf_level = 0.90)
  
  cols_ok <- all(c("group_interaction","mean","sd","n","se","ci_low","ci_high") %in%
                   names(r95$summary_stats))
  w95 <- with(r95$summary_stats, ci_high - ci_low)
  w90 <- with(r90$summary_stats, ci_high - ci_low)
  width_ok <- mean(w90, na.rm = TRUE) < mean(w95, na.rm = TRUE)
  
  print_and_store_result("summary_stats: expected columns present", cols_ok)
  print_and_store_result("conf_level: 90% CI narrower than 95%", width_ok)
}

# Test: means plot returned (ggplot) -------------------------------------------
test_means_plot_is_ggplot <- function() {
  set.seed(909)
  df <- gen_between(n_per_cell = 20, levels1 = c("A","B","C"))
  res <- anova_welch(df, "DV", "Group1")
  ok  <- inherits(res$means_plot, "ggplot")
  print_and_store_result("means_plot: returns ggplot object", ok)
}

# Test: input immutability (no column leakage into original data) -------------
test_no_input_mutation <- function() {
  set.seed(1001)
  df <- gen_between(n_per_cell = 15, levels1 = c("A","B","C"))
  df_copy <- df
  invisible(anova_welch(df_copy, "DV", "Group1"))
  ok <- !("group_interaction" %in% names(df))
  print_and_store_result("immutability: original data not mutated with helper columns", ok)
}

# Test: one-level grouping should error ---------------------------------------
test_one_level_error <- function() {
  df <- gen_between(n_per_cell = 20, levels1 = c("A"), mu1 = 0)
  res <- tryCatch(anova_welch(df, "DV", "Group1"), error = identity)
  ok  <- inherits(res, "error")
  print_and_store_result("grouping: error when only one distinct group", ok)
}

# Test: end-to-end scenario ----------------------------------------------------
test_end_to_end_welch <- function() {
  set.seed(1111)
  df <- gen_between(n_per_cell = 25,
                    levels1 = c("A","B"),
                    levels2 = c("X","Y"),
                    mu1 = c(0, 0.8),
                    mu2 = c(0, 0.3))
  res <- anova_welch(df, "DV", c("Group1","Group2"), conf_level = 0.95, posthoc_adjust = "bonferroni")
  needed <- c("assumptions","welch_anova","posthoc","summary_stats","effect_sizes","means_plot")
  t1 <- is.list(res) && all(needed %in% names(res))
  t2 <- inherits(res$welch_anova, "htest")
  t3 <- is.null(res$effect_sizes) || (is.data.frame(res$effect_sizes) &&
                                        all(c("group1","group2","cohen_d","magnitude") %in% names(res$effect_sizes)))
  t4 <- inherits(res$means_plot, "ggplot")
  print_and_store_result("end_to_end: full suite of outputs present & typed", t1 && t2 && t3 && t4)
}

###############################################################################
# Test Runner
###############################################################################

run_all_tests_anova_welch <- function() {
  cat("========== Running Comprehensive anova_welch UAT ==========\n")
  
  test_input_validation_welch()
  test_na_handling_welch()
  test_two_group_no_posthoc()
  test_posthoc_presence_class()
  test_multiple_group_vars()
  test_effect_sizes_presence()
  test_ad_test_behavior()
  test_summary_stats_and_conflevel()
  test_means_plot_is_ggplot()
  test_no_input_mutation()
  test_one_level_error()
  test_end_to_end_welch()
  
  cat("========== Tests Completed ==========\n\n")
  cat("Test Summary:\n"); print(table(test_results$Result))
  cat("\nDetailed Results:\n"); print(test_results)
}

# Uncomment to run all tests:
# run_all_tests_anova_welch()
