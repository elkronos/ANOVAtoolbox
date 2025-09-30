###############################################################################
# Testing Infrastructure - anova_rm UAT (Repeated Measures ANOVA)
###############################################################################

# Package management for tests -------------------------------------------------
required_packages <- c(
  "afex","emmeans","car","effectsize","ggplot2"
)
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
suppressPackageStartupMessages({
  library(afex)
  library(emmeans)
  library(car)
  library(effectsize)
  library(ggplot2)
})

# NOTE: This UAT assumes `anova_rm()` is defined in the session.
# If a generator isn't present, define a simple one here so tests are runnable.
if (!exists("gen_rm_one_within", mode = "function")) {
  #' Generate long-format repeated-measures data (one within factor, optional between)
  #' @param n Integer subjects
  #' @param within_levels Character vector of within levels (default c("T1","T2","T3"))
  #' @param means_by_within Numeric means per within level (same length as within_levels)
  #' @param between_levels Optional character vector levels for between factor "Group"
  #' @param between_shift Numeric shift added to the last between level relative to first
  #' @param sd Residual SD
  gen_rm_one_within <- function(
    n,
    within_levels   = c("T1","T2","T3"),
    means_by_within = c(0, 0.6, 1.2),
    between_levels  = NULL,
    between_shift   = 0.8,
    sd = 1
  ) {
    stopifnot(length(within_levels) == length(means_by_within))
    k <- length(within_levels)
    ID   <- factor(rep(seq_len(n), each = k))
    Time <- factor(rep(within_levels, times = n), levels = within_levels)
    
    # Optional between-subject factor
    if (!is.null(between_levels)) {
      grp_by_id <- factor(sample(between_levels, n, replace = TRUE), levels = between_levels)
      Group <- factor(rep(grp_by_id, each = k), levels = between_levels)
      grp_shift_vec <- if (length(between_levels) >= 2) {
        if (length(between_levels) == 2) c(0, between_shift) else seq(0, between_shift, length.out = length(between_levels))
      } else 0
      group_mean <- grp_shift_vec[as.integer(Group)]
    } else {
      Group <- NULL
      group_mean <- 0
    }
    
    # Subject random intercept + mild subject trend across Time to induce correlation
    subj_intercept <- rnorm(n, 0, 0.6)
    subj_slope     <- rnorm(n, 0, 0.2)
    subj_i <- as.integer(ID)
    time_num <- as.numeric(Time) - 1
    
    mu_time <- means_by_within[as.integer(Time)]
    DV <- subj_intercept[subj_i] + subj_slope[subj_i] * time_num + mu_time + group_mean +
      rnorm(n * k, 0, sd)
    
    out <- data.frame(ID = ID, DV = DV, Time = Time, stringsAsFactors = FALSE)
    if (!is.null(Group)) out$Group <- Group
    out
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

# Test: input validation & coercion -------------------------------------------
test_input_validation_rm <- function() {
  df <- gen_rm_one_within(20)
  
  # Pass a non-existent subject column name -> should error
  res1 <- tryCatch(anova_rm(df, "Subject", "DV", "Time", verbose = FALSE), error = identity)
  ok1  <- inherits(res1, "error")
  
  # Pass a non-existent DV column name -> should error
  res2 <- tryCatch(anova_rm(df, "ID", "Response", "Time", verbose = FALSE), error = identity)
  ok2  <- inherits(res2, "error")
  
  # DV coerced from character
  df_chr <- df; df_chr$DV <- as.character(df_chr$DV)
  res3 <- tryCatch(anova_rm(df_chr, "ID", "DV", "Time", verbose = FALSE), error = identity)
  ok3  <- is.list(res3)
  
  print_and_store_result("input_validation: error when required columns missing (subject)", ok1)
  print_and_store_result("input_validation: error when required columns missing (dv)", ok2)
  print_and_store_result("input_validation: DV character coerced (result structure present)", ok3)
}

# Test: NA handling ------------------------------------------------------------
test_na_handling_rm <- function() {
  set.seed(2)
  df <- gen_rm_one_within(30)
  df$DV[sample.int(nrow(df), 10)] <- NA
  res <- tryCatch(anova_rm(df, "ID", "DV", "Time", verbose = FALSE), error = identity)
  ok <- is.list(res) && inherits(res$anova_obj, "afex_aov")
  print_and_store_result("na_handling: model runs and returns afex_aov", ok)
}

# Test: two-level within design ------------------------------------------------
test_two_level_within <- function() {
  df <- gen_rm_one_within(40,
                          within_levels   = c("T1","T2"),
                          means_by_within = c(0, 0.5))
  res <- tryCatch(anova_rm(df, "ID", "DV", "Time", verbose = FALSE), error = identity)
  ok <- is.list(res) && inherits(res$anova_obj, "afex_aov")
  print_and_store_result("two_level_within: model runs cleanly with 2 levels", ok)
}

# Test: effect sizes present ---------------------------------------------------
test_effect_sizes_rm <- function() {
  df <- gen_rm_one_within(30)
  res <- anova_rm(df, "ID", "DV", "Time", verbose = FALSE)
  ok <- !is.null(res$effect_sizes) && is.data.frame(res$effect_sizes)
  print_and_store_result("effect_sizes: eta_squared table returned", ok)
}

# Test: EMM plot creation ------------------------------------------------------
test_emm_plot_rm <- function() {
  df <- gen_rm_one_within(30)
  res <- anova_rm(df, "ID", "DV", "Time", emm_specs = "Time", verbose = FALSE)
  ok <- is.null(res$emm_plot) || inherits(res$emm_plot, "ggplot")
  print_and_store_result("emm_plot: returns NULL or ggplot object", ok)
}

# Test: normality diagnostics --------------------------------------------------
test_normality_rm <- function() {
  df <- gen_rm_one_within(30)
  res <- anova_rm(df, "ID", "DV", "Time", verbose = FALSE)
  t1 <- !is.null(res$normality$residuals)
  t2 <- is.null(res$normality$shapiro_test) || inherits(res$normality$shapiro_test, "htest")
  print_and_store_result("normality: residuals present", t1)
  print_and_store_result("normality: shapiro.test returns htest (or NULL if residuals unavailable)", t2)
}

# Test: sphericity diagnostics (Mauchly) --------------------------------------
test_sphericity_rm <- function() {
  df <- gen_rm_one_within(60,
                          within_levels   = c("T1","T2","T3"),
                          means_by_within = c(0, 0.5, 1.0))
  res <- anova_rm(df, "ID", "DV", "Time", verbose = FALSE)
  # Returns result of car::Mauchly.test(...) or NULL if not applicable
  ok <- is.null(res$mauchly) || inherits(res$mauchly, "anova") || is.data.frame(res$mauchly) || is.list(res$mauchly)
  print_and_store_result("sphericity: NULL (not applicable) or Mauchly test object returned", ok)
}

# Test: between-subject effect sizes ------------------------------------------
test_between_effect_sizes <- function() {
  df <- gen_rm_one_within(40,
                          within_levels   = c("T1","T2"),
                          means_by_within = c(0, 1),
                          between_levels  = c("Ctl","Trt"),
                          between_shift   = 1.0)
  res <- anova_rm(df, "ID", "DV", "Time", between = "Group", verbose = FALSE)
  es_txt <- paste(capture.output(print(res$effect_sizes)), collapse = "\n")
  ok <- grepl("Group", es_txt, fixed = TRUE)
  print_and_store_result("between: effect_sizes include Group term", ok)
}

# Test: diagnostics toggle -----------------------------------------------------
test_diagnostics_toggle <- function() {
  df <- gen_rm_one_within(25, within_levels = c("T1","T2","T3"))
  res_off <- anova_rm(df, "ID", "DV", "Time", diagnostics = FALSE, verbose = FALSE)
  res_on  <- anova_rm(df, "ID", "DV", "Time", diagnostics = TRUE,  verbose = FALSE)
  t1 <- is.null(res_off$diagnostic_plots)
  t2 <- is.null(res_on$diagnostic_plots) || inherits(res_on$diagnostic_plots, "recordedplot")
  print_and_store_result("diagnostics_toggle: diagnostic_plots NULL when disabled", t1)
  print_and_store_result("diagnostics_toggle: diagnostic_plots returned when available", t2)
}

# Test: verbose toggle ---------------------------------------------------------
test_verbose_toggle <- function() {
  df <- gen_rm_one_within(20)
  captured <- utils::capture.output({
    res <- anova_rm(df, "ID", "DV", "Time", verbose = FALSE)
  })
  ok <- exists("res") && is.list(res) && length(captured) == 0
  print_and_store_result("verbose_toggle: no console output when verbose=FALSE", ok)
}

# Test: end-to-end scenario ----------------------------------------------------
test_end_to_end_rm <- function() {
  set.seed(10)
  df <- gen_rm_one_within(50,
                          within_levels   = c("T1","T2","T3"),
                          means_by_within = c(0, 0.6, 1.2))
  res <- anova_rm(df, "ID", "DV", "Time", emm_specs = "Time", diagnostics = TRUE, verbose = FALSE)
  needed <- c("anova_obj","effect_sizes","normality","mauchly","diagnostic_plots","emm_plot")
  ok <- is.list(res) && all(needed %in% names(res))
  print_and_store_result("end_to_end: full suite of outputs present", ok)
}

###############################################################################
# Test Runner
###############################################################################

run_all_tests_anova_rm <- function() {
  cat("========== Running Comprehensive anova_rm UAT ==========\n")
  # Ensure sum contrasts to avoid afex contrast chatter
  old_contr <- getOption("contrasts")
  on.exit(options(contrasts = old_contr), add = TRUE)
  options(contrasts = c("contr.sum","contr.poly"))
  
  test_input_validation_rm()
  test_na_handling_rm()
  test_two_level_within()
  test_effect_sizes_rm()
  test_emm_plot_rm()
  test_normality_rm()
  test_sphericity_rm()
  test_between_effect_sizes()
  test_diagnostics_toggle()
  test_verbose_toggle()
  test_end_to_end_rm()
  cat("========== Tests Completed ==========\n\n")
  cat("Test Summary:\n"); print(table(test_results$Result))
  cat("\nDetailed Results:\n"); print(test_results)
}

# Uncomment to run all tests:
# run_all_tests_anova_rm()
