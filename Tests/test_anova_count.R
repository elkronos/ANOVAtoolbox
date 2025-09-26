###############################################################################
# Testing Infrastructure - anova_count UAT (Poisson/NB factorial "ANOVA")
###############################################################################

# Package management for tests -------------------------------------------------
required_packages <- c(
  "MASS","car","emmeans","multcomp","sandwich","ggplot2","dplyr","broom",
  "withr","rlang","scales","magrittr"
)
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
suppressPackageStartupMessages({
  library(MASS)
  library(car)
  library(emmeans)
  library(multcomp)
  library(sandwich)
  library(ggplot2)
  library(dplyr)
  library(broom)
  library(withr)
  library(rlang)
  library(scales)
  library(magrittr)
})

# NOTE: This UAT assumes `anova_count()` is already defined in the session.

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
# Synthetic Data Generators
###############################################################################

# Poisson data: 1–2 factors; optional interaction effect on log-rate
gen_pois_data <- function(n = 400, include_g2 = TRUE,
                          rates = matrix(c(5,10,15,20), nrow = 2, byrow = TRUE)) {
  if (include_g2) {
    g1 <- factor(rep(c("A","B"), each = n/2))
    g2 <- factor(rep(rep(c("X","Y"), each = n/4), 2))
    lam <- ifelse(g1=="A" & g2=="X", rates[1,1],
                  ifelse(g1=="A" & g2=="Y", rates[1,2],
                         ifelse(g1=="B" & g2=="X", rates[2,1], rates[2,2])))
    count <- rpois(n, lam)
    data.frame(count = count, g1 = g1, g2 = g2)
  } else {
    g1 <- factor(rep(c("A","B"), each = n/2))
    lam <- ifelse(g1=="A", rates[1,1], rates[2,2])
    count <- rpois(n, lam)
    data.frame(count = count, g1 = g1)
  }
}

# Negative Binomial (overdispersed) data with strong dispersion
gen_nb_data <- function(n = 400, theta = 0.7) {
  g1 <- factor(rep(c("A","B"), each = n/2))
  g2 <- factor(rep(rep(c("X","Y"), each = n/4), 2))
  lam <- ifelse(g1=="A" & g2=="X", 6,
                ifelse(g1=="A" & g2=="Y", 12,
                       ifelse(g1=="B" & g2=="X", 18, 24)))
  count <- rnbinom(n, size = theta, mu = lam)
  data.frame(count = count, g1 = g1, g2 = g2)
}

# Offset data: exposure/time per row with constant rate across rows within cell
gen_offset_data <- function(n = 400, base_rate = 2.5) {
  g1 <- factor(rep(c("A","B"), each = n/2))
  g2 <- factor(rep(rep(c("X","Y"), each = n/4), 2))
  exposure <- runif(n, min = 0.5, max = 5) # strictly > 0
  mult <- ifelse(g1=="A" & g2=="X", 1.0,
                 ifelse(g1=="A" & g2=="Y", 1.4,
                        ifelse(g1=="B" & g2=="X", 1.8, 2.2)))
  lambda <- base_rate * mult * exposure
  count <- rpois(n, lambda)
  data.frame(count = count, g1 = g1, g2 = g2, exposure = exposure)
}

###############################################################################
# Unit Tests
###############################################################################

# Test: input validation -------------------------------------------------------
test_input_validation <- function() {
  df <- gen_pois_data(200, include_g2 = FALSE)
  names(df)[names(df) == "count"] <- "y"
  # Missing response
  e1 <- inherits(tryCatch(anova_count(df, "y_missing", "g1", plot = FALSE), error = identity), "error")
  # Missing group var
  e2 <- inherits(tryCatch(anova_count(df, "y", "g_missing", plot = FALSE), error = identity), "error")
  # Non-integer (force doubles)
  df_bad <- df; df_bad$y <- df_bad$y + 0.5
  e3 <- inherits(tryCatch(anova_count(df_bad, "y", "g1", plot = FALSE), error = identity), "error")
  # Negative counts
  df_bad2 <- df; df_bad2$y[1] <- -1
  e4 <- inherits(tryCatch(anova_count(df_bad2, "y", "g1", plot = FALSE), error = identity), "error")
  # Group with < 2 levels
  df_bad3 <- df; df_bad3$g1 <- factor("A")
  e5 <- inherits(tryCatch(anova_count(df_bad3, "y", "g1", plot = FALSE), error = identity), "error")
  # Offset invalid (nonpositive)
  df_off <- gen_offset_data(100); df_off$exposure[1] <- 0
  e6 <- inherits(tryCatch(anova_count(df_off, "count", c("g1","g2"), offset_var = "exposure", plot = FALSE), error = identity), "error")
  # Unknown offset column
  e7 <- inherits(tryCatch(anova_count(df_off, "count", c("g1","g2"), offset_var = "nope", plot = FALSE), error = identity), "error")
  
  print_and_store_result("input_validation: errors when response missing", e1)
  print_and_store_result("input_validation: errors when group var missing", e2)
  print_and_store_result("input_validation: errors for non-integer counts", e3)
  print_and_store_result("input_validation: errors for negative counts", e4)
  print_and_store_result("input_validation: errors when group has < 2 levels", e5)
  print_and_store_result("input_validation: errors when offset <= 0", e6)
  print_and_store_result("input_validation: errors when offset column missing", e7)
}

# Test: NA handling ------------------------------------------------------------
test_na_handling <- function() {
  df <- gen_pois_data(300, include_g2 = TRUE)
  idx_y <- sample.int(nrow(df), 15); idx_g <- sample.int(nrow(df), 11)
  df$count[idx_y] <- NA
  df$g1[idx_g] <- NA
  
  res <- anova_count(df, "count", c("g1","g2"), plot = FALSE)
  used_n <- length(res$model$y)
  expect_n <- sum(stats::complete.cases(df[, c("count","g1","g2")]))
  t1 <- identical(used_n, expect_n)
  t2 <- is.numeric(res$n_removed_na) && res$n_removed_na == (nrow(df) - expect_n)
  
  print_and_store_result("na_handling: model uses only complete cases", t1)
  print_and_store_result("na_handling: n_removed_na equals dropped rows", t2)
}

# Test: basic Poisson fit, emmeans, and plot ----------------------------------
test_basic_pois <- function() {
  df <- gen_pois_data(400, include_g2 = TRUE)
  res <- anova_count(df, "count", c("g1","g2"), plot = FALSE)
  
  need <- c("model","model_type","overdispersion_statistic","anova_table",
            "emm_grid","posthoc_pairs","emmeans_table","plot","messages")
  t1 <- is.list(res) && all(need %in% names(res)) && inherits(res$model, "glm")
  t2 <- identical(res$model_type, "poisson")
  t3 <- !is.null(res$anova_table)
  
  # emmeans_table normalized to response/lower.CL/upper.CL
  has_cols <- is.data.frame(res$emmeans_table) &&
    all(c("response","lower.CL","upper.CL") %in% names(res$emmeans_table))
  
  # posthoc_pairs structure guaranteed; may be empty but columns must exist
  has_pairs <- is.data.frame(res$posthoc_pairs) &&
    all(c("contrast","IRR","lower.CL","upper.CL","p.value") %in% names(res$posthoc_pairs))
  # if non-empty, IRRs must be finite
  pairs_finite <- nrow(res$posthoc_pairs) == 0 || all(is.finite(res$posthoc_pairs$IRR))
  
  print_and_store_result("pois_fit: returns all key objects", t1)
  print_and_store_result("pois_fit: model_type is 'poisson'", t2)
  print_and_store_result("pois_fit: anova_table is present", t3)
  print_and_store_result("pois_fit: emmeans_table has response & CIs", has_cols)
  print_and_store_result("pois_fit: posthoc_pairs has IRR & CIs (finite if present)", has_pairs && pairs_finite)
}

# Test: overdispersion detection and NB fallback -------------------------------
test_overdispersion_nb <- function() {
  set.seed(123)
  df <- gen_nb_data(500, theta = 0.6)
  res <- anova_count(df, "count", c("g1","g2"),
                     overdispersion_threshold = 1.2,
                     use_nb_if_overdispersed = TRUE,
                     plot = FALSE)
  is_nb <- identical(res$model_type, "negbin")
  msg_has_od <- any(grepl("Overdispersion detected", res$messages))
  t1 <- isTRUE(is_nb) || isTRUE(msg_has_od)
  t2 <- !is.null(res$anova_table)
  
  print_and_store_result("overdispersion: flagged and/or NB fallback engaged", t1)
  print_and_store_result("overdispersion: anova_table present under NB/OD", t2)
}

# Test: robust vcov pathway ----------------------------------------------------
test_robust_vcov <- function() {
  df <- gen_pois_data(400, include_g2 = TRUE)
  res <- anova_count(df, "count", c("g1","g2"), vcov_type = "robust", plot = FALSE)
  
  t1 <- is.data.frame(res$posthoc_pairs) &&
    all(c("contrast","IRR","lower.CL","upper.CL","p.value") %in% names(res$posthoc_pairs))
  t2 <- nrow(res$posthoc_pairs) == 0 || all(is.finite(res$posthoc_pairs$IRR))
  t3 <- isTRUE(is.character(res$messages))
  
  print_and_store_result("robust_vcov: posthoc_pairs returned with IRRs", t1 && t2)
  print_and_store_result("robust_vcov: messages vector returned", t3)
}

# Test: ANOVA test statistic option (LR/Wald) ----------------------------------
test_anova_modes <- function() {
  df <- gen_pois_data(400, include_g2 = TRUE)
  res_lr <- anova_count(df, "count", c("g1","g2"), type_anova = "LR", plot = FALSE)
  res_wald <- anova_count(df, "count", c("g1","g2"), type_anova = "Wald", plot = FALSE)
  
  t1 <- !is.null(res_lr$anova_table)
  t2 <- !is.null(res_wald$anova_table)
  print_and_store_result("anova_modes: LR returns an anova_table", t1)
  print_and_store_result("anova_modes: Wald returns an anova_table", t2)
}

# Test: offset handling --------------------------------------------------------
test_offset <- function() {
  df <- gen_offset_data(500, base_rate = 3.0)
  res <- anova_count(df, "count", c("g1","g2"), offset_var = "exposure", plot = FALSE)
  
  t1 <- is.data.frame(res$emmeans_table) &&
    all(c("response","lower.CL","upper.CL") %in% names(res$emmeans_table)) &&
    all(res$emmeans_table$response >= 0)
  t2 <- is.data.frame(res$posthoc_pairs) &&
    all(c("contrast","IRR","lower.CL","upper.CL","p.value") %in% names(res$posthoc_pairs)) &&
    (nrow(res$posthoc_pairs) == 0 || all(is.finite(res$posthoc_pairs$IRR)))
  
  print_and_store_result("offset: emmeans_table present with nonnegative rates", t1)
  print_and_store_result("offset: posthoc_pairs contains finite IRRs", t2)
}

# Test: sparse / tiny cells warning -------------------------------------------
test_sparse_cells_warning <- function() {
  df <- gen_pois_data(200, include_g2 = TRUE)
  # Deliberately reduce some cells
  keep_idx <- c(1:5, 6:100, 101:110, 111:200)
  df_small <- df[keep_idx, ]
  res <- anova_count(df_small, "count", c("g1","g2"), plot = FALSE)
  msg <- paste(res$messages, collapse = " | ")
  t1 <- any(grepl("small counts", msg)) || any(grepl("zero observations", msg))
  print_and_store_result("sparse_cells: messages indicate small/zero cells", t1)
}

# Test: plot object structure --------------------------------------------------
test_plot_object <- function() {
  df <- gen_pois_data(320, include_g2 = TRUE)
  res <- anova_count(df, "count", c("g1","g2"), plot = TRUE)
  t1 <- inherits(res$plot, "ggplot") || inherits(res$plot, "gg")
  # pick first non-empty layer to count x positions
  gb <- ggplot_build(res$plot)
  layer_idx <- which(lengths(gb$data) > 0)[1]
  n_found <- length(unique(gb$data[[layer_idx]]$x))
  n_expected <- nlevels(df$g1) * nlevels(df$g2)
  t2 <- n_found == n_expected
  
  print_and_store_result("plot: returns ggplot object", t1)
  print_and_store_result("plot: x categories match group combinations", t2)
}

# Test: messages and bookkeeping ----------------------------------------------
test_messages_and_flags <- function() {
  df <- gen_nb_data(500, theta = 0.8)
  res <- anova_count(df, "count", c("g1","g2"),
                     overdispersion_threshold = 1.1,
                     use_nb_if_overdispersed = FALSE,
                     plot = FALSE)
  t1 <- isTRUE(res$overdispersion_flagged)
  t2 <- identical(res$model_type, "poisson")
  t3 <- any(grepl("Overdispersion detected", res$messages))
  t4 <- is.numeric(res$overdispersion_statistic) && is.finite(res$overdispersion_statistic)
  
  print_and_store_result("messages: overdispersion_flagged is TRUE", t1)
  print_and_store_result("messages: model_type remains 'poisson' when NB disabled", t2)
  print_and_store_result("messages: overdispersion note present", t3)
  print_and_store_result("messages: dispersion statistic is finite", t4)
}

# Test: normalization of emmeans column names ---------------------------------
test_emmeans_name_normalization <- function() {
  # This confirms that any `rate/asymp.LCL/asymp.UCL` are mapped to response/lower.CL/upper.CL
  df <- gen_pois_data(300, include_g2 = TRUE)
  res <- anova_count(df, "count", c("g1","g2"), plot = FALSE)
  nm <- names(res$emmeans_table)
  t1 <- all(c("response","lower.CL","upper.CL") %in% nm)
  print_and_store_result("emmeans: names normalized to response/lower.CL/upper.CL", t1)
}

# End-to-end scenarios ---------------------------------------------------------
test_end_to_end <- function() {
  # Scenario 1: One factor Poisson, robust vcov
  df1 <- gen_pois_data(500, include_g2 = FALSE)
  res1 <- anova_count(df1, "count", "g1", vcov_type = "robust", plot = FALSE)
  need1 <- c("model","anova_table","emmeans_table","posthoc_pairs","messages")
  t1 <- is.list(res1) && all(need1 %in% names(res1)) && inherits(res1$model, "glm") &&
    is.data.frame(res1$emmeans_table) && is.data.frame(res1$posthoc_pairs)
  
  # Scenario 2: Two factors, NB fallback, LR test
  df2 <- gen_nb_data(600, theta = 0.6)
  res2 <- anova_count(df2, "count", c("g1","g2"),
                      overdispersion_threshold = 1.2,
                      use_nb_if_overdispersed = TRUE,
                      type_anova = "LR",
                      plot = FALSE)
  t2 <- (identical(res2$model_type, "negbin") || any(grepl("Overdispersion detected", res2$messages))) &&
    !is.null(res2$anova_table)
  
  # Scenario 3: Two factors with offset, Wald test, robust vcov
  df3 <- gen_offset_data(700, base_rate = 2.0)
  res3 <- anova_count(df3, "count", c("g1","g2"),
                      offset_var = "exposure",
                      type_anova = "Wald",
                      vcov_type = "robust",
                      plot = TRUE)
  t3 <- is.data.frame(res3$emmeans_table) && is.data.frame(res3$posthoc_pairs) &&
    (inherits(res3$plot, "ggplot") || inherits(res3$plot, "gg"))
  
  print_and_store_result("end_to_end: Scenario 1 returns core objects (robust Poisson)", t1)
  print_and_store_result("end_to_end: Scenario 2 NB fallback + LR anova", t2)
  print_and_store_result("end_to_end: Scenario 3 offset + Wald + robust + plot", t3)
}

###############################################################################
# Test Runner
###############################################################################

#' Run All anova_count Tests
#'
#' Execute a comprehensive suite of tests for the `anova_count()` function.
#'
#' @details
#' All tests write a PASS/FAIL line and populate a global \code{test_results} table.
run_all_tests_anova_count <- function() {
  cat("========== Running Comprehensive anova_count UAT ==========\n")
  test_input_validation()
  test_na_handling()
  test_basic_pois()
  test_overdispersion_nb()
  test_robust_vcov()
  test_anova_modes()
  test_offset()
  test_sparse_cells_warning()
  test_plot_object()
  test_messages_and_flags()
  test_emmeans_name_normalization()
  test_end_to_end()
  cat("========== Tests Completed ==========\n\n")
  cat("Test Summary:\n")
  print(table(test_results$Result))
  cat("\nDetailed Results:\n")
  print(test_results)
}

# Uncomment to run all tests:
# run_all_tests_anova_count()
