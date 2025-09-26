###############################################################################
# Testing Infrastructure - anova_kw UAT (Kruskal–Wallis + Dunn post-hoc)
###############################################################################

# Package management for tests -------------------------------------------------
required_packages <- c(
  "ggplot2","data.table","nortest","dunn.test"
)
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(nortest)
  library(dunn.test)
})

# NOTE: This UAT assumes `anova_kw()` is already defined in the session.

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
# Synthetic Data Generators (KW-friendly) -- robust for any n
###############################################################################

# One-factor numeric response with known location shifts across >=3 levels
gen_kw_one_factor <- function(n = 450, levels = c("A","B","C"), means = c(0, 0.7, 1.2), sd = 1) {
  stopifnot(length(levels) == length(means))
  g <- factor(rep(levels, length.out = n), levels = levels)
  mu <- setNames(means, levels)[as.character(g)]
  y <- rnorm(n, mu, sd)
  data.frame(y = y, g = g, stringsAsFactors = FALSE)
}

# Two-factor numeric response with combined-location effects
gen_kw_two_factor <- function(n = 600,
                              g1_levels = c("A","B"),
                              g2_levels = c("X","Y"),
                              cell_means = matrix(c(0.0, 0.6, 1.0, 1.4), nrow = 2, byrow = TRUE),
                              sd = 1.0) {
  stopifnot(nrow(cell_means) == length(g1_levels), ncol(cell_means) == length(g2_levels))
  cells <- expand.grid(g1 = g1_levels, g2 = g2_levels, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  idx <- rep(seq_len(nrow(cells)), length.out = n)
  g1 <- factor(cells$g1[idx], levels = g1_levels)
  g2 <- factor(cells$g2[idx], levels = g2_levels)
  mu <- ifelse(g1==g1_levels[1] & g2==g2_levels[1], cell_means[1,1],
               ifelse(g1==g1_levels[1] & g2==g2_levels[2], cell_means[1,2],
                      ifelse(g1==g1_levels[2] & g2==g2_levels[1], cell_means[2,1], cell_means[2,2])))
  y <- rnorm(n, mu, sd)
  data.frame(y = y, g1 = g1, g2 = g2, stringsAsFactors = FALSE)
}

# Heavy-tailed generator to stress-test non-normality diagnostics
gen_kw_heavy_tailed <- function(n = 480, levels = c("A","B","C"), loc = c(0, 0.8, 1.4), df = 3) {
  g <- factor(rep(levels, length.out = n), levels = levels)
  loc_vec <- setNames(loc, levels)[as.character(g)]
  # t distribution centered at 'loc'
  y <- loc_vec + stats::rt(n, df = df)
  data.frame(y = y, g = g, stringsAsFactors = FALSE)
}

###############################################################################
# Unit Tests
###############################################################################

# Test: input validation & factor coercion ------------------------------------
test_input_validation_kw <- function() {
  df <- gen_kw_one_factor(200); names(df)[1] <- "resp"
  e1 <- inherits(tryCatch(anova_kw(df, "missing_resp", "g"), error = identity), "error")
  e2 <- inherits(tryCatch(anova_kw(df, "resp", "missing_g"), error = identity), "error")
  
  # Non-numeric response should error
  df_bad <- df; df_bad$resp <- as.character(df_bad$resp)
  e3 <- inherits(tryCatch(anova_kw(df_bad, "resp", "g"), error = identity), "error")
  
  # Character group coerced to factor inside function (verify in returned data_used)
  df_chr <- df; df_chr$g <- as.character(df_chr$g)
  res <- anova_kw(df_chr, "resp", "g", plot_qq = FALSE)
  coerced <- is.factor(res$data_used$g)
  
  print_and_store_result("input_validation: errors when response missing", e1)
  print_and_store_result("input_validation: errors when group var missing", e2)
  print_and_store_result("input_validation: error when response non-numeric", e3)
  print_and_store_result("input_validation: character group coerced to factor (data_used)", coerced)
}

# Test: NA handling (rows dropped & notes reflect it) -------------------------
test_na_handling_kw <- function() {
  set.seed(1)
  df <- gen_kw_two_factor(360)
  idx_y <- sample.int(nrow(df), 20)
  idx_g <- sample.int(nrow(df), 15)
  df$y[idx_y]  <- NA
  df$g1[idx_g] <- NA
  
  res <- anova_kw(df, "y", c("g1","g2"), plot_qq = FALSE)
  expected_drop <- sum(stats::complete.cases(df[, c("y","g1","g2")]) == FALSE)
  t1 <- identical(res$n_dropped, expected_drop)
  t2 <- any(grepl("Dropped", paste(res$notes, collapse = " ")))
  t3 <- all(!is.na(res$kruskal_test$statistic)) && is.numeric(res$kruskal_test$statistic)
  
  print_and_store_result("na_handling: n_dropped equals expected", t1)
  print_and_store_result("na_handling: notes mention rows dropped", t2)
  print_and_store_result("na_handling: kruskal_test returns statistic", t3)
}

# Test: insufficient groups -> error ------------------------------------------
test_insufficient_groups_error <- function() {
  df <- gen_kw_one_factor(150, levels = c("A"), means = c(0))
  e <- inherits(tryCatch(anova_kw(df, "y", "g"), error = identity), "error")
  print_and_store_result("insufficient_groups: error when only 1 group", e)
}

# Test: KW + Dunn structure and sizes (one factor, 3+ groups) -----------------
test_kw_posthoc_structure <- function() {
  set.seed(2)
  df <- gen_kw_one_factor(600, levels = c("A","B","C","D"), means = c(0, 0.5, 1.0, 1.5))
  res <- anova_kw(df, "y", "g", plot_qq = FALSE, posthoc_method = "bh")
  
  # KW test object sanity
  t_kw <- inherits(res$kruskal_test, "htest") && res$kruskal_test$parameter >= 1
  
  # Posthoc tidy frame present with required columns and expected number of comparisons
  tidy_ok <- !is.null(res$posthoc) && is.list(res$posthoc) &&
    all(c("raw","tidy") %in% names(res$posthoc)) &&
    is.data.frame(res$posthoc$tidy) &&
    all(c("comparison","Z","p_unadjusted","p_adjusted","group1","group2") %in% names(res$posthoc$tidy))
  k <- nlevels(df$g)
  n_expected <- choose(k, 2)
  n_rows_ok <- nrow(res$posthoc$tidy) == n_expected
  padj_bounds <- all(is.finite(res$posthoc$tidy$p_adjusted)) &&
    all(res$posthoc$tidy$p_adjusted >= 0 & res$posthoc$tidy$p_adjusted <= 1)
  
  print_and_store_result("kw_posthoc: kruskal.test returns valid htest object", t_kw)
  print_and_store_result("kw_posthoc: Dunn tidy table present with required columns", tidy_ok)
  print_and_store_result("kw_posthoc: Dunn number of comparisons is choose(k,2)", n_rows_ok)
  print_and_store_result("kw_posthoc: adjusted p-values in [0,1]", padj_bounds)
}

# Test: post-hoc adjustment method variants -----------------------------------
test_posthoc_methods <- function() {
  df <- gen_kw_one_factor(480)
  res_bh  <- anova_kw(df, "y", "g", plot_qq = FALSE, posthoc_method = "bh")
  res_bon <- anova_kw(df, "y", "g", plot_qq = FALSE, posthoc_method = "bonferroni")
  
  ok_bh  <- all(res_bh$posthoc$tidy$p_adjusted >= 0 & res_bh$posthoc$tidy$p_adjusted <= 1)
  ok_bon <- all(res_bon$posthoc$tidy$p_adjusted >= 0 & res_bon$posthoc$tidy$p_adjusted <= 1)
  
  # Usually Bonferroni is more conservative: p_adj_bon >= p_adj_bh for many rows
  comp_vec <- res_bon$posthoc$tidy$p_adjusted >= res_bh$posthoc$tidy$p_adjusted
  conservative_often <- mean(comp_vec) > 0.7  # allow noise
  
  print_and_store_result("posthoc_methods: BH adjusted p in [0,1]", ok_bh)
  print_and_store_result("posthoc_methods: Bonferroni adjusted p in [0,1]", ok_bon)
  print_and_store_result("posthoc_methods: Bonferroni typically >= BH", conservative_often)
}

# Test: QQ-plot toggle and AD test presence -----------------------------------
test_qqplot_toggle <- function() {
  df <- gen_kw_heavy_tailed(420)
  res_off <- anova_kw(df, "y", "g", plot_qq = FALSE)
  res_on  <- anova_kw(df, "y", "g", plot_qq = TRUE)
  
  t1 <- is.null(res_off$assumptions$qqplot)
  t2 <- inherits(res_on$assumptions$qqplot, "ggplot")
  t3 <- inherits(res_off$assumptions$ad_test, "htest") &&
    inherits(res_on$assumptions$ad_test, "htest")
  
  print_and_store_result("qqplot: absent when disabled", t1)
  print_and_store_result("qqplot: ggplot object when enabled", t2)
  print_and_store_result("qqplot: AD test objects always present", t3)
}

# Test: boxplot object & x categories / faceting ------------------------------
test_boxplot_structure <- function() {
  df <- gen_kw_two_factor(500)
  res <- anova_kw(df, "y", c("g1","g2"), plot_qq = FALSE)
  
  t1 <- is.list(res$plots) && inherits(res$plots$boxplot, "ggplot")
  # Extract number of x categories from built plot
  gb <- ggplot_build(res$plots$boxplot)
  # Heuristic: look at first non-empty layer
  layer_idx <- which(lengths(gb$data) > 0)[1]
  n_x_found <- length(unique(round(gb$data[[layer_idx]]$x)))
  n_x_expected <- nlevels(df$g1)  # x-axis is first group by design
  
  # Faceting presence: two-factor input -> facetting by second factor bundle
  facet_ok <- !is.null(res$data_used$facet_var)
  print_and_store_result("boxplot: returns ggplot object", t1)
  print_and_store_result("boxplot: x categories match first grouping levels", n_x_found == n_x_expected)
  print_and_store_result("boxplot: facet_var present for 2+ grouping variables", facet_ok)
}

# Test: reorder by median yields increasing order of medians on x -------------
test_boxplot_reorder_by_median <- function() {
  set.seed(3)
  # Construct groups with clearly different medians
  df <- gen_kw_one_factor(450, levels = c("L1","L2","L3"),
                          means = c(1.5, 0.0, 0.8), sd = 0.4)
  res <- anova_kw(df, "y", "g", plot_qq = FALSE)
  p <- res$plots$boxplot
  # Reconstruct x-level order from the built object
  gb <- ggplot_build(p)
  data_layer <- gb$data[[which(lengths(gb$data) > 0)[1]]]
  # x positions map to the order of the discrete scale; recover labels via scale mapping
  dlabs <- ggplot_build(p)$layout$panel_params[[1]]$x$get_labels()
  # Strip counts "(n=...)" if present
  dlabs_clean <- gsub("\\s*\\(n=\\d+\\)$", "", dlabs)
  # Compute medians by original factor and compare in plotted order
  med <- tapply(df$y, df$g, median, na.rm = TRUE)
  med_in_order <- med[dlabs_clean]
  nondecreasing <- all(diff(as.numeric(med_in_order)) >= -1e-8)
  
  print_and_store_result("boxplot_reorder: x levels ordered by increasing median", nondecreasing)
}

# Test: no printing side-effects (capture console output) ---------------------
test_no_print_side_effects <- function() {
  # dunn.test() prints to console; capture output to ensure function is quiet when captured.
  df <- gen_kw_one_factor(200)
  res <- NULL
  captured <- utils::capture.output({
    res <- anova_kw(df, "y", "g", plot_qq = FALSE)
  })
  # Assert we got a proper result object and nothing escaped our capture
  ok <- is.list(res) && all(c("data_used","kruskal_test","posthoc","plots") %in% names(res))
  print_and_store_result("side_effects: function runs cleanly with console capture", ok)
}

# End-to-end scenarios ---------------------------------------------------------
test_end_to_end_kw <- function() {
  # Scenario 1: One factor, clear differences, BH post-hoc
  set.seed(10)
  df1 <- gen_kw_one_factor(600, levels = c("A","B","C"), means = c(0, 0.9, 1.6))
  res1 <- anova_kw(df1, "y", "g", plot_qq = TRUE, posthoc_method = "bh")
  need1 <- c("data_used","assumptions","kruskal_test","posthoc","plots","n_dropped","notes")
  t1 <- is.list(res1) && all(need1 %in% names(res1)) && inherits(res1$plots$boxplot, "ggplot")
  
  # Scenario 2: Two factors, faceted boxplot, Bonferroni post-hoc
  df2 <- gen_kw_two_factor(700)
  res2 <- anova_kw(df2, "y", c("g1","g2"), plot_qq = FALSE, posthoc_method = "bonferroni")
  t2 <- is.list(res2$posthoc) && is.data.frame(res2$posthoc$tidy) &&
    inherits(res2$plots$boxplot, "ggplot") && !is.null(res2$data_used$facet_var)
  
  # Scenario 3: Heavy-tailed data; AD likely non-normal
  df3 <- gen_kw_heavy_tailed(500)
  res3 <- anova_kw(df3, "y", "g", plot_qq = TRUE)
  ad_p <- res3$assumptions$ad_test$p.value
  t3 <- is.finite(ad_p) && ad_p >= 0 && ad_p <= 1
  
  print_and_store_result("end_to_end: one-factor BH + qqplot returns full suite", t1)
  print_and_store_result("end_to_end: two-factor faceting + Bonferroni OK", t2)
  print_and_store_result("end_to_end: heavy-tailed diagnostic AD p-value numeric", t3)
}

###############################################################################
# Test Runner
###############################################################################

#' Run All anova_kw Tests
#'
#' Execute a comprehensive suite of tests for the `anova_kw()` function.
#'
#' @details
#' All tests write a PASS/FAIL line and populate a global \code{test_results} table.
run_all_tests_anova_kw <- function() {
  cat("========== Running Comprehensive anova_kw UAT ==========\n")
  test_input_validation_kw()
  test_na_handling_kw()
  test_insufficient_groups_error()
  test_kw_posthoc_structure()
  test_posthoc_methods()
  test_qqplot_toggle()
  test_boxplot_structure()
  test_boxplot_reorder_by_median()
  test_no_print_side_effects()
  test_end_to_end_kw()
  cat("========== Tests Completed ==========\n\n")
  cat("Test Summary:\n")
  print(table(test_results$Result))
  cat("\nDetailed Results:\n")
  print(test_results)
}

# Uncomment to run all tests:
# run_all_tests_anova_kw()
