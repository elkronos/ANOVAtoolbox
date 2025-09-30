# =========================
# UAT: Full Test Suite Using testthat for manova_analysis()
# =========================

# ---- Lightweight synthetic data generator for MANOVA/MANCOVA (no extra deps) ----
.gen_manova_data <- function(n = 120,
                             groups = c("A","B","C"),
                             dv_names = c("DV1","DV2"),
                             covar = FALSE,
                             group_effect = NULL,   # auto-fits to length(groups) if NULL/length 1
                             sigma = NULL,         # auto PD matrix if NULL
                             unbalanced = FALSE,
                             seed = 42) {
  set.seed(seed)
  G <- length(groups)
  K <- length(dv_names)
  
  # group membership (balanced/unbalanced)
  if (unbalanced) {
    probs <- runif(G); probs <- probs / sum(probs)
    g <- factor(sample(groups, n, TRUE, prob = probs), levels = groups)
  } else {
    g <- factor(rep(groups, length.out = n), levels = groups)
  }
  
  # group_effect handling (robust)
  if (is.null(group_effect)) {
    # gentle monotone shift across groups
    group_effect <- seq(0, 0.8, length.out = G)
  } else if (length(group_effect) == 1L) {
    group_effect <- rep(group_effect, G)
  } else if (length(group_effect) != G) {
    stop("`group_effect` must be length 1 or length(groups).")
  }
  
  # covariance matrix for DVs
  if (is.null(sigma)) {
    # moderate positive correlations
    sigma <- matrix(0.3, nrow = K, ncol = K); diag(sigma) <- 1
  } else {
    if (!is.matrix(sigma) || any(dim(sigma) != K)) {
      stop("`sigma` must be a K x K matrix where K = length(dv_names).")
    }
  }
  
  # multivariate normal errors using base R (no MASS)
  Z <- matrix(rnorm(n * K), nrow = n, ncol = K)
  L <- tryCatch(chol(sigma), error = function(e) diag(K))
  E <- Z %*% t(L)
  
  # group mean shift applied to all DVs (plus tiny across-DV differentiation)
  shift_g <- group_effect[as.integer(g)]
  shift_mat <- matrix(shift_g, nrow = n, ncol = K)
  shift_mat <- sweep(shift_mat, 2, seq(0, 0.3, length.out = K), `+`)
  
  Y <- shift_mat + E
  colnames(Y) <- dv_names
  
  df <- data.frame(Y, Group = g, check.names = FALSE)
  if (isTRUE(covar)) {
    df$Covar <- scale(rnorm(n) + 0.2 * as.numeric(g))[, 1]
  }
  df
}

# =========================
# Test Runner (call this)
# =========================
run_all_tests_manova <- function() {
  if (!requireNamespace("testthat", quietly = TRUE)) {
    message("Package 'testthat' not installed; skipping UAT tests.")
    return(invisible(NULL))
  }
  
  message("\n==================== Running UAT Tests for manova_analysis() ====================")
  library(testthat)
  
  # Helper to quiet function chatter (plots/messages/warnings)
  .quiet <- function(expr) {
    suppressWarnings(suppressMessages(capture.output(result <- force(expr))))
    result
  }
  
  # ---- Test 1: Basic functionality with iris (2 DVs, single grouping) ----
  test_that("manova_analysis works on iris dataset", {
    data(iris)
    res <- .quiet(manova_analysis(cbind(Sepal.Length, Sepal.Width) ~ Species, data = iris))
    
    expect_s3_class(res$manova_fit, "manova")
    expect_true(!is.null(res$manova_summary))
    expect_true(is.list(res$aov_results))
    expect_true(is.list(res$effect_sizes))
    if (requireNamespace("emmeans", quietly = TRUE)) {
      expect_true(length(res$ci_results) > 0)
    }
  })
  
  # ---- Test 2: Explicit group_var when multiple predictors exist ----------
  test_that("manova_analysis handles explicit group_var when multiple predictors exist", {
    set.seed(123)
    dat <- data.frame(
      DV1 = rnorm(100),
      DV2 = rnorm(100),
      Group = sample(letters[1:3], 100, replace = TRUE),
      Covar = rnorm(100)
    )
    res <- .quiet(manova_analysis(cbind(DV1, DV2) ~ Group + Covar, data = dat, group_var = "Group"))
    expect_s3_class(res$manova_fit, "manova")
  })
  
  # ---- Test 3: Errors when group_var missing ------------------------------
  test_that("manova_analysis errors when group_var is not in data", {
    dat <- iris
    expect_error(
      manova_analysis(cbind(Sepal.Length, Sepal.Width) ~ Species, data = dat, group_var = "Nonexistent")
    )
  })
  
  # ---- Test 4: Diagnostic plots generated when requested -------------------
  test_that("Diagnostic plots are generated when requested", {
    data(iris)
    res <- .quiet(manova_analysis(cbind(Sepal.Length, Sepal.Width) ~ Species,
                                  data = iris, plotDiagnostics = TRUE))
    expect_true(is.list(res$diagnostic_plots))
    for (dv in names(res$diagnostic_plots)) {
      expect_true(all(c("residual_plot", "qq_plot") %in% names(res$diagnostic_plots[[dv]])))
    }
  })
  
  # ---- Test 5: Diagnostic plots skipped when disabled ----------------------
  test_that("Diagnostic plots are skipped when plotDiagnostics = FALSE", {
    data(iris)
    res <- .quiet(manova_analysis(cbind(Sepal.Length, Sepal.Width) ~ Species,
                                  data = iris, plotDiagnostics = FALSE))
    expect_equal(length(res$diagnostic_plots), 0)
  })
  
  # ---- Test 6: NA handling across DVs -------------------------------------
  test_that("manova_analysis handles NA rows by omission (complete-cases)", {
    set.seed(77)
    dat <- .gen_manova_data(n = 90, dv_names = c("DV1","DV2"), groups = c("A","B","C"))
    idx <- sample.int(nrow(dat), 10)
    dat$DV1[idx[1:5]] <- NA
    dat$DV2[idx[6:10]] <- NA
    res <- .quiet(manova_analysis(cbind(DV1, DV2) ~ Group, data = dat))
    expect_s3_class(res$manova_fit, "manova")
    expect_true(is.list(res$aov_results))
  })
  
  # ---- Test 7: Non-numeric DV should error ---------------------------------
  test_that("manova_analysis errors when any DV is non-numeric", {
    dat <- .gen_manova_data(n = 50, dv_names = c("DV1","DV2"), groups = c("A","B"))
    dat$DV2 <- as.character(dat$DV2)
    expect_error(manova_analysis(cbind(DV1, DV2) ~ Group, data = dat))
  })
  
  # ---- Test 8: Unbalanced groups ------------------------------------------
  test_that("manova_analysis runs with unbalanced groups", {
    set.seed(101)
    dat <- .gen_manova_data(n = 150, dv_names = c("DV1","DV2"),
                            groups = c("A","B","C"), unbalanced = TRUE)
    res <- .quiet(manova_analysis(cbind(DV1, DV2) ~ Group, data = dat))
    expect_s3_class(res$manova_fit, "manova")
  })
  
  # ---- Test 9: More than two DVs ------------------------------------------
  test_that("manova_analysis supports 3+ DVs", {
    set.seed(202)
    dat <- .gen_manova_data(n = 120, dv_names = c("Y1","Y2","Y3"), groups = c("A","B","C"))
    res <- .quiet(manova_analysis(cbind(Y1, Y2, Y3) ~ Group, data = dat))
    expect_s3_class(res$manova_fit, "manova")
    expect_true(is.list(res$aov_results))
  })
  
  # ---- Test 10: MANCOVA with covariate ------------------------------------
  test_that("manova_analysis supports covariates (MANCOVA)", {
    set.seed(303)
    dat <- .gen_manova_data(n = 140, dv_names = c("DV1","DV2"),
                            groups = c("Ctl","Trt"), covar = TRUE)
    res <- .quiet(manova_analysis(cbind(DV1, DV2) ~ Group + Covar, data = dat, group_var = "Group"))
    expect_s3_class(res$manova_fit, "manova")
    expect_true(!is.null(res$manova_summary))
  })
  
  # ---- Test 11: Effect sizes container present -----------------------------
  test_that("effect_sizes is a list (may contain multivariate and per-DV entries)", {
    dat <- .gen_manova_data(n = 120, dv_names = c("A","B"), groups = c("G1","G2","G3"))
    res <- .quiet(manova_analysis(cbind(A, B) ~ Group, data = dat))
    expect_true(is.list(res$effect_sizes))
  })
  
  # ---- Test 12: CI results behavior w/ emmeans ------------------------------
  test_that("CI results follow emmeans availability", {
    dat <- .gen_manova_data(n = 100, dv_names = c("A","B"), groups = c("G1","G2","G3"))
    res <- .quiet(manova_analysis(cbind(A, B) ~ Group, data = dat))
    if (requireNamespace("emmeans", quietly = TRUE)) {
      expect_true(length(res$ci_results) > 0)
    } else {
      expect_true(length(res$ci_results) == 0 || is.null(res$ci_results))
    }
  })
  
  # ---- Test 13: Single-DV behavior (error OR valid fallback both OK) -------
  test_that("Single-DV behavior is acceptable (clear error OR valid fallback)", {
    dat <- .gen_manova_data(n = 90, dv_names = c("OnlyY"), groups = c("A","B"))
    got <- tryCatch(
      .quiet(manova_analysis(OnlyY ~ Group, data = dat)),
      error = function(e) e
    )
    if (inherits(got, "error")) {
      # If implementation errors, ensure message is clear
      expect_match(conditionMessage(got), "need multiple responses|at least.*two|single", ignore.case = TRUE)
    } else {
      # If implementation falls back to aov(), that is also acceptable
      expect_true(
        (is.list(got) && !is.null(got$manova_fit)) ||
          inherits(got$manova_fit, "manova") ||
          inherits(got$manova_fit, "aov")
      )
    }
  })
  
  # ---- Test 14: Clear error for malformed formula ---------------------------
  test_that("malformed formula produces a clear error", {
    dat <- .gen_manova_data(n = 60, dv_names = c("Y1","Y2"), groups = c("A","B"))
    expect_error(
      manova_analysis(Y1 + Y2 ~ Group, data = dat),
      regexp = "cbind|DVs|formula|need multiple responses",  # allow legacy/base error, too
      ignore.case = TRUE
    )
  })
  
  # ---- Test 15: End-to-end, rich scenario -----------------------------------
  test_that("End-to-end scenario returns the full suite", {
    set.seed(909)
    dat <- .gen_manova_data(n = 160, dv_names = c("DV1","DV2","DV3"),
                            groups = c("A","B","C"), covar = TRUE, unbalanced = TRUE)
    res <- .quiet(manova_analysis(cbind(DV1, DV2, DV3) ~ Group + Covar,
                                  data = dat, group_var = "Group", plotDiagnostics = TRUE))
    needed <- c("manova_fit","manova_summary","aov_results","effect_sizes","ci_results","diagnostic_plots")
    expect_true(is.list(res) && all(needed %in% names(res)))
    expect_true(is.list(res$diagnostic_plots))
  })
  
  message("All UAT tests run completed.")
  invisible(TRUE)
}

# Example (uncomment to run):
# run_all_tests_manova()
