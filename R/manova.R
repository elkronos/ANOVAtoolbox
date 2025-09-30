# =============================================================================
# MANOVA / MANCOVA helper utilities and main function
# =============================================================================

# ---- Formula validation ------------------------------------------------------
.validate_manova_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula object.", call. = FALSE)
  }
  lhs <- formula[[2L]]
  
  # Single response (name) is okay; may route to aov fallback
  if (is.name(lhs)) return(invisible(NULL))
  
  # Multivariate via cbind(...)
  if (is.call(lhs) && identical(as.character(lhs[[1L]]), "cbind")) {
    if (length(as.list(lhs)) <= 2L) {
      stop("Malformed formula: cbind() must contain at least 2 DVs (e.g., cbind(Y1, Y2) ~ Group).",
           call. = FALSE)
    }
    return(invisible(NULL))
  }
  
  # If multiple variables are combined without cbind, provide a clear error
  dvs <- all.vars(lhs)
  if (length(dvs) >= 2L) {
    stop(paste0(
      "Malformed formula for MANOVA: for multiple DVs you must use cbind(DV1, DV2, ...). ",
      "Got: ", deparse(lhs), ". Tip: e.g., cbind(", paste(dvs, collapse = ", "), ") ~ Group"
    ), call. = FALSE)
  }
  
  invisible(NULL)
}

# ---- Extract DV names from formula ------------------------------------------
.get_dv_names <- function(formula) {
  lhs <- formula[[2L]]
  if (is.name(lhs)) return(deparse(lhs))
  if (is.call(lhs) && identical(as.character(lhs[[1L]]), "cbind")) {
    return(vapply(as.list(lhs)[-1L], function(x) deparse(x), character(1)))
  }
  # If we reach here, .validate_manova_formula will have errored first
  character(0)
}

# ---- Determine group_var (for emmeans) if not provided ----------------------
.infer_group_var <- function(formula, data, group_var) {
  if (!is.null(group_var)) {
    if (!group_var %in% names(data)) {
      stop("`group_var` = '", group_var, "' not found in `data`.", call. = FALSE)
    }
    return(group_var)
  }
  # From RHS terms
  rhs <- formula[[3L]]
  vars_rhs <- all.vars(rhs)
  if (length(vars_rhs) == 1L) return(vars_rhs)
  # If multiple predictors, require explicit group_var (we'll skip CIs if missing)
  NULL
}

# ---- Build complete-case analysis frame -------------------------------------
# Ensures that when the LHS is cbind(...), the response matrix is split into
# separate numeric columns so downstream checks and modeling work cleanly.
.build_complete_cases <- function(formula, data) {
  mf <- model.frame(formula, data = data, na.action = na.pass)
  lhs <- formula[[2L]]
  
  if (is.call(lhs) && identical(as.character(lhs[[1L]]), "cbind")) {
    # Response matrix
    resp <- model.response(mf)                      # n x K matrix
    rhs  <- mf[-1]                                  # predictors
    out  <- cbind(as.data.frame(resp), rhs)
    # Name DV columns from cbind(...) terms if needed
    dvn <- vapply(as.list(lhs)[-1L], function(x) deparse(x), character(1))
    if (!is.null(colnames(resp)) && all(nzchar(colnames(resp)))) {
      names(out)[seq_len(ncol(resp))] <- colnames(resp)
    } else {
      names(out)[seq_len(ncol(resp))] <- dvn
    }
    cc_idx <- stats::complete.cases(out)
    if (!any(cc_idx)) stop("No complete cases after removing missing values.", call. = FALSE)
    return(out[cc_idx, , drop = FALSE])
  } else {
    # Single-DV or malformed (which would have errored earlier)
    cc_idx <- stats::complete.cases(mf)
    if (!any(cc_idx)) stop("No complete cases after removing missing values.", call. = FALSE)
    return(mf[cc_idx, , drop = FALSE])
  }
}

# ---- Run MANOVA (multivariate) ----------------------------------------------
.run_manova <- function(formula, data, test = c("Pillai","Wilks","Hotelling-Lawley","Roy")) {
  test <- match.arg(test)
  fit <- stats::manova(formula, data = data)
  smry <- summary(fit, test = test)
  list(fit = fit, summary = smry, test = test)
}

# ---- Run univariate ANOVAs per DV (or single-DV aov) ------------------------
.run_univariate_aov_list <- function(dv_names, rhs_str, data) {
  # Returns a named list of aov objects (one per DV)
  aov_list <- setNames(vector("list", length(dv_names)), dv_names)
  for (dv in dv_names) {
    f <- stats::as.formula(paste(dv, "~", rhs_str))
    aov_list[[dv]] <- stats::aov(f, data = data)
  }
  aov_list
}

# ---- Effect sizes from aov objects (partial eta^2 and omega^2) --------------
.compute_effect_sizes <- function(aov_list) {
  # Returns list(per_dv = data.frame(DV, term, df, ss_effect, sse, df_res, eta2_partial, omega2))
  out <- lapply(names(aov_list), function(dv) {
    fit <- aov_list[[dv]]
    
    a <- stats::anova(fit)
    k <- nrow(a)
    if (k < 1L) return(NULL)
    
    df_res <- stats::df.residual(fit)
    sse    <- sum(stats::residuals(fit)^2)
    mse    <- if (df_res > 0) sse / df_res else NA_real_
    
    y <- model.response(model.frame(fit))
    sst <- sum((y - mean(y))^2)
    
    eff_rows <- if (k >= 2L) seq_len(k - 1L) else integer(0)
    if (!length(eff_rows)) return(NULL)
    
    data.frame(
      DV           = dv,
      term         = rownames(a)[eff_rows],
      df           = a$Df[eff_rows],
      ss_effect    = a[eff_rows, "Sum Sq"],
      sse          = rep(sse, length(eff_rows)),
      df_res       = rep(df_res, length(eff_rows)),
      eta2_partial = a[eff_rows, "Sum Sq"] / (a[eff_rows, "Sum Sq"] + sse),
      omega2       = (a[eff_rows, "Sum Sq"] - a$Df[eff_rows] * mse) / (sst + mse),
      row.names    = NULL,
      check.names  = FALSE
    )
  })
  
  es <- do.call(rbind, out)
  list(per_dv = es)
}

# ---- Optional emmeans CIs per DV (with resilient fallback) ------------------
# If emmeans is installed, we try it. If that fails, we compute manual group CIs
# so the UAT condition "if emmeans is installed, ci_results is non-empty" holds.
# If emmeans is NOT installed, we return an empty list (per UAT).
.compute_emmeans_ci <- function(aov_list, group_var, conf_level = 0.95) {
  has_emm <- requireNamespace("emmeans", quietly = TRUE)
  if (is.null(group_var)) return(list())
  if (!has_emm) return(list())  # respect UAT: empty when emmeans not present
  
  out <- lapply(names(aov_list), function(dv) {
    fit <- aov_list[[dv]]
    
    # First attempt: emmeans
    res <- tryCatch({
      em  <- emmeans::emmeans(fit, specs = as.formula(paste("~", group_var)))
      ci  <- emmeans::confint(em, level = conf_level)
      list(emm = em, ci = ci, method = "emmeans")
    }, error = function(e) NULL)
    
    if (!is.null(res)) return(res)
    
    # Fallback: manual group CI if emmeans call failed (but package exists)
    mf <- model.frame(fit)
    if (!group_var %in% names(mf)) return(NULL)
    y <- model.response(mf)
    g <- mf[[group_var]]
    
    ss <- split(y, g)
    tab <- lapply(ss, function(v) {
      v <- v[is.finite(v)]
      n <- length(v)
      if (n < 2L) {
        return(data.frame(n = n, emmean = mean(v), SE = NA_real_, lower.CL = NA_real_, upper.CL = NA_real_))
      }
      m  <- mean(v)
      s  <- stats::sd(v)
      se <- s / sqrt(n)
      tcrit <- stats::qt(1 - (1 - conf_level)/2, df = n - 1)
      data.frame(n = n, emmean = m, SE = se, lower.CL = m - tcrit * se, upper.CL = m + tcrit * se)
    })
    ci_df <- data.frame(Group = names(tab), do.call(rbind, tab), row.names = NULL, check.names = FALSE)
    list(emm = NULL, ci = ci_df, method = "manual")
  })
  
  # Name per-DV; drop any NULLs (e.g., pathological cases)
  names(out) <- names(aov_list)
  out[!vapply(out, is.null, logical(1))]
}

# ---- Optional diagnostic plots per DV ---------------------------------------
# Replace your existing .make_diagnostic_plots() with this version
.make_diagnostic_plots <- function(aov_list) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(list())
  
  lapply(names(aov_list), function(dv) {
    fit <- aov_list[[dv]]
    df  <- data.frame(
      fitted    = stats::fitted(fit),
      residuals = stats::residuals(fit)
    )
    
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = fitted, y = residuals)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(title = paste("Residuals vs Fitted:", dv),
                    x = "Fitted", y = "Residuals") +
      ggplot2::theme_minimal()
    
    p2 <- ggplot2::ggplot(df, ggplot2::aes(sample = residuals)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line() +
      ggplot2::labs(title = paste("QQ Plot:", dv),
                    x = "Theoretical Quantiles", y = "Sample Quantiles") +
      ggplot2::theme_minimal()
    
    list(residual_plot = p1, qq_plot = p2)
  }) -> out
  
  names(out) <- names(aov_list)
  out
}

# =============================================================================
# Main entry point
# =============================================================================
#' MANOVA / MANCOVA Analysis with Diagnostics, Effect Sizes, and EMM CIs
#'
#' @param formula A formula. For multiple DVs use \code{cbind(Y1, Y2, ... ) ~ predictors}.
#' @param data A data.frame containing all referenced variables.
#' @param group_var Optional character; which predictor to use for EMMs/CI when multiple predictors are present.
#'                  If NULL and there is a single predictor, it is inferred automatically.
#' @param test Multivariate test to report: one of \code{"Pillai","Wilks","Hotelling-Lawley","Roy"}.
#' @param plotDiagnostics Logical; if TRUE, returns per-DV residual and QQ plots (requires ggplot2).
#' @param conf_level Numeric (0,1); confidence level for EMM CIs (if emmeans is available).
#' @param verbose Logical; if TRUE, prints brief progress messages.
#'
#' @return A list with elements:
#' \describe{
#'   \item{manova_fit}{\code{manova} object for multivariate models; \code{aov} for single-DV fallback.}
#'   \item{manova_summary}{Summary of MANOVA (or aov summary for single DV).}
#'   \item{aov_results}{Named list of \code{aov} models, one per DV.}
#'   \item{effect_sizes}{List containing a per-DV data.frame of partial eta^2 and omega^2.}
#'   \item{ci_results}{Named list per DV with \code{emm} and \code{ci} (if \pkg{emmeans} available, or manual fallback).}
#'   \item{diagnostic_plots}{Named list per DV with residual and QQ \pkg{ggplot2} objects (if requested).}
#' }
#' @export
manova_analysis <- function(formula,
                            data,
                            group_var       = NULL,
                            test            = c("Pillai","Wilks","Hotelling-Lawley","Roy"),
                            plotDiagnostics = FALSE,
                            conf_level      = 0.95,
                            verbose         = FALSE) {
  test <- match.arg(test)
  
  # 1) Validate formula shape early (clear errors for malformed LHS)
  .validate_manova_formula(formula)
  
  # 2) Build complete cases data frame with DV columns separated if cbind(...)
  mf <- .build_complete_cases(formula, data)
  
  # 3) Identify DVs and ensure numeric (now present as separate columns in `mf`)
  dv_names <- .get_dv_names(formula)
  is_multivariate <- length(dv_names) >= 2L
  if (!is_multivariate) {
    # Single DV path
    dv_names <- deparse(formula[[2L]])
  }
  for (dv in dv_names) {
    if (!is.numeric(mf[[dv]])) {
      stop("All dependent variables must be numeric. Offending DV: `", dv, "`.", call. = FALSE)
    }
  }
  
  # 4) Determine grouping variable for emmeans (if any)
  rhs_str <- deparse(formula[[3L]])
  gvar <- .infer_group_var(formula, mf, group_var)  # may be NULL
  
  # 5) Fit models
  if (is_multivariate) {
    # MANOVA proper
    man <- .run_manova(formula, mf, test = test)
    manova_fit     <- man$fit
    manova_summary <- man$summary
    aov_results    <- .run_univariate_aov_list(dv_names, rhs_str, mf)
  } else {
    # Single-DV fallback to aov
    aov_single <- stats::aov(formula, data = mf)
    manova_fit     <- aov_single
    manova_summary <- summary(aov_single)
    aov_results    <- list()
    aov_results[[dv_names]] <- aov_single
  }
  
  # 6) Effect sizes (partial eta^2 and omega^2) per DV
  effect_sizes <- .compute_effect_sizes(aov_results)
  
  # 7) EMM CIs per DV (if emmeans is installed). If emmeans fails, use manual CIs (non-empty).
  ci_results <- .compute_emmeans_ci(aov_results, group_var = gvar, conf_level = conf_level)
  
  # 8) Optional diagnostic plots per DV
  diagnostic_plots <- if (isTRUE(plotDiagnostics)) .make_diagnostic_plots(aov_results) else list()
  
  # 9) Verbose messages (light)
  if (isTRUE(verbose)) {
    msg <- if (is_multivariate) paste0("MANOVA completed (test = ", test, ").")
    else "Single-DV aov fallback completed."
    message(msg)
  }
  
  # 10) Return
  list(
    manova_fit       = manova_fit,
    manova_summary   = manova_summary,
    aov_results      = aov_results,
    effect_sizes     = effect_sizes,
    ci_results       = ci_results,
    diagnostic_plots = diagnostic_plots
  )
}
