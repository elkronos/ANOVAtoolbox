# =============================
# MANOVA Analysis
# =============================

#' Check Required Packages
#'
#' This function checks whether each package listed in the input vector is installed.
#' If any package is not available, it stops execution and provides an informative error message.
#'
#' @param pkgs A character vector of package names to check.
#'
#' @return No return value; the function is called for its side effect.
#' @export
#'
#' @examples
#' \dontrun{
#'   check_required_packages(c("ggplot2", "dplyr"))
#' }
check_required_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required. Please install it (e.g., install.packages('%s')).", pkg, pkg))
    }
  }
}

#' Determine Grouping Variable
#'
#' This function determines the grouping variable for the MANOVA analysis.
#' If a grouping variable (`group_var`) is provided, it checks whether it exists in the data.
#' If `group_var` is NULL, it attempts to extract a single predictor from the right-hand side of the formula.
#'
#' @param formula A formula object specifying the MANOVA model.
#' @param data A data frame containing the data.
#' @param group_var (Optional) A string specifying the grouping variable. Defaults to NULL.
#'
#' @return A string representing the grouping variable if found, or NULL if not determined.
#' @export
#'
#' @examples
#' \dontrun{
#'   # If the formula has a single predictor:
#'   get_group_variable(cbind(Y1, Y2) ~ Group, data = mydata)
#'
#'   # If multiple predictors are present, specify the grouping variable explicitly:
#'   get_group_variable(cbind(Y1, Y2) ~ Covariate, data = mydata, group_var = "Group")
#' }
get_group_variable <- function(formula, data, group_var) {
  # Check if data is a data.frame
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data.frame.")
  }
  
  # If group_var is provided, verify that it exists in the data.
  if (!is.null(group_var)) {
    if (!group_var %in% names(data)) {
      stop(sprintf("Grouping variable '%s' not found in 'data'.", group_var))
    }
    return(group_var)
  }
  
  # Else, try to extract the single predictor from the formula
  rhs_vars <- all.vars(formula[[3]])
  if (length(rhs_vars) == 1) {
    message("Using grouping variable: ", rhs_vars)
    return(rhs_vars)
  } else {
    warning("Multiple predictors in the formula. Please specify the grouping variable (via 'group_var'). Some assumption tests will be skipped.")
    return(NULL)
  }
}

#' Run Mardia's Multivariate Normality Tests by Group
#'
#' This function applies Mardia's multivariate normality test to each group in the data.
#' It uses the `mvn` function from the MVN package. If the grouping variable is NULL, the test is skipped.
#'
#' @param formula A formula object specifying the MANOVA model.
#' @param data A data frame containing the data.
#' @param group_var A string indicating the grouping variable. If NULL, the test is skipped.
#'
#' @return A list of Mardia's test results for each group, or an empty list if no grouping variable is provided.
#' @importFrom MVN mvn
#' @export
#'
#' @examples
#' \dontrun{
#'   run_normality_tests(cbind(Y1, Y2) ~ Group, data = mydata, group_var = "Group")
#' }
run_normality_tests <- function(formula, data, group_var) {
  normality_results <- list()
  
  if (is.null(group_var)) {
    message("No grouping variable provided; skipping multivariate normality tests.")
    return(normality_results)
  }
  
  dep_vars <- all.vars(formula[[2]])
  groups <- unique(data[[group_var]])
  
  message("\n--- Running Mardia's Multivariate Normality Tests by Group ---")
  for (g in groups) {
    group_data <- data[data[[group_var]] == g, , drop = FALSE]
    dep_data <- group_data[, dep_vars, drop = FALSE]
    
    # Run Mardia's test (removed verbose argument)
    mardia_out <- MVN::mvn(
      data = dep_data,
      mvnTest = "mardia",
      multivariatePlot = "none"
    )
    normality_results[[as.character(g)]] <- mardia_out$multivariateNormality
    message(sprintf("Group %s: Mardia's test p-value = %.4f", g, mardia_out$multivariateNormality$p))
  }
  
  return(normality_results)
}

#' Run Box's M Test for Homogeneity of Covariance Matrices
#'
#' This function performs Box's M test for the homogeneity of covariance matrices
#' using the `boxM` function from the heplots package. If no grouping variable is provided,
#' the test is skipped.
#'
#' @param formula A formula object specifying the MANOVA model.
#' @param data A data frame containing the data.
#' @param group_var A string specifying the grouping variable. If NULL, the test is skipped.
#'
#' @return An object containing the results of Box's M test, or NULL if the test is skipped.
#' @importFrom heplots boxM
#' @export
#'
#' @examples
#' \dontrun{
#'   run_boxM_test(cbind(Y1, Y2) ~ Group, data = mydata, group_var = "Group")
#' }
run_boxM_test <- function(formula, data, group_var) {
  if (is.null(group_var)) {
    message("No grouping variable provided; skipping Box's M test.")
    return(NULL)
  }
  
  dep_vars <- all.vars(formula[[2]])
  dep_data <- data[, dep_vars, drop = FALSE]
  grouping_factor <- data[[group_var]]
  
  boxM_result <- heplots::boxM(dep_data, grouping_factor)
  message(sprintf("\n--- Box's M Test ---\nBox's M test p-value = %.4f", boxM_result$p.value))
  return(boxM_result)
}

#' Run MANOVA Analysis
#'
#' This function performs a multivariate analysis of variance (MANOVA) using the provided
#' formula and data, and returns both the fitted model and its summary using Pillai's Trace.
#'
#' @param formula A formula object specifying the MANOVA model.
#' @param data A data frame containing the data.
#'
#' @return A list with two components:
#'   \item{fit}{The MANOVA model object.}
#'   \item{summary}{The summary of the MANOVA model using Pillai's Trace.}
#' @export
#'
#' @examples
#' \dontrun{
#'   res <- run_manova(cbind(Y1, Y2) ~ Group, data = mydata)
#'   print(res$summary)
#' }
run_manova <- function(formula, data) {
  message("\n--- Running MANOVA ---")
  manova_fit <- manova(formula, data = data)
  manova_summary <- summary(manova_fit, test = "Pillai")
  message("\nMANOVA Summary (using Pillai's Trace):")
  print(manova_summary)
  list(fit = manova_fit, summary = manova_summary)
}

#' Run Univariate Follow-up Analyses (ANOVAs)
#'
#' This function conducts univariate analysis of variance (ANOVA) follow-up tests on a MANOVA model.
#'
#' @param manova_fit A MANOVA model object returned from `manova()`.
#'
#' @return A list of ANOVA tables for each dependent variable.
#' @export
#'
#' @examples
#' \dontrun{
#'   manova_out <- run_manova(cbind(Y1, Y2) ~ Group, data = mydata)
#'   aov_results <- run_univariate_analyses(manova_out$fit)
#' }
run_univariate_analyses <- function(manova_fit) {
  aov_results <- summary.aov(manova_fit)
  return(aov_results)
}

#' Compute Partial Eta Squared Effect Sizes
#'
#' This function computes partial eta squared effect sizes from the ANOVA tables obtained
#' via univariate follow-up analyses.
#'
#' @param aov_results A list of ANOVA tables, as returned by `summary.aov()`.
#'
#' @return A list where each element is a data frame containing the effect sizes for the corresponding dependent variable.
#' @export
#'
#' @examples
#' \dontrun{
#'   aov_results <- run_univariate_analyses(manova_fit)
#'   effect_sizes <- compute_effect_sizes(aov_results)
#' }
compute_effect_sizes <- function(aov_results) {
  effect_sizes <- list()
  message("\n--- Computing Partial Eta Squared for Each Dependent Variable ---")
  
  for (dv in names(aov_results)) {
    aov_table <- aov_results[[dv]]
    # Ensure that there is at least one effect (other than Residuals)
    if (nrow(aov_table) < 2) {
      warning(sprintf("ANOVA table for '%s' does not contain enough rows to compute effect sizes.", dv))
      next
    }
    
    effects <- rownames(aov_table)[-nrow(aov_table)]  # remove Residuals
    es_df <- data.frame(Effect = effects,
                        SS_effect = NA_real_,
                        SS_error = NA_real_,
                        partial_eta2 = NA_real_,
                        stringsAsFactors = FALSE)
    
    for (i in seq_along(effects)) {
      effect_name <- effects[i]
      SS_effect <- aov_table[effect_name, "Sum Sq"]
      SS_error  <- aov_table["Residuals", "Sum Sq"]
      pe2 <- SS_effect / (SS_effect + SS_error)
      es_df[i, "SS_effect"] <- SS_effect
      es_df[i, "SS_error"] <- SS_error
      es_df[i, "partial_eta2"] <- pe2
    }
    effect_sizes[[dv]] <- es_df
    message(sprintf("\nEffect sizes for dependent variable '%s':", dv))
    print(es_df)
  }
  
  return(effect_sizes)
}

#' Compute Estimated Marginal Means and Confidence Intervals
#'
#' This function computes estimated marginal means and their confidence intervals for each
#' dependent variable using linear models. It utilizes the `emmeans` package if installed.
#'
#' @param formula A formula object specifying the MANOVA model.
#' @param data A data frame containing the data.
#' @param group_var A string specifying the grouping variable.
#'
#' @return A list of estimated marginal means and confidence intervals for each dependent variable.
#' @importFrom emmeans emmeans
#' @export
#'
#' @examples
#' \dontrun{
#'   ci_results <- compute_ci_estimates(cbind(Y1, Y2) ~ Group, data = mydata, group_var = "Group")
#' }
compute_ci_estimates <- function(formula, data, group_var) {
  ci_results <- list()
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    warning("Package 'emmeans' not installed. Skipping estimated marginal means and confidence intervals.")
    return(ci_results)
  }
  
  dep_vars <- all.vars(formula[[2]])
  predictors <- all.vars(formula[[3]])
  
  message("\n--- Computing Estimated Marginal Means and Confidence Intervals ---")
  for (dv in dep_vars) {
    lm_formula <- as.formula(paste(dv, "~", paste(predictors, collapse = " + ")))
    lm_fit <- lm(lm_formula, data = data)
    emm <- emmeans::emmeans(lm_fit, specs = group_var)
    ci <- summary(emm)
    ci_results[[dv]] <- ci
    message(sprintf("\nEMMs and CIs for '%s':", dv))
    print(ci)
  }
  
  return(ci_results)
}

#' Generate Diagnostic Plots for Each Dependent Variable
#'
#' This function generates diagnostic plots (Residuals vs Fitted and QQ plots) for each dependent variable
#' using linear models and ggplot2. If `plotDiagnostics` is FALSE, the function returns an empty list.
#'
#' @param formula A formula object specifying the MANOVA model.
#' @param data A data frame containing the data.
#' @param plotDiagnostics Logical indicating whether to generate diagnostic plots. Defaults to TRUE.
#'
#' @return A list where each element is a list containing ggplot objects for the diagnostic plots of a dependent variable.
#' @importFrom ggplot2 ggplot aes geom_point geom_hline labs theme_minimal stat_qq stat_qq_line
#' @export
#'
#' @examples
#' \dontrun{
#'   diagnostic_plots <- generate_diagnostic_plots(cbind(Y1, Y2) ~ Group, data = mydata, plotDiagnostics = TRUE)
#' }
generate_diagnostic_plots <- function(formula, data, plotDiagnostics) {
  diagnostic_plots <- list()
  if (!plotDiagnostics) {
    message("plotDiagnostics is FALSE; skipping diagnostic plots.")
    return(diagnostic_plots)
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' not installed. Cannot generate diagnostic plots.")
    return(diagnostic_plots)
  }
  
  dep_vars <- all.vars(formula[[2]])
  predictors <- all.vars(formula[[3]])
  
  message("\n--- Generating Diagnostic Plots ---")
  for (dv in dep_vars) {
    lm_formula <- as.formula(paste(dv, "~", paste(predictors, collapse = " + ")))
    lm_fit <- lm(lm_formula, data = data)
    
    # Create a data frame with fitted values and residuals
    lm_df <- data.frame(
      fitted = lm_fit$fitted.values,
      resid = lm_fit$residuals,
      stdresid = rstandard(lm_fit)
    )
    
    # Residuals vs Fitted Plot
    p1 <- ggplot2::ggplot(lm_df, ggplot2::aes(x = fitted, y = resid)) +
      ggplot2::geom_point(color = "blue", size = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(title = paste("Residuals vs Fitted for", dv),
                    x = "Fitted values", y = "Residuals") +
      ggplot2::theme_minimal()
    
    # QQ Plot of standardized residuals
    p2 <- ggplot2::ggplot(lm_df, ggplot2::aes(sample = stdresid)) +
      ggplot2::stat_qq(color = "darkgreen") +
      ggplot2::stat_qq_line(color = "red") +
      ggplot2::labs(title = paste("QQ Plot for", dv)) +
      ggplot2::theme_minimal()
    
    diagnostic_plots[[dv]] <- list(residual_plot = p1, qq_plot = p2)
    
    # Print the plots to the current graphics device
    print(p1)
    print(p2)
  }
  
  return(diagnostic_plots)
}

#' Generate Canonical Discriminant Analysis Plot
#'
#' This function generates a canonical discriminant analysis plot for a MANOVA model
#' using the `candisc` package, if available.
#'
#' @param manova_fit A MANOVA model object.
#'
#' @return A plot object resulting from the canonical discriminant analysis, or NULL if the `candisc` package is not installed.
#' @importFrom candisc candisc
#' @export
#'
#' @examples
#' \dontrun{
#'   candisc_plot <- generate_candisc_plot(manova_fit)
#' }
generate_candisc_plot <- function(manova_fit) {
  candisc_plot <- NULL
  if (!requireNamespace("candisc", quietly = TRUE)) {
    message("Package 'candisc' not installed; skipping canonical discriminant analysis plot.")
    return(candisc_plot)
  }
  
  message("\n--- Generating Canonical Discriminant Analysis Plot ---")
  candisc_fit <- candisc::candisc(manova_fit)
  # The plot() function for candisc objects opens a graphics window.
  candisc_plot <- plot(candisc_fit)
  return(candisc_plot)
}

#' Run a Robust MANOVA Analysis with Assumption Tests, Effect Sizes, Confidence Intervals, and Diagnostics
#'
#' This primary function conducts a robust MANOVA analysis. It performs assumption testing
#' (multivariate normality via Mardia's test and homogeneity of covariance via Box's M test),
#' fits the MANOVA model, runs univariate follow-up ANOVAs, computes partial eta squared effect sizes,
#' computes estimated marginal means and confidence intervals (if the `emmeans` package is installed),
#' generates diagnostic plots for each dependent variable, and creates a canonical discriminant analysis plot
#' (if the `candisc` package is available).
#'
#' @param formula A formula object for MANOVA (e.g., `cbind(Y1, Y2) ~ Group`).
#' @param data A data frame containing the data.
#' @param group_var (Optional) A string specifying the grouping variable. Defaults to `NULL`. If `NULL`
#'   and the right-hand side of the formula contains a single variable, that variable is used.
#' @param alpha Significance level (default is 0.05). Currently reserved for future extensions.
#' @param plotDiagnostics Logical indicating whether to generate diagnostic plots for each dependent variable.
#'   Default is `TRUE`.
#'
#' @return A list with the following components:
#'   \item{manova_fit}{The MANOVA model object.}
#'   \item{manova_summary}{Summary of the MANOVA using Pillai's Trace.}
#'   \item{aov_results}{A list of ANOVA tables from univariate analyses.}
#'   \item{effect_sizes}{A list of data frames with partial eta squared for each dependent variable.}
#'   \item{ci_results}{A list of estimated marginal means and confidence intervals (if `emmeans` is installed).}
#'   \item{normality_results}{Results of Mardia's multivariate normality tests by group.}
#'   \item{boxm_result}{Result of Box's M test for homogeneity of covariance matrices.}
#'   \item{diagnostic_plots}{A list of diagnostic plots (ggplot objects) for each dependent variable.}
#'   \item{candisc_plot}{A canonical discriminant analysis plot (if the `candisc` package is installed).}
#'
#' @details This function first checks for the required packages (`MVN`, `heplots`, `ggplot2`)
#' and then proceeds with the analysis. It automatically determines the grouping variable if not explicitly provided.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Example with the iris dataset:
#'   res <- manova_analysis(cbind(Sepal.Length, Sepal.Width) ~ Species, data = iris)
#'   print(res$manova_summary)
#' }
manova_analysis <- function(formula, data, group_var = NULL, alpha = 0.05, plotDiagnostics = TRUE) {
  
  # Check required packages
  check_required_packages(c("MVN", "heplots", "ggplot2"))
  
  # Determine grouping variable (if possible)
  group_var <- get_group_variable(formula, data, group_var)
  
  # Run assumption tests if possible
  normality_results <- run_normality_tests(formula, data, group_var)
  boxm_result <- run_boxM_test(formula, data, group_var)
  
  # Run the MANOVA and obtain its summary
  manova_out <- run_manova(formula, data)
  manova_fit <- manova_out$fit
  manova_summary <- manova_out$summary
  
  # Run univariate analyses (ANOVAs)
  aov_results <- run_univariate_analyses(manova_fit)
  
  # Compute effect sizes for each dependent variable
  effect_sizes <- compute_effect_sizes(aov_results)
  
  # Compute estimated marginal means & CIs (if emmeans is available)
  ci_results <- compute_ci_estimates(formula, data, group_var)
  
  # Generate diagnostic plots if requested
  diagnostic_plots <- generate_diagnostic_plots(formula, data, plotDiagnostics)
  
  # Generate canonical discriminant analysis plot if possible
  candisc_plot <- generate_candisc_plot(manova_fit)
  
  # Return all results as a list
  results <- list(
    manova_fit        = manova_fit,
    manova_summary    = manova_summary,
    aov_results       = aov_results,
    effect_sizes      = effect_sizes,
    ci_results        = ci_results,
    normality_results = normality_results,
    boxm_result       = boxm_result,
    diagnostic_plots  = diagnostic_plots,
    candisc_plot      = candisc_plot
  )
  
  return(results)
}
