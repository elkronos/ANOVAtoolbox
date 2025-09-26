##############################################
# ANCOVA Analysis Toolkit
##############################################

# Package management -----------------------------------------------------------
required_packages <- c("car", "ggplot2", "emmeans", "withr", "sandwich", "lmtest", "rlang")
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
  library(rlang)
})

#' Validate the presence of variables in a data frame
#'
#' @title Validate Variables
#' @description
#' Check that specified variables are present in a data frame. If any are missing,
#' stop with an informative error message.
#'
#' @param data A data frame to inspect.
#' @param vars A character vector of variable names expected in \code{data}.
#'
#' @return Invisibly returns \code{NULL} if all variables are present; otherwise
#' an error is thrown.
#'
#' @examples
#' df <- data.frame(a = 1:3, b = 4:6)
#' validate_vars(df, c("a", "b"))
#' # validate_vars(df, c("a", "c")) # errors
#'
#' @export
validate_vars <- function(data, vars) {
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("The following variables are missing from the data: ",
         paste(missing_vars, collapse = ", "))
  }
  invisible(NULL)
}

#' Ensure a grouping variable is a factor with at least two levels
#'
#' @title Ensure Factor Grouping Variable
#' @description
#' Convert a grouping variable to a factor (if needed) and verify it has at least two levels.
#'
#' @param data A data frame containing the grouping variable.
#' @param iv A single character string giving the name of the grouping variable.
#'
#' @return The input \code{data} with \code{iv} coerced to a factor if not already.
#'
#' @examples
#' df <- data.frame(y = rnorm(6), g = c(1,1,1,2,2,2))
#' df2 <- ensure_factor_iv(df, "g")
#' is.factor(df2$g)
#'
#' @export
ensure_factor_iv <- function(data, iv) {
  if (!is.factor(data[[iv]])) data[[iv]] <- as.factor(data[[iv]])
  if (nlevels(data[[iv]]) < 2) {
    stop("Independent variable (", iv, ") must have at least two levels.")
  }
  data
}

#' Test homogeneity of regression slopes for ANCOVA
#'
#' @title Homogeneity of Regression Slopes
#' @description
#' Fit two linear models—one with and one without the interaction between a covariate
#' and a grouping variable—and compare them using a nested ANOVA to evaluate whether
#' regression slopes are homogeneous across groups.
#'
#' @param data A data frame containing the variables.
#' @param dv Character string for the dependent variable (numeric).
#' @param iv Character string for the independent (grouping) variable (factorable).
#' @param covariate Character string for the covariate (numeric).
#' @param alpha Numeric significance level used to flag homogeneity. Default is \code{0.05}.
#'
#' @details
#' Models are fit under sum-to-zero contrasts (\code{contr.sum}, \code{contr.poly})
#' to align with Type-III logic for effect assessment when factors are unbalanced.
#'
#' @return A list with elements:
#' \describe{
#'   \item{homogeneous}{Logical indicating whether slopes are considered homogeneous (\code{p >= alpha}).}
#'   \item{p_value}{P-value for the interaction comparison in the nested ANOVA.}
#'   \item{model_interaction}{The fitted model including the covariate-by-group interaction.}
#'   \item{model_no_interaction}{The fitted model without the interaction.}
#' }
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   y = rnorm(100),
#'   g = rep(letters[1:3], length.out = 100),
#'   x = rnorm(100)
#' )
#' out <- check_homogeneity_slopes(d, "y", "g", "x")
#' out$p_value
#'
#' @export
check_homogeneity_slopes <- function(data, dv, iv, covariate, alpha = 0.05) {
  validate_vars(data, c(dv, iv, covariate))
  if (!is.numeric(data[[dv]]))       stop("Dependent variable (", dv, ") must be numeric.")
  if (!is.numeric(data[[covariate]])) stop("Covariate (", covariate, ") must be numeric.")
  data <- ensure_factor_iv(data, iv)
  
  f_int <- as.formula(paste(dv, "~", covariate, "*", iv))
  f_no  <- as.formula(paste(dv, "~", covariate, "+", iv))
  
  an_cmp <- withr::with_options(
    list(contrasts = c("contr.sum", "contr.poly")),
    {
      m_int <- lm(f_int, data = data)
      m_no  <- lm(f_no,  data = data)
      anova(m_no, m_int)
    }
  )
  
  p_value <- suppressWarnings(an_cmp$`Pr(>F)`[2])
  if (is.na(p_value)) {
    warning("Could not compute a p-value for the homogeneity of slopes test. Check model specification.")
    return(list(homogeneous = NA, p_value = NA,
                model_interaction = m_int, model_no_interaction = m_no))
  }
  
  list(
    homogeneous = (p_value >= alpha),
    p_value = p_value,
    model_interaction = m_int,
    model_no_interaction = m_no
  )
}

#' Assess normality of model residuals
#'
#' @title Residual Normality Check
#' @description
#' Perform a Shapiro–Wilk test (when applicable) on residuals from a fitted model and
#' produce a Q–Q plot using \pkg{ggplot2}.
#'
#' @param model A fitted model object (e.g., from \code{lm}).
#'
#' @details
#' The Shapiro–Wilk test is carried out only for sample sizes from 3 to 5000 (inclusive).
#' For other sizes, the p-value is not computed and the plot is still returned.
#'
#' @return A list with elements:
#' \describe{
#'   \item{shapiro}{An object returned by \code{shapiro.test} or \code{NULL} if not applicable.}
#'   \item{qq_plot}{A \pkg{ggplot2} object showing the Q–Q plot of residuals.}
#' }
#'
#' @examples
#' m <- lm(mpg ~ wt, data = mtcars)
#' nchk <- check_normality(m)
#' nchk$shapiro
#' # print(nchk$qq_plot)
#'
#' @export
check_normality <- function(model) {
  res <- residuals(model)
  n <- length(res)
  
  if (n < 3 || n > 5000) {
    if (n < 3) {
      warning("Too few residuals (n = ", n, ") for the Shapiro–Wilk test.")
    } else {
      warning("Sample size (n = ", n, ") exceeds 5000; the Shapiro–Wilk test is not applied.")
    }
    shapiro_result <- NULL
  } else {
    shapiro_result <- shapiro.test(res)
  }
  
  qq_plot <- ggplot(data.frame(res = res), aes(sample = res)) +
    stat_qq() +
    stat_qq_line(linetype = "dashed") +
    ggtitle("Q–Q Plot of Residuals") +
    theme_minimal()
  
  list(shapiro = shapiro_result, qq_plot = qq_plot)
}

#' Test equality of variances across groups using Levene's test
#'
#' @title Levene's Test on Model Residuals
#' @description
#' Apply Levene's test on residuals from a fitted model across levels of a grouping variable.
#'
#' @param data A data frame containing the grouping variable.
#' @param model A fitted model object (e.g., from \code{lm}).
#' @param iv Character string for the independent (grouping) variable.
#'
#' @return The \code{car::leveneTest} result object.
#'
#' @examples
#' m <- lm(mpg ~ factor(cyl), data = mtcars)
#' check_homogeneity_variance(mtcars, m, "cyl")
#'
#' @export
check_homogeneity_variance <- function(data, model, iv) {
  validate_vars(data, iv)
  tmp <- data.frame(res = residuals(model), group = as.factor(data[[iv]]))
  car::leveneTest(res ~ group, data = tmp)
}

#' Type-III ANOVA with optional HC covariance inference
#'
#' @title Type-III ANOVA Table
#' @description
#' Compute a Type-III ANOVA table for a linear model and, optionally, a
#' heteroscedasticity-consistent (HC) adjusted Type-III table using
#' \code{car::Anova}'s \code{white.adjust} argument.
#'
#' @param model A fitted linear model.
#' @param use_hc Logical; if \code{TRUE}, compute a Type-III table with HC
#' covariance adjustment using \code{white.adjust}.
#' @param hc_type Character string for the HC estimator type. Typical values are
#' \code{"HC0"}, \code{"HC1"}, \code{"HC2"}, \code{"HC3"}; case-insensitive.
#'
#' @details
#' The classical Type-III table is obtained from \code{car::Anova(model, type = "III")}.
#' When \code{use_hc = TRUE}, the function additionally computes
#' \code{car::Anova(model, type = "III", white.adjust = <hc>)} where
#' \code{<hc>} is the lowercase version of \code{hc_type} (e.g., \code{"hc3"}).
#'
#' @return A list with elements:
#' \describe{
#'   \item{classical}{The \code{car::Anova} Type-III table.}
#'   \item{hc_tests}{A data-frame version of the HC-adjusted Type-III table, or
#'   \code{NULL} if \code{use_hc = FALSE}.}
#' }
#'
#' @examples
#' m <- withr::with_options(
#'   list(contrasts = c("contr.sum", "contr.poly")),
#'   lm(mpg ~ wt + factor(cyl), data = mtcars)
#' )
#' a3 <- anova_type3(m, use_hc = TRUE, hc_type = "HC3")
#' a3$classical
#' a3$hc_tests
#'
#' @export
anova_type3 <- function(model, use_hc = FALSE, hc_type = "HC3") {
  # Classical Type-III ANOVA
  aov3 <- car::Anova(model, type = "III")
  
  if (!use_hc) {
    return(list(classical = aov3, hc_tests = NULL))
  }
  
  # HC-adjusted Type-III ANOVA via white.adjust (expects lowercase)
  hc_arg <- tolower(hc_type)
  # car::Anova will error if an unsupported white.adjust is provided
  aov3_hc <- car::Anova(model, type = "III", white.adjust = hc_arg)
  
  list(classical = aov3, hc_tests = as.data.frame(aov3_hc))
}

#' Compute partial eta squared effect sizes from Type-III sums of squares
#'
#' @title Partial Eta Squared (Type-III)
#' @description
#' Calculate partial eta squared values for model effects using Type-III sums of squares
#' for the effects and the residual sum of squares from the classical ANOVA table.
#'
#' @param model A fitted linear model.
#'
#' @details
#' Effect sums of squares are taken from \code{car::Anova(model, type = "III")}.
#' The residual sum of squares is taken from \code{stats::anova(model)}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{Effect}{Effect name.}
#'   \item{Partial_Eta_Squared}{Partial eta squared value for the effect.}
#' }
#'
#' @examples
#' m <- withr::with_options(
#'   list(contrasts = c("contr.sum", "contr.poly")),
#'   lm(mpg ~ wt + factor(cyl), data = mtcars)
#' )
#' calculate_effect_sizes(m)
#'
#' @export
calculate_effect_sizes <- function(model) {
  aov3 <- car::Anova(model, type = "III")
  ss_eff <- aov3$`Sum Sq`
  names(ss_eff) <- rownames(aov3)
  
  aov1 <- stats::anova(model)
  if (!"Residuals" %in% rownames(aov1)) {
    stop("Residuals row not found in stats::anova(model).")
  }
  ss_err <- aov1["Residuals", "Sum Sq"]
  
  if ("(Intercept)" %in% names(ss_eff)) {
    ss_eff <- ss_eff[setdiff(names(ss_eff), "(Intercept)")]
  }
  
  pes <- ss_eff / (ss_eff + ss_err)
  data.frame(Effect = names(pes),
             Partial_Eta_Squared = as.numeric(pes),
             row.names = NULL)
}

#' Create diagnostic plots for an ANCOVA model
#'
#' @title ANCOVA Diagnostics
#' @description
#' Generate three diagnostic plots: Residuals vs Fitted, Q–Q plot of residuals,
#' and a scatterplot of the dependent variable versus the covariate colored by group
#' with fitted lines from separate regressions.
#'
#' @param model A fitted linear model.
#' @param data The original data frame used to fit the model.
#' @param dv Character string naming the dependent variable.
#' @param iv Character string naming the grouping variable.
#' @param covariate Character string naming the covariate.
#'
#' @return A list containing:
#' \describe{
#'   \item{residuals_vs_fitted}{A \pkg{ggplot2} plot of residuals versus fitted values.}
#'   \item{qq_plot}{A \pkg{ggplot2} Q–Q plot of residuals.}
#'   \item{ancova_plot}{A \pkg{ggplot2} scatterplot of \code{dv} vs \code{covariate} by \code{iv} with linear fits.}
#' }
#'
#' @examples
#' m <- lm(mpg ~ wt + factor(cyl), data = mtcars)
#' p <- plot_diagnostics(m, mtcars, "mpg", "cyl", "wt")
#' # print(p$residuals_vs_fitted); print(p$qq_plot); print(p$ancova_plot)
#'
#' @export
plot_diagnostics <- function(model, data, dv, iv, covariate) {
  diag1 <- ggplot(data.frame(Fitted = fitted(model), Residuals = residuals(model)),
                  aes(x = Fitted, y = Residuals)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle("Residuals vs Fitted") +
    xlab("Fitted values") +
    ylab("Residuals") +
    theme_minimal()
  
  norm_check <- check_normality(model)
  qq_plot <- norm_check$qq_plot
  
  ancova_plot <- ggplot(data, aes_string(x = covariate, y = dv, color = iv)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    ggtitle("Dependent Variable vs Covariate by Group") +
    xlab(covariate) +
    ylab(dv) +
    theme_minimal()
  
  list(
    residuals_vs_fitted = diag1,
    qq_plot = qq_plot,
    ancova_plot = ancova_plot
  )
}

#' Compute estimated marginal means for the grouping variable
#'
#' @title Estimated Marginal Means (HC Option)
#' @description
#' Compute estimated marginal means (EMMs) for the grouping factor, optionally using
#' a heteroscedasticity-consistent covariance matrix when forming standard errors and intervals.
#'
#' @param model A fitted linear model.
#' @param iv Character string for the grouping variable.
#' @param cov_reduce A function used by \code{emmeans} to reduce covariates; default is \code{mean}.
#' @param use_hc Logical; if \code{TRUE}, pass an HC covariance matrix to \code{emmeans} via \code{vcov.}.
#' @param hc_type Character string for the HC estimator type passed to \code{vcovHC}. Default is \code{"HC3"}.
#'
#' @return An \code{emmGrid} object containing estimated marginal means and interval estimates.
#'
#' @examples
#' m <- withr::with_options(
#'   list(contrasts = c("contr.sum", "contr.poly")),
#'   lm(mpg ~ wt + factor(cyl), data = mtcars)
#' )
#' compute_emm(m, "cyl", use_hc = TRUE)
#'
#' @export
compute_emm <- function(model, iv, cov_reduce = mean, use_hc = FALSE, hc_type = "HC3") {
  specs_formula <- as.formula(paste("~", iv))
  if (!use_hc) {
    return(emmeans::emmeans(model, specs = specs_formula, cov.reduce = cov_reduce))
  }
  vc <- sandwich::vcovHC(model, type = hc_type)
  emmeans::emmeans(model, specs = specs_formula, cov.reduce = cov_reduce, vcov. = vc)
}

#' Conduct an ANCOVA workflow with diagnostics, effect sizes, and EMMs
#'
#' @title ANCOVA Workflow
#' @description
#' Perform an ANCOVA workflow that includes variable validation, homogeneity of regression slopes,
#' model fitting under contrasts suitable for Type-III inference, residual diagnostics,
#' Levene's test for equality of variances, Type-III ANOVA, partial eta squared effect sizes,
#' and estimated marginal means for the grouping factor.
#'
#' @param data A data frame containing the variables.
#' @param dv Character string for the dependent variable (numeric).
#' @param iv Character string for the independent (grouping) variable.
#' @param covariate Character string for the covariate (numeric).
#' @param alpha Numeric significance level used for the slopes homogeneity decision. Default \code{0.05}.
#' @param use_hc Logical; if \code{TRUE}, perform per-effect tests and EMMs with an HC covariance matrix.
#' @param hc_type Character string for the HC estimator type passed to \code{vcovHC}. Default \code{"HC3"}.
#' @param drop_na Logical; if \code{TRUE}, rows with missing values in the required variables are removed.
#' @param plots Logical; if \code{TRUE}, create diagnostic plots.
#' @param verbose Logical; if \code{TRUE}, print progress messages and summaries.
#'
#' @details
#' The final model includes the covariate-by-group interaction when the slopes homogeneity
#' test indicates differential slopes; otherwise, the interaction is omitted.
#' Models are fit under sum-to-zero contrasts for Type-III inference.
#'
#' @return A list containing:
#' \describe{
#'   \item{final_model}{The fitted linear model selected by the workflow.}
#'   \item{homogeneity_slopes}{Output from \code{check_homogeneity_slopes}.}
#'   \item{normality}{Output from \code{check_normality}.}
#'   \item{levene_test}{Output from \code{check_homogeneity_variance}.}
#'   \item{anova_type3}{List with \code{$classical} Type-III table and optional \code{$hc_tests}.}
#'   \item{effect_sizes}{Data frame of partial eta squared values.}
#'   \item{emm}{Estimated marginal means for the grouping factor.}
#'   \item{plots}{List of diagnostic plots, or \code{NULL} if \code{plots = FALSE}.}
#' }
#'
#' @examples
#' set.seed(123)
#' d <- data.frame(
#'   dv = rnorm(120, 50, 10),
#'   iv = rep(c("Group1", "Group2", "Group3"), length.out = 120),
#'   covariate = rnorm(120, 5, 2)
#' )
#' res <- ancova_analysis(d, "dv", "iv", "covariate",
#'                        use_hc = TRUE, plots = TRUE, verbose = TRUE)
#' summary(res$final_model)
#' res$effect_sizes
#' res$anova_type3$classical
#' res$anova_type3$hc_tests
#'
#' @export
ancova_analysis <- function(data, dv, iv, covariate,
                            alpha = 0.05,
                            use_hc = FALSE,
                            hc_type = "HC3",
                            drop_na = TRUE,
                            plots = TRUE,
                            verbose = TRUE) {
  
  if (!is.data.frame(data)) stop("Input 'data' must be a data.frame.")
  validate_vars(data, c(dv, iv, covariate))
  if (!is.numeric(data[[dv]]))       stop("Dependent variable (", dv, ") must be numeric.")
  if (!is.numeric(data[[covariate]])) stop("Covariate (", covariate, ") must be numeric.")
  data <- ensure_factor_iv(data, iv)
  
  if (drop_na) {
    initial_rows <- nrow(data)
    data <- stats::na.omit(data[, c(dv, iv, covariate)])
    if (verbose && nrow(data) < initial_rows) {
      message("Removed ", initial_rows - nrow(data), " rows with missing values.")
    }
  }
  
  if (verbose) message("=== Checking homogeneity of regression slopes ===")
  hs <- check_homogeneity_slopes(data, dv, iv, covariate, alpha)
  if (verbose) {
    message("Interaction p-value: ", ifelse(is.na(hs$p_value), "NA", signif(hs$p_value, 4)))
    message(if (isTRUE(hs$homogeneous))
      "Using additive model (no interaction)."
      else
        "Using interaction model.")
    message("")
  }
  
  f_final <- if (isTRUE(hs$homogeneous)) {
    as.formula(paste(dv, "~", covariate, "+", iv))
  } else {
    as.formula(paste(dv, "~", covariate, "*", iv))
  }
  
  model_final <- withr::with_options(
    list(contrasts = c("contr.sum", "contr.poly")),
    lm(f_final, data = data)
  )
  
  if (verbose) {
    message("=== Model summary ===")
    print(summary(model_final))
    message("")
  }
  
  norm_result <- check_normality(model_final)
  if (verbose) {
    message("=== Shapiro–Wilk test for residuals ===")
    if (is.null(norm_result$shapiro)) {
      message("Not applied (n < 3 or n > 5000).")
    } else {
      print(norm_result$shapiro)
    }
    message("")
  }
  
  lev_result <- check_homogeneity_variance(data, model_final, iv)
  if (verbose) {
    message("=== Levene's test ===")
    print(lev_result)
    message("")
  }
  
  if (verbose) message("=== Type-III ANOVA ===")
  aov_out <- anova_type3(model_final, use_hc = use_hc, hc_type = hc_type)
  if (verbose) {
    message("--- Classical (Type III) ---")
    print(aov_out$classical)
    if (!is.null(aov_out$hc_tests)) {
      message("\n--- HC covariance per-effect tests ---")
      print(aov_out$hc_tests)
    }
    message("")
  }
  
  eff_sizes <- calculate_effect_sizes(model_final)
  if (verbose) {
    message("=== Partial eta squared ===")
    print(eff_sizes)
    message("")
  }
  
  emm_results <- tryCatch({
    compute_emm(model_final, iv, cov_reduce = mean, use_hc = use_hc, hc_type = hc_type)
  }, error = function(e) {
    warning("Could not compute estimated marginal means. Error: ", e$message)
    NULL
  })
  if (verbose && !is.null(emm_results)) {
    message("=== Estimated marginal means ===")
    print(emm_results)
    message("")
  }
  
  plots_list <- NULL
  if (plots) {
    plots_list <- plot_diagnostics(model_final, data, dv, iv, covariate)
    if (verbose) {
      message("=== Diagnostic plots ===")
      print(plots_list$residuals_vs_fitted)
      print(plots_list$qq_plot)
      print(plots_list$ancova_plot)
      message("")
    }
  }
  
  list(
    final_model        = model_final,
    homogeneity_slopes = hs,
    normality          = norm_result,
    levene_test        = lev_result,
    anova_type3        = aov_out,
    effect_sizes       = eff_sizes,
    emm                = emm_results,
    plots              = plots_list
  )
}

# Example (commented) ----------------------------------------------------------
# set.seed(123)
# d <- data.frame(
#   dv = rnorm(120, 50, 10),
#   iv = rep(c("Group1", "Group2", "Group3"), length.out = 120),
#   covariate = rnorm(120, 5, 2)
# )
# res <- ancova_analysis(d, "dv", "iv", "covariate",
#                        use_hc = TRUE, plots = TRUE, verbose = TRUE)
# summary(res$final_model)
# res$effect_sizes
# res$anova_type3$classical
# res$anova_type3$hc_tests
