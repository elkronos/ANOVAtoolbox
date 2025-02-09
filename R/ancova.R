##############################################
# Enhanced Robust ANCOVA Analysis Code
##############################################

# 1. Load (or install if needed) required packages:
required_packages <- c("car", "ggplot2", "emmeans")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
suppressPackageStartupMessages({
  library(car)
  library(ggplot2)
  library(emmeans)
})

##############################################
# Helper Function: Validate Variables in Data
##############################################
#' Validate Variables in Data
#'
#' This helper function checks whether the specified variables exist in the given data frame.
#'
#' @param data A data frame containing the data.
#' @param vars A character vector of variable names to check for.
#'
#' @return If any variables in \code{vars} are missing from \code{data}, the function stops with an error message.
#' Otherwise, it returns \code{NULL}.
#'
#' @examples
#' data <- data.frame(a = 1:5, b = 6:10)
#' validate_vars(data, c("a", "b"))
#' # validate_vars(data, c("a", "c"))  # This will throw an error.
#'
#' @export
validate_vars <- function(data, vars) {
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("The following variables are missing from the data: ", 
         paste(missing_vars, collapse = ", "))
  }
}

##############################################
# Function: Check Homogeneity of Regression Slopes
##############################################
#' Check Homogeneity of Regression Slopes
#'
#' This function tests the assumption of homogeneity of regression slopes in an ANCOVA by fitting two models:
#' one with and one without the interaction between the covariate and the independent variable. It then compares
#' the models using a nested ANOVA to determine if the slopes are equal across groups.
#'
#' @param data A data frame containing the data.
#' @param dv A character string specifying the dependent variable.
#' @param iv A character string specifying the independent (grouping) variable.
#' @param covariate A character string specifying the covariate.
#' @param alpha Significance level for the test. Default is 0.05.
#'
#' @return A list containing:
#' \describe{
#'   \item{homogeneous}{Logical. \code{TRUE} if slopes are homogeneous (i.e., p-value >= \code{alpha}), \code{FALSE} otherwise.}
#'   \item{p_value}{The p-value from the nested ANOVA comparing the models.}
#'   \item{model_interaction}{The fitted linear model including the interaction term.}
#'   \item{model_no_interaction}{The fitted linear model without the interaction term.}
#' }
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   dv = rnorm(100),
#'   iv = rep(c("A", "B"), each = 50),
#'   covariate = rnorm(100)
#' )
#' result <- check_homogeneity_slopes(data, "dv", "iv", "covariate")
#' print(result$p_value)
#' }
#'
#' @export
check_homogeneity_slopes <- function(data, dv, iv, covariate, alpha = 0.05) {
  # Validate that variables exist
  validate_vars(data, c(dv, iv, covariate))
  
  # Validate that DV and covariate are numeric
  if (!is.numeric(data[[dv]])) {
    stop("Dependent variable (", dv, ") must be numeric.")
  }
  if (!is.numeric(data[[covariate]])) {
    stop("Covariate (", covariate, ") must be numeric.")
  }
  
  # Ensure independent variable is a factor (if not, coerce)
  if (!is.factor(data[[iv]])) {
    data[[iv]] <- as.factor(data[[iv]])
  }
  if (length(levels(data[[iv]])) < 2) {
    stop("Independent variable (", iv, ") must have at least two levels.")
  }
  
  # Create model formulas
  formula_interaction <- as.formula(paste(dv, "~", covariate, "*", iv))
  formula_no_interaction <- as.formula(paste(dv, "~", covariate, "+", iv))
  
  # Fit both models
  model_int <- lm(formula_interaction, data = data)
  model_no_int <- lm(formula_no_interaction, data = data)
  
  # Compare models using nested ANOVA
  anova_comp <- anova(model_no_int, model_int)
  p_value <- anova_comp$`Pr(>F)`[2]
  
  if (is.na(p_value)) {
    warning("Could not compute p-value for the homogeneity of slopes test. Please check your model specification.")
    return(list(homogeneous = NA, p_value = NA,
                model_interaction = model_int, model_no_interaction = model_no_int))
  }
  
  homogeneous <- (p_value >= alpha)
  
  return(list(
    homogeneous = homogeneous,
    p_value = p_value,
    model_interaction = model_int,
    model_no_interaction = model_no_int
  ))
}

##############################################
# Function: Check Normality of Residuals
##############################################
#' Check Normality of Model Residuals
#'
#' This function performs a Shapiro-Wilk test on the residuals of a fitted model to assess normality,
#' and generates a Q-Q plot using ggplot2.
#'
#' @param model A fitted model object (e.g., from \code{lm} or \code{aov}).
#'
#' @return A list containing:
#' \describe{
#'   \item{shapiro}{The result of the Shapiro-Wilk test (or \code{NULL} if not applicable).}
#'   \item{qq_plot}{A ggplot2 object displaying the Q-Q plot of the residuals.}
#' }
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ wt, data = mtcars)
#' norm_check <- check_normality(model)
#' print(norm_check$shapiro)
#' print(norm_check$qq_plot)
#' }
#'
#' @export
check_normality <- function(model) {
  res <- residuals(model)
  n <- length(res)
  if (n < 3) {
    warning("Too few residuals (n=", n, ") for Shapiro-Wilk normality test.")
    shapiro_result <- NULL
  } else if (n > 5000) {
    warning("Sample size (n=", n, ") exceeds 5000; Shapiro-Wilk test may not be appropriate.")
    shapiro_result <- shapiro.test(res)
  } else {
    shapiro_result <- shapiro.test(res)
  }
  
  qq_plot <- ggplot(data.frame(res = res), aes(sample = res)) +
    stat_qq(color = "blue", size = 2) +
    stat_qq_line(color = "red", linetype = "dashed") +
    ggtitle("Q-Q Plot of Residuals") +
    theme_minimal()
  
  return(list(shapiro = shapiro_result, qq_plot = qq_plot))
}

##############################################
# Function: Check Homogeneity of Variance (Levene’s Test)
##############################################
#' Check Homogeneity of Variance (Levene's Test)
#'
#' This function tests the assumption of equal variances across groups using Levene's Test.
#' It uses the residuals from a fitted model and a specified grouping variable.
#'
#' @param data A data frame containing the data.
#' @param model A fitted model object (e.g., from \code{lm}).
#' @param iv A character string specifying the independent (grouping) variable.
#'
#' @return The result of Levene's Test as returned by \code{car::leveneTest}.
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ factor(cyl), data = mtcars)
#' levene_result <- check_homogeneity_variance(mtcars, model, "cyl")
#' print(levene_result)
#' }
#'
#' @export
check_homogeneity_variance <- function(data, model, iv) {
  validate_vars(data, iv)
  
  res <- residuals(model)
  temp_df <- data.frame(res = res, group = data[[iv]])
  temp_df$group <- as.factor(temp_df$group)
  
  levene_result <- car::leveneTest(res ~ group, data = temp_df)
  return(levene_result)
}

##############################################
# Function: Calculate Effect Sizes (Partial Eta Squared)
##############################################
#' Calculate Partial Eta Squared Effect Sizes
#'
#' This function calculates partial eta squared effect sizes from a Type-III ANOVA table of a fitted model.
#'
#' @param model A fitted model object (e.g., from \code{lm}).
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{Effect}{The name of the effect (independent variable or interaction).}
#'   \item{Partial_Eta_Squared}{The partial eta squared effect size.}
#' }
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ wt + factor(cyl), data = mtcars)
#' effect_sizes <- calculate_effect_sizes(model)
#' print(effect_sizes)
#' }
#'
#' @export
calculate_effect_sizes <- function(model) {
  aov_table <- car::Anova(model, type = "III")
  ss <- aov_table$`Sum Sq`
  names(ss) <- rownames(aov_table)
  
  if (!"Residuals" %in% names(ss)) {
    stop("Residuals not found in the ANOVA table.")
  }
  
  ss_error <- ss["Residuals"]
  ss_effects <- ss[setdiff(names(ss), "Residuals")]
  
  partial_eta_squared <- ss_effects / (ss_effects + ss_error)
  
  return(data.frame(
    Effect = names(partial_eta_squared),
    Partial_Eta_Squared = partial_eta_squared,
    row.names = NULL
  ))
}

##############################################
# Function: Create Diagnostic Plots
##############################################
#' Create Diagnostic Plots for ANCOVA
#'
#' This function generates diagnostic plots for an ANCOVA analysis including:
#'   - A Residuals vs. Fitted values plot.
#'   - A Q-Q plot of model residuals.
#'   - A scatterplot of the dependent variable versus the covariate, colored by group with regression lines.
#'
#' @param model A fitted model object (e.g., from \code{lm}).
#' @param data A data frame containing the original data.
#' @param dv A character string specifying the dependent variable.
#' @param iv A character string specifying the independent (grouping) variable.
#' @param covariate A character string specifying the covariate.
#'
#' @return A list of ggplot2 objects:
#' \describe{
#'   \item{residuals_vs_fitted}{The Residuals vs. Fitted values plot.}
#'   \item{qq_plot}{The Q-Q plot of residuals.}
#'   \item{ancova_plot}{A scatterplot of the dependent variable versus the covariate, with group-specific regression lines.}
#' }
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ wt + factor(cyl), data = mtcars)
#' plots <- plot_diagnostics(model, mtcars, "mpg", "cyl", "wt")
#' print(plots$residuals_vs_fitted)
#' print(plots$qq_plot)
#' print(plots$ancova_plot)
#' }
#'
#' @export
plot_diagnostics <- function(model, data, dv, iv, covariate) {
  # Residuals vs. Fitted Plot
  diag1 <- ggplot(data.frame(Fitted = fitted(model), Residuals = residuals(model)),
                  aes(x = Fitted, y = Residuals)) +
    geom_point(color = "darkgreen") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggtitle("Residuals vs. Fitted") +
    xlab("Fitted values") +
    ylab("Residuals") +
    theme_minimal()
  
  # Q-Q Plot (using check_normality)
  norm_check <- check_normality(model)
  qq_plot <- norm_check$qq_plot
  
  # ANCOVA Scatterplot: DV vs. Covariate by Group
  ancova_plot <- ggplot(data, aes_string(x = covariate, y = dv, color = iv)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    ggtitle("ANCOVA: DV vs. Covariate by Group") +
    xlab(covariate) +
    ylab(dv) +
    theme_minimal()
  
  return(list(
    residuals_vs_fitted = diag1,
    qq_plot = qq_plot,
    ancova_plot = ancova_plot
  ))
}

##############################################
# Main Function: Run the ANCOVA Analysis
##############################################
#' Perform Enhanced Robust ANCOVA Analysis
#'
#' This function orchestrates a robust ANCOVA analysis by performing the following steps:
#' \enumerate{
#'   \item Validates the input data and variables.
#'   \item Checks homogeneity of regression slopes between the covariate and the dependent variable across groups.
#'   \item Fits the final ANCOVA model (with or without interaction, based on slope homogeneity).
#'   \item Assesses the normality of model residuals using a Shapiro-Wilk test and Q-Q plot.
#'   \item Tests the homogeneity of variances using Levene's Test.
#'   \item Calculates partial eta squared effect sizes.
#'   \item Computes estimated marginal means (EMMs) with 95% confidence intervals.
#'   \item Generates diagnostic plots including Residuals vs. Fitted, Q-Q plot, and an ANCOVA scatterplot.
#' }
#'
#' @param data A data frame containing the data for analysis.
#' @param dv A character string specifying the dependent variable.
#' @param iv A character string specifying the independent (grouping) variable.
#' @param covariate A character string specifying the covariate.
#' @param alpha Significance level for testing homogeneity of slopes. Default is 0.05.
#' @param plots Logical; if \code{TRUE}, diagnostic plots are generated and printed. Default is \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, progress messages and summaries are printed. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{final_model}{The final fitted ANCOVA model.}
#'   \item{homogeneity_slopes}{The result from the homogeneity of slopes test (including p-value and fitted models).}
#'   \item{normality}{The results of the normality test and Q-Q plot for model residuals.}
#'   \item{levene_test}{The result of Levene's Test for homogeneity of variances.}
#'   \item{effect_sizes}{A data frame of partial eta squared effect sizes.}
#'   \item{emm}{The estimated marginal means (EMMs) with 95% confidence intervals, if computed; otherwise, \code{NULL}.}
#'   \item{plots}{A list of diagnostic plots (Residuals vs. Fitted, Q-Q plot, and ANCOVA scatterplot).}
#' }
#'
#' @examples
#' \dontrun{
#' # Example data:
#' set.seed(123)
#' data_example <- data.frame(
#'   dv = rnorm(100, 50, 10),
#'   iv = rep(c("Group1", "Group2", "Group3"), length.out = 100),
#'   covariate = rnorm(100, 5, 2)
#' )
#'
#' # Run the ANCOVA analysis:
#' result <- ancova_analysis(data_example, "dv", "iv", "covariate")
#'
#' # Display the final model summary:
#' summary(result$final_model)
#'
#' # View diagnostic plots:
#' print(result$plots$residuals_vs_fitted)
#' print(result$plots$qq_plot)
#' print(result$plots$ancova_plot)
#' }
#'
#' @export
ancova_analysis <- function(data, dv, iv, covariate, alpha = 0.05,
                            plots = TRUE, verbose = TRUE) {
  # Ensure data is a data.frame
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame.")
  }
  
  # Validate that required variables are present
  validate_vars(data, c(dv, iv, covariate))
  
  # Validate numeric type for DV and covariate
  if (!is.numeric(data[[dv]])) {
    stop("Dependent variable (", dv, ") must be numeric.")
  }
  if (!is.numeric(data[[covariate]])) {
    stop("Covariate (", covariate, ") must be numeric.")
  }
  
  # Convert independent variable to factor if not already
  if (!is.factor(data[[iv]])) {
    data[[iv]] <- as.factor(data[[iv]])
  }
  if (length(levels(data[[iv]])) < 2) {
    stop("Independent variable (", iv, ") must have at least two levels.")
  }
  
  # Optionally remove rows with missing values in the required columns
  initial_rows <- nrow(data)
  data <- na.omit(data[, c(dv, iv, covariate)])
  if (nrow(data) < initial_rows && verbose) {
    message("Removed ", initial_rows - nrow(data), " rows with missing values.")
  }
  
  if (verbose) {
    message("=== Checking Homogeneity of Regression Slopes ===")
  }
  hs_result <- check_homogeneity_slopes(data, dv, iv, covariate, alpha)
  if (verbose) {
    message("p-value for interaction term: ", signif(hs_result$p_value, 4))
    if (hs_result$homogeneous) {
      message("Assumption met: slopes are homogeneous across groups.\n")
    } else {
      message("Assumption violated: slopes differ across groups.")
      message("Consider using the model with the interaction term.\n")
    }
  }
  
  # Choose final model: include interaction if slopes are heterogeneous
  if (isTRUE(hs_result$homogeneous)) {
    formula_final <- as.formula(paste(dv, "~", covariate, "+", iv))
  } else {
    formula_final <- as.formula(paste(dv, "~", covariate, "*", iv))
  }
  model_final <- lm(formula_final, data = data)
  
  if (verbose) {
    message("=== Final ANCOVA Model Summary ===")
    print(summary(model_final))
    message("")
  }
  
  # Check normality of residuals
  norm_result <- check_normality(model_final)
  if (verbose && !is.null(norm_result$shapiro)) {
    message("=== Shapiro-Wilk Normality Test ===")
    print(norm_result$shapiro)
    message("")
  }
  
  # Check homogeneity of variance using Levene's test
  levene_result <- check_homogeneity_variance(data, model_final, iv)
  if (verbose) {
    message("=== Levene's Test for Homogeneity of Variance ===")
    print(levene_result)
    message("")
  }
  
  # Calculate effect sizes (partial eta-squared)
  effect_sizes <- calculate_effect_sizes(model_final)
  if (verbose) {
    message("=== Effect Sizes (Partial Eta-Squared) ===")
    print(effect_sizes)
    message("")
  }
  
  # Compute Estimated Marginal Means (EMMs) with 95% CIs
  emm_results <- tryCatch({
    emmeans::emmeans(model_final, specs = iv, cov.keep = covariate)
  }, error = function(e) {
    warning("Could not compute estimated marginal means. Error: ", e$message)
    NULL
  })
  if (verbose && !is.null(emm_results)) {
    message("=== Estimated Marginal Means (with 95% CIs) ===")
    print(emm_results)
    message("")
  }
  
  # Generate diagnostic plots if requested
  plots_list <- NULL
  if (plots) {
    plots_list <- plot_diagnostics(model_final, data, dv, iv, covariate)
    if (verbose) {
      message("=== Diagnostic Plots ===")
      print(plots_list$residuals_vs_fitted)
      print(plots_list$qq_plot)
      print(plots_list$ancova_plot)
    }
  }
  
  # Return all analysis results as a list
  return(list(
    final_model = model_final,
    homogeneity_slopes = hs_result,
    normality = norm_result,
    levene_test = levene_result,
    effect_sizes = effect_sizes,
    emm = emm_results,
    plots = plots_list
  ))
}
