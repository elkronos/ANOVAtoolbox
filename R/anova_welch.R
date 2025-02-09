#' Perform Welch's ANOVA with Assumption Checks, Post-hoc Tests, and Effect Sizes
#'
#' This function performs Welch's Analysis of Variance (ANOVA) for comparing a numeric response variable
#' across groups defined by one or more grouping variables. In addition, it evaluates the assumption of
#' normality of the residuals (using a ggplot2-based QQ plot and the Anderson-Darling test if sample size permits),
#' computes summary statistics (mean, standard deviation, standard error, and 95% confidence intervals),
#' performs pairwise post-hoc tests with Bonferroni correction (if applicable), and calculates Cohen's d for
#' pairwise group comparisons.
#'
#' **Note:** If the sample size (i.e., number of model residuals) is less than 8, the Anderson-Darling test
#' will be skipped because the test requires at least 8 observations.
#'
#' @param data A data frame containing the data for analysis.
#' @param response_var A character string specifying the name of the numeric response variable in \code{data}.
#' @param group_vars A character vector specifying the names of one or more grouping variables in \code{data}.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{assumptions}{A list containing:
#'       \describe{
#'         \item{qq_plot}{A ggplot2 object of the QQ plot for model residuals.}
#'         \item{ad_test}{The result of the Anderson-Darling normality test (or \code{NA} if sample size is too small).}
#'       }
#'     }
#'     \item{welch_anova}{The result of Welch's ANOVA (using \code{coin::oneway_test}).}
#'     \item{posthoc}{The result of post-hoc pairwise comparisons (if more than two groups are present).}
#'     \item{summary_stats}{A data.table containing the group-wise summary statistics (mean, SD, N, SE, and 95% CI).}
#'     \item{effect_sizes}{A list of Cohen's d effect sizes for all pairwise group comparisons.}
#'     \item{means_plot}{A ggplot2 barplot showing group means with 95% confidence intervals.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Example dataset
#' set.seed(123)
#' group <- rep(c("A", "B", "C"), each = 30)
#' value <- c(rnorm(30, 10, 2), rnorm(30, 12, 2.5), rnorm(30, 9, 1.5))
#' data_example <- data.frame(group, value)
#'
#' # Run Welch's ANOVA with the function
#' results <- anova_welch(data_example, "value", "group")
#'
#' # Display the QQ plot
#' print(results$assumptions$qq_plot)
#'
#' # Display Welch's ANOVA result
#' print(results$welch_anova)
#'
#' # Display group summary statistics
#' print(results$summary_stats)
#'
#' # Display the means plot
#' print(results$means_plot)
#' }
#'
#' @export
#' @import coin
#' @import ggplot2
#' @import data.table
#' @import nortest
#' @import effsize
anova_welch <- function(data, response_var, group_vars) {
  ## --- Input Validation ---
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data frame.")
  }
  if (!(response_var %in% names(data))) {
    stop(paste("Response variable", response_var, "not found in data."))
  }
  if (!all(group_vars %in% names(data))) {
    missing_vars <- group_vars[!group_vars %in% names(data)]
    stop(paste("The following grouping variables were not found in data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Ensure response_var is numeric
  if (!is.numeric(data[[response_var]])) {
    stop("The response variable must be numeric.")
  }
  
  # Remove rows with missing values in response_var or any group_vars
  data <- data[complete.cases(data[, c(response_var, group_vars)]), ]
  
  # Convert to data.table (if not already)
  data <- as.data.table(data)
  
  ## --- Create Combined Grouping Variable ---
  # Use a unique column name (here "group_interaction") to store the interaction of grouping variables.
  data[, group_interaction := do.call(interaction, c(.SD, sep = " : ")), .SDcols = group_vars]
  
  ## --- Model Fitting and Residuals ---
  model_formula <- as.formula(paste(response_var, "~ group_interaction"))
  aov_model <- aov(model_formula, data = data)
  model_residuals <- resid(aov_model)
  
  ## --- Assumption Checking: QQ Plot and Anderson-Darling Test ---
  # Create a ggplot2-based QQ plot for the residuals
  qq_data <- data.frame(residuals = model_residuals)
  qq_plot <- ggplot(qq_data, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line() +
    labs(title = "QQ Plot of Model Residuals", 
         x = "Theoretical Quantiles", 
         y = "Sample Quantiles") +
    theme_minimal()
  
  # Perform the Anderson-Darling test if sample size is at least 8; otherwise skip and issue a warning.
  if (length(model_residuals) < 8) {
    warning("Sample size too small for Anderson-Darling test. AD test not performed.")
    ad_test_result <- NA
  } else {
    ad_test_result <- ad.test(model_residuals)
  }
  
  ## --- Welch's ANOVA ---
  welch_anova_result <- oneway_test(model_formula, data = data, distribution = "asymptotic")
  
  ## --- Post-hoc Pairwise Comparisons ---
  # Only perform post-hoc tests if more than 2 groups are present.
  n_groups <- length(unique(data$group_interaction))
  posthoc_result <- NULL
  if (n_groups > 2) {
    posthoc_result <- pairwise.t.test(data[[response_var]], data$group_interaction,
                                      p.adjust.method = "bonferroni",
                                      pool.sd = FALSE,
                                      alternative = "two.sided")
  }
  
  ## --- Group-wise Summary Statistics ---
  # Compute mean, standard deviation, sample size, standard error, and 95% confidence intervals for each group.
  summary_stats <- data[, .(
    mean = mean(get(response_var), na.rm = TRUE),
    sd = sd(get(response_var), na.rm = TRUE),
    n = .N
  ), by = group_interaction]
  
  summary_stats[, se := sd / sqrt(n)]
  summary_stats[, `:=`(
    ci_low = mean - qt(0.975, n - 1) * se,
    ci_high = mean + qt(0.975, n - 1) * se
  )]
  
  ## --- Effect Sizes (Cohen's d) ---
  # Compute pairwise Cohen's d for all unique group combinations.
  groups <- unique(data$group_interaction)
  effect_sizes <- list()
  if (length(groups) > 1) {
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        grp1_values <- data[group_interaction == groups[i], get(response_var)]
        grp2_values <- data[group_interaction == groups[j], get(response_var)]
        d_val <- effsize::cohen.d(grp1_values, grp2_values)$estimate
        comp_name <- paste0(groups[i], " vs ", groups[j])
        effect_sizes[[comp_name]] <- d_val
      }
    }
  }
  
  ## --- Plot: Barplot of Group Means with 95% Confidence Intervals ---
  means_plot <- ggplot(summary_stats, aes(x = group_interaction, y = mean, fill = group_interaction)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2,
                  position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(mean, 2)), vjust = -0.5,
              position = position_dodge(width = 0.9), size = 3) +
    labs(x = "Group", y = "Mean", title = "Group Means with 95% Confidence Intervals") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ## --- Return Results ---
  results <- list(
    assumptions = list(
      qq_plot = qq_plot,
      ad_test = ad_test_result
    ),
    welch_anova = welch_anova_result,
    posthoc = posthoc_result,
    summary_stats = summary_stats,
    effect_sizes = effect_sizes,
    means_plot = means_plot
  )
  
  return(results)
}
