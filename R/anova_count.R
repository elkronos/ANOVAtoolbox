library(data.table)
library(stats)
library(ggplot2)
library(broom)
library(rlang)
library(multcomp)

#' Perform Deviance ANOVA for Count Data Using Poisson Regression
#'
#' This function fits a Poisson regression model to count data using an interaction of grouping variables,
#' performs a deviance analysis of variance (ANOVA), checks for overdispersion, runs Tukey‐adjusted post‐hoc tests,
#' calculates effect sizes (incidence rate ratios with confidence intervals), and generates a bar plot of total counts
#' per group.
#'
#' @param data A data frame or data.table containing the variables of interest.
#' @param response_var A character string specifying the name of the response variable (count data).
#' @param group_vars_list A character vector specifying the names of the grouping variables.
#' @param overdispersion_threshold A numeric value specifying the threshold above which overdispersion is flagged.
#'   Default is 1.5.
#' @param plot Logical. If \code{TRUE} (default), the count plot is printed.
#'
#' @return A list with the following elements:
#'   \item{model}{The fitted Poisson regression model object.}
#'   \item{overdispersion_statistic}{The overdispersion statistic (Pearson chi-square / df).}
#'   \item{overdispersion_check_message}{A message regarding overdispersion.}
#'   \item{deviance_anova}{The deviance ANOVA table.}
#'   \item{posthoc}{The post-hoc test results (Tukey contrasts).}
#'   \item{effect_size}{A tidy data frame of incidence rate ratios with confidence intervals.}
#'   \item{count_plot}{A ggplot2 object showing the total count by group.}
#'
#' @examples
#' \dontrun{
#' # Create an example dataset for a 2x2 factorial design
#' set.seed(1)
#' data_example <- data.frame(
#'   group1 = factor(rep(c("A", "B"), each = 100)),
#'   group2 = factor(rep(c("X", "Y"), times = 100)),
#'   count = c(rpois(100, 5), rpois(100, 10), rpois(100, 15), rpois(100, 20))
#' )
#'
#' # Run the analysis
#' results <- anova_count(data_example, response_var = "count", group_vars_list = c("group1", "group2"))
#'
#' # Check the post-hoc results
#' summary(results$posthoc)
#' }
anova_count <- function(data, response_var, group_vars_list, 
                        overdispersion_threshold = 1.5, plot = TRUE) {
  # ----- Input Validation -----------------------------------------------------
  # Check that data is provided and non-empty
  if (missing(data) || nrow(data) == 0) {
    stop("The data provided is empty. Please supply a non-empty data frame or data.table.")
  }
  
  # Check that the response variable exists
  if (!response_var %in% names(data)) {
    stop(sprintf("Response variable '%s' not found in the data.", response_var))
  }
  
  # Check that grouping variables exist
  if (!all(group_vars_list %in% names(data))) {
    missing_vars <- group_vars_list[!group_vars_list %in% names(data)]
    stop(sprintf("The following grouping variable(s) are missing from the data: %s",
                 paste(missing_vars, collapse = ", ")))
  }
  
  # Check that response variable is numeric and contains non-negative integers
  if (!is.numeric(data[[response_var]])) {
    stop(sprintf("Response variable '%s' must be numeric.", response_var))
  }
  if (any(data[[response_var]] < 0)) {
    stop(sprintf("Response variable '%s' must be non-negative.", response_var))
  }
  if (!all(round(data[[response_var]]) == data[[response_var]])) {
    stop(sprintf("Response variable '%s' must contain integer values (count data).", response_var))
  }
  
  # Convert to data.table if not already
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }
  
  # Ensure each grouping variable is a factor with at least 2 levels
  for (gv in group_vars_list) {
    data[[gv]] <- as.factor(data[[gv]])
    if (nlevels(data[[gv]]) < 2) {
      stop(sprintf("Grouping variable '%s' must have at least 2 levels.", gv))
    }
  }
  
  # ----- Create Interaction Term ----------------------------------------------
  # Create a new column 'group_interaction' representing the interaction of all grouping variables
  data[, group_interaction := do.call(interaction, .SD), .SDcols = group_vars_list]
  
  if (nlevels(data$group_interaction) < 2) {
    stop("The interaction of the grouping variables must have at least 2 levels.")
  }
  
  # ----- Fit Poisson Regression Model -----------------------------------------
  formula_str <- paste(response_var, "~ group_interaction")
  poisson_model <- glm(as.formula(formula_str), data = data, family = poisson())
  
  # ----- Check for Overdispersion ---------------------------------------------
  pearson_resid_sq <- sum(resid(poisson_model, type = "pearson")^2)
  df_resid <- df.residual(poisson_model)
  overdispersion_stat <- pearson_resid_sq / df_resid
  
  overdispersion_check_message <- if (overdispersion_stat > overdispersion_threshold) {
    sprintf("Warning: Overdispersion detected (dispersion statistic = %.2f). Consider using Negative Binomial regression.", 
            overdispersion_stat)
  } else {
    sprintf("No overdispersion detected (dispersion statistic = %.2f).", overdispersion_stat)
  }
  
  # ----- Deviance ANOVA -------------------------------------------------------
  deviance_anova <- anova(poisson_model, test = "Chisq")
  
  # ----- Post-hoc Tests ---------------------------------------------------------
  # Using Tukey's HSD for pairwise comparisons of the interaction levels.
  posthoc <- tryCatch({
    glht(poisson_model, linfct = mcp(group_interaction = "Tukey"))
  }, error = function(e) {
    warning("Post-hoc test failed: ", e$message)
    NULL
  })
  
  # ----- Calculate Effect Sizes -----------------------------------------------
  # Compute incidence rate ratios (IRRs) and their confidence intervals
  effect_size <- tidy(poisson_model, exponentiate = TRUE, conf.int = TRUE)
  
  # ----- Generate Count Plot --------------------------------------------------
  # Summarize total counts for each interaction group
  summary_data <- data[, .(total_count = sum(get(response_var))), by = group_interaction]
  
  count_plot <- ggplot(summary_data, aes(x = group_interaction, y = total_count, fill = group_interaction)) +
    geom_bar(stat = "identity", color = "black") +
    labs(x = paste(group_vars_list, collapse = " * "), 
         y = paste("Total", response_var),
         title = paste("Total", response_var, "by", paste(group_vars_list, collapse = " * "))) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold")
    )
  
  if (plot) {
    print(count_plot)
  }
  
  # ----- Return Results -------------------------------------------------------
  return(list(
    model = poisson_model,
    overdispersion_statistic = overdispersion_stat,
    overdispersion_check_message = overdispersion_check_message,
    deviance_anova = deviance_anova,
    posthoc = posthoc,
    effect_size = effect_size,
    count_plot = count_plot
  ))
}
