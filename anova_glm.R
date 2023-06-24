# Load required packages
library(data.table)
library(ggplot2)
library(rlang)
library(multcomp)
library(broom)
library(car)

#' Perform GLM and ANOVA Analysis
#'
#' This function performs a Generalized Linear Model (GLM) analysis with ANOVA and provides various outputs such as residuals, ANOVA test results, effect size, confidence intervals, and plots.
#'
#' @param data A data frame or data table containing the data for analysis.
#' @param response_var The name of the response variable in the data.
#' @param group_vars_vec A character vector specifying the names of the grouping variables in the data.
#' @param family A family object specifying the type of GLM model to be fitted (default is "gaussian()").
#' @param sum_squares_type The type of sum of squares to be used in the ANOVA test ("II" or "III") (default is "II").
#' @param plot_residuals A logical value indicating whether to generate a Residuals vs Fitted plot (default is FALSE).
#'
#' @return A list of results including the residuals plot, ANOVA test results, posthoc analysis (if applicable), effect size, confidence intervals, and a plot.
#'
#' @importFrom data.table data.table
#' @importFrom ggplot2 ggplot, geom_boxplot, scale_fill_brewer, labs, theme_minimal, theme, element_text, facet_grid
#' @importFrom rlang sym
#' @importFrom multcomp glht, mcp, summary
#' @importFrom broom tidy
#' @importFrom car Anova
#'
#' @examples
#' # Single factor example
#' group <- rep(c("Group A", "Group B", "Group C"), each = 50)
#' value <- c(rnorm(50, mean = 5, sd = 1),
#'            rnorm(50, mean = 7, sd = 1.5),
#'            rnorm(50, mean = 4, sd = 0.8))
#'
#' data_single <- data.frame(group, value)
#' data_single$group <- as.factor(data_single$group)  # Convert to factor
#'
#' # Two-factor example
#' group1 <- rep(c("Group A", "Group B"), each = 50)
#' group2 <- rep(c("Group X", "Group Y"), times = 50)
#' value <- c(rnorm(50, mean = 5, sd = 1),
#'            rnorm(50, mean = 7, sd = 1.5),
#'            rnorm(50, mean = 4, sd = 0.8),
#'            rnorm(50, mean = 6, sd = 1.2))
#'
#' data_double <- data.frame(group1, group2, value)
#' data_double$group1 <- as.factor(data_double$group1)  # Convert group1 to factor
#' data_double$group2 <- as.factor(data_double$group2)  # Convert group2 to factor
#'
#' # Single factor test
#' results_single <- anova_glm(data_single, "value", c("group"), family = gaussian(), plot_residuals = TRUE)
#'
#' # Two-factor test
#' results_double <- anova_glm(data_double, "value", c("group1", "group2"), family = gaussian(), plot_residuals = TRUE)
#'
#' # Access the plot
#' print(results_single$plot)  # Print the single-factor plot
#' print(results_double$plot)  # Print the two-factor plot
#' 
#' @export
anova_glm <- function(data, response_var, group_vars_vec, family = gaussian(), sum_squares_type = "II", plot_residuals = FALSE) {
  
  # Check if response_var exists in the data
  if (!(response_var %in% names(data))) stop("response_var not found in data")
  
  # Check if all group_vars_vec elements exist in the data
  not_found_vars <- setdiff(group_vars_vec, names(data))
  if (length(not_found_vars) > 0) stop(paste("The following group_vars_vec elements are not found in data:", paste(not_found_vars, collapse = ", ")))
  
  # Convert data to data.table
  data <- data.table(data)
  
  # Create interaction term based on group variables
  interaction_vars <- data[, lapply(.SD, as.character), .SDcols = group_vars_vec]
  data[, 'interaction_term' := do.call(interaction, interaction_vars)]
  
  # Build a model formula and fit the GLM
  model_formula <- as.formula(paste(response_var, "~", "interaction_term"))
  model <- glm(model_formula, data = data, family = family)
  
  # Calculate the residuals
  residuals <- resid(model)
  
  # Generate Residuals vs Fitted plot if requested
  if (plot_residuals) {
    plot(model$fitted.values, residuals, 
         main="Residuals vs Fitted Values", 
         xlab="Fitted Values", 
         ylab="Residuals", 
         pch=19, frame=F, col="steelblue")
    abline(h=0, col="red", lwd=2)
  }
  
  # Conduct ANOVA on GLM
  if (sum_squares_type == "III") {
    anova_test <- Anova(model, type = sum_squares_type)
  } else {
    anova_test <- anova(model, test="Chisq")
  }
  
  # Effect size (pseudo R-squared)
  eff_size <- 1 - (model$deviance / model$null.deviance)
  
  # Conduct post-hoc test if necessary
  posthoc <- NULL
  if (length(unique(data[[group_vars_vec[1]]])) > 2 || length(group_vars_vec) > 1) {
    posthoc_res <- glht(model, linfct=mcp(interaction_term="Tukey"))
    posthoc <- summary(posthoc_res)  # Extract the summary dataframe from the output
  }
  
  # Confidence intervals
  conf_intervals <- as.data.frame(confint(model))
  conf_intervals <- setNames(conf_intervals, c("2.5 %", "97.5 %")) # rename the column names
  
  # Create a boxplot for each interaction term group, color-coded by the group
  plot <- ggplot(data, aes(x = !!rlang::sym(group_vars_vec[1]), y = !!rlang::sym(response_var))) +
    geom_boxplot(aes(fill = !!rlang::sym(group_vars_vec[1])), show.legend = FALSE) +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      title = "Boxplot of Response by Group", 
      subtitle = paste("Response Variable:", response_var, " Group Variables:", paste(group_vars_vec, collapse = ", ")),
      x = "Interaction Group",
      y = "Response Variable Value"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5, size = 15),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  # Add facetting for each additional grouping variable
  if (length(group_vars_vec) > 1) {
    for (var in group_vars_vec[-1]) {
      plot <- plot + facet_grid(paste(var, "~ ."))
    }
  }
  
  # Return a list of results including the residuals plot, ANOVA test results, posthoc analysis (if applicable), effect size and confidence intervals
  results <- list(
    assumptions = list(
      residuals_plot = residuals  # Note: this is not a plot but the residuals from the model
    ),
    anova_test = as.data.frame(anova_test),  # changed tidy(anova_test) to as.data.frame(anova_test)
    posthoc = posthoc,
    effect_size = eff_size,
    confidence_intervals = conf_intervals,
    plot = plot  # Add the plot to the results
  )
  
  return(results)
  
}