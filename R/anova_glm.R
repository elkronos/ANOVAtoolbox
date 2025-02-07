# Load required packages
library(data.table)
library(ggplot2)
library(rlang)
library(multcomp)
library(broom)
library(car)

#' Perform GLM and ANOVA Analysis
#'
#' This function fits a generalized linear model (GLM) using one or more grouping variables,
#' performs an ANOVA on the model, computes key model statistics including a pseudo R-squared,
#' and optionally produces diagnostic plots. If more than one grouping variable is provided,
#' they are combined into an interaction term.
#'
#' @param data A data frame or data table containing the data for analysis.
#' @param response_var A character string specifying the name of the response variable.
#' @param group_vars_vec A character vector of one or more grouping variable names.
#' @param family A family object specifying the GLM family (default is \code{gaussian()}).
#' @param sum_squares_type A character string indicating the sum of squares type: "II" (default)
#'   or "III". (Note: Type III requires appropriate contrasts, so ensure factors use \code{contr.sum}
#'   if needed.)
#' @param plot_residuals Logical. If \code{TRUE}, a ggplot residuals vs fitted values plot is generated.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{model}}{The fitted GLM object.}
#'     \item{\code{model_stats}}{A list of additional model statistics (AIC, BIC, etc.).}
#'     \item{\code{anova_test}}{A data frame of ANOVA test results.}
#'     \item{\code{posthoc}}{A summary of post-hoc Tukey comparisons (if applicable), or \code{NULL}.}
#'     \item{\code{effect_size}}{The pseudo R-squared (1 - deviance/null.deviance).}
#'     \item{\code{confidence_intervals}}{A data frame of confidence intervals for model coefficients.}
#'     \item{\code{residuals}}{The model residuals.}
#'     \item{\code{residuals_plot}}{A ggplot object of residuals vs fitted values (if requested).}
#'     \item{\code{boxplot}}{A ggplot boxplot of the response variable by group.}
#'   }
#'
#' @examples
#' # Example with a single grouping variable
#' set.seed(123)
#' group <- rep(c("Group A", "Group B", "Group C"), each = 50)
#' value <- c(rnorm(50, mean = 5, sd = 1),
#'            rnorm(50, mean = 7, sd = 1.5),
#'            rnorm(50, mean = 4, sd = 0.8))
#' data_single <- data.frame(group, value)
#'
#' results_single <- anova_glm(data_single, "value", c("group"),
#'                             family = gaussian(), sum_squares_type = "II", plot_residuals = TRUE)
#'
#' # Example with two grouping variables
#' group1 <- rep(c("Group A", "Group B"), each = 50)
#' group2 <- rep(c("Group X", "Group Y"), times = 50)
#' value <- c(rnorm(50, mean = 5, sd = 1),
#'            rnorm(50, mean = 7, sd = 1.5),
#'            rnorm(50, mean = 4, sd = 0.8),
#'            rnorm(50, mean = 6, sd = 1.2))
#' data_double <- data.frame(group1, group2, value)
#'
#' results_double <- anova_glm(data_double, "value", c("group1", "group2"),
#'                             family = gaussian(), sum_squares_type = "III", plot_residuals = TRUE)
#'
#' @export
anova_glm <- function(data, response_var, group_vars_vec,
                      family = gaussian(), sum_squares_type = "II",
                      plot_residuals = FALSE) {
  
  # --- Input Checks ---
  if (!response_var %in% names(data)) {
    stop("The response variable '", response_var, "' is not found in the data.")
  }
  
  missing_groups <- setdiff(group_vars_vec, names(data))
  if (length(missing_groups) > 0) {
    stop("The following grouping variable(s) are not found in the data: ",
         paste(missing_groups, collapse = ", "))
  }
  
  # Convert data to data.table for fast processing
  data <- as.data.table(data)
  
  # Ensure that each grouping variable is a factor; if not, convert and message the user.
  for (var in group_vars_vec) {
    if (!is.factor(data[[var]])) {
      data[[var]] <- as.factor(data[[var]])
      message(sprintf("Variable '%s' was converted to a factor.", var))
    }
  }
  
  # --- Create Predictor Variable ---
  # If there is a single grouping variable, use it directly.
  # Otherwise, combine the grouping variables into an interaction term.
  if (length(group_vars_vec) == 1) {
    data[, interaction_term := get(group_vars_vec[1])]
  } else {
    data[, interaction_term := do.call(interaction, .SD), .SDcols = group_vars_vec]
  }
  
  # --- Fit the GLM ---
  model_formula <- as.formula(paste(response_var, "~ interaction_term"))
  model <- glm(model_formula, data = data, family = family)
  
  # --- Diagnostics: Residuals ---
  model_resid <- resid(model)
  fitted_vals <- model$fitted.values
  
  # Create a ggplot residuals vs fitted values plot if requested
  residuals_plot <- NULL
  if (plot_residuals) {
    resid_df <- data.frame(Fitted = fitted_vals, Residuals = model_resid)
    residuals_plot <- ggplot(resid_df, aes(x = Fitted, y = Residuals)) +
      geom_point(color = "steelblue") +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      labs(title = "Residuals vs Fitted Values",
           x = "Fitted Values",
           y = "Residuals") +
      theme_minimal()
  }
  
  # --- ANOVA Analysis ---
  if (toupper(sum_squares_type) == "III") {
    # Note: For type III, ensure that contrasts are set appropriately.
    anova_test <- car::Anova(model, type = 3)
  } else if (toupper(sum_squares_type) == "II") {
    anova_test <- anova(model, test = "Chisq")
  } else {
    stop("Unsupported sum_squares_type. Please choose either 'II' or 'III'.")
  }
  
  # --- Compute Effect Size and Other Model Statistics ---
  pseudo_r2 <- 1 - (model$deviance / model$null.deviance)
  model_stats <- list(
    AIC = AIC(model),
    BIC = BIC(model),
    NullDeviance = model$null.deviance,
    ResidualDeviance = model$deviance,
    DF_null = model$df.null,
    DF_residual = model$df.residual
  )
  
  # --- Post-hoc Analysis ---
  posthoc <- NULL
  if (length(unique(data$interaction_term)) > 2) {
    posthoc_model <- tryCatch({
      multcomp::glht(model, linfct = mcp(interaction_term = "Tukey"))
    }, error = function(e) {
      warning("Post-hoc analysis could not be performed: ", e$message)
      return(NULL)
    })
    if (!is.null(posthoc_model)) {
      posthoc <- summary(posthoc_model)
    }
  }
  
  # --- Confidence Intervals ---
  conf_intervals <- as.data.frame(confint(model))
  colnames(conf_intervals) <- c("2.5 %", "97.5 %")
  
  # --- Create Boxplot of Response by Group ---
  # For a single grouping variable: x axis is that variable.
  # For multiple grouping variables: use the first variable for x and facet by the rest.
  if (length(group_vars_vec) == 1) {
    plot_box <- ggplot(data, aes_string(x = group_vars_vec[1], y = response_var,
                                        fill = group_vars_vec[1])) +
      geom_boxplot(show.legend = FALSE) +
      scale_fill_brewer(palette = "Dark2") +
      labs(title = "Boxplot of Response by Group",
           subtitle = paste("Response Variable:", response_var,
                            " | Group Variable:", group_vars_vec[1]),
           x = group_vars_vec[1],
           y = response_var) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12)
      )
  } else {
    facet_formula <- as.formula(paste(paste(group_vars_vec[-1], collapse = " + "), "~ ."))
    plot_box <- ggplot(data, aes_string(x = group_vars_vec[1], y = response_var,
                                        fill = group_vars_vec[1])) +
      geom_boxplot(show.legend = FALSE) +
      scale_fill_brewer(palette = "Dark2") +
      labs(title = "Boxplot of Response by Group",
           subtitle = paste("Response Variable:", response_var,
                            " | Group Variables:", paste(group_vars_vec, collapse = ", ")),
           x = group_vars_vec[1],
           y = response_var) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12)
      ) +
      facet_grid(facet_formula)
  }
  
  # --- Return Results ---
  results <- list(
    model = model,
    model_stats = model_stats,
    anova_test = as.data.frame(anova_test),
    posthoc = posthoc,
    effect_size = pseudo_r2,
    confidence_intervals = conf_intervals,
    residuals = model_resid,
    residuals_plot = residuals_plot,
    boxplot = plot_box
  )
  
  return(results)
}
