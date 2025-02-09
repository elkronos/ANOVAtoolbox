#' Perform ANOVA for a Generalized Linear Model with Grouping Variables
#'
#' This function fits a generalized linear model (GLM) to the provided dataset using one or more
#' grouping variables and performs an analysis of variance (ANOVA). It constructs the model formula
#' based on the provided grouping variables and whether a full factorial design is desired. Additionally,
#' the function computes model statistics, generates diagnostic plots (residuals and boxplot), and
#' optionally conducts post-hoc analyses using Tukey comparisons.
#'
#' @param data A \code{data.frame} (or \code{data.table}) containing the dataset.
#' @param response_var A character string specifying the name of the response variable column in \code{data}.
#' @param group_vars_vec A character vector specifying the names of the grouping variables in \code{data}.
#' @param family A family object specifying the error distribution and link function for the GLM.
#'   Defaults to \code{gaussian()}.
#' @param sum_squares_type A character string indicating the type of sum of squares to use in the ANOVA.
#'   Supported options are \code{"II"} (default) and \code{"III"}. When \code{"III"} is specified, the function
#'   employs \code{car::Anova}.
#' @param plot_residuals Logical. If \code{TRUE}, a residuals vs fitted values plot is generated using ggplot2.
#'   Defaults to \code{FALSE}.
#' @param full_factorial Logical. If \code{TRUE} and multiple grouping variables are provided, a full factorial
#'   model (with all interactions) is fitted. Otherwise, an interaction term is created to combine the groups.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' The function carries out the following steps:
#' \enumerate{
#'   \item \strong{Input Checks:} Validates that the \code{response_var} and each element in \code{group_vars_vec}
#'   exist as columns in \code{data}. If a grouping variable is not a factor, it is converted to one.
#'   \item \strong{Data Preparation:} Converts \code{data} to a \code{data.table} for efficient processing.
#'   \item \strong{Model Formula Construction:} Depending on the number of grouping variables and the
#'   \code{full_factorial} flag, constructs an appropriate model formula. For a single grouping variable, the
#'   model is simple; for multiple groups, either a full factorial or an interaction term is used.
#'   \item \strong{Model Fitting:} Fits the GLM using the constructed formula and the specified \code{family}.
#'   \item \strong{Diagnostics:} Computes residuals and fitted values. If \code{plot_residuals} is \code{TRUE},
#'   a residuals vs fitted values plot is created.
#'   \item \strong{ANOVA Analysis:} Performs ANOVA using either type II or type III sum of squares based on
#'   the \code{sum_squares_type} parameter.
#'   \item \strong{Model Statistics:} Calculates pseudo R-squared and other model statistics (AIC, BIC, deviance, etc.).
#'   \item \strong{Post-hoc Analysis:} When appropriate (and if not using a full factorial model), conducts a
#'   post-hoc Tukey comparison.
#'   \item \strong{Plotting:} Generates a boxplot of the response variable by group.
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{model}}{The fitted GLM object.}
#'   \item{\code{model_stats}}{A list of model statistics including AIC, BIC, null and residual deviances, and degrees of freedom.}
#'   \item{\code{anova_test}}{A data frame of the ANOVA table results.}
#'   \item{\code{posthoc}}{A summary of the post-hoc analysis (if performed), or \code{NULL} otherwise.}
#'   \item{\code{effect_size}}{The pseudo R-squared value calculated as \code{1 - (deviance / null.deviance)}.}
#'   \item{\code{confidence_intervals}}{A data frame of confidence intervals for the model coefficients.}
#'   \item{\code{residuals}}{The residuals from the GLM.}
#'   \item{\code{residuals_plot}}{A ggplot object showing residuals vs fitted values if \code{plot_residuals} is \code{TRUE}, else \code{NULL}.}
#'   \item{\code{boxplot}}{A ggplot object displaying a boxplot of the response variable by group.}
#' }
#'
#' @seealso \code{\link[stats]{glm}}, \code{\link[car]{Anova}}, \code{\link[multcomp]{glht}},
#'   \code{\link[data.table]{as.data.table}}, \code{\link[ggplot2]{ggplot}}
#'
#' @examples
#' \dontrun{
#'   # Load required libraries
#'   library(data.table)
#'   library(ggplot2)
#'   library(car)
#'   library(multcomp)
#'
#'   # Generate example data
#'   set.seed(123)
#'   example_data <- data.frame(
#'     response = rnorm(100),
#'     group1 = sample(letters[1:4], 100, replace = TRUE),
#'     group2 = sample(LETTERS[1:3], 100, replace = TRUE)
#'   )
#'
#'   # Perform the ANOVA GLM analysis with type II sum of squares and residuals plot
#'   result <- anova_glm(data = example_data,
#'                       response_var = "response",
#'                       group_vars_vec = c("group1", "group2"),
#'                       family = gaussian(),
#'                       sum_squares_type = "II",
#'                       plot_residuals = TRUE,
#'                       full_factorial = FALSE)
#'
#'   # Display the ANOVA table
#'   print(result$anova_test)
#'
#'   # If generated, display the residuals plot
#'   if (!is.null(result$residuals_plot)) {
#'     print(result$residuals_plot)
#'   }
#'
#'   # Display the boxplot of the response by group
#'   print(result$boxplot)
#' }
#'
#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom car Anova
#' @importFrom multcomp glht
#' @importFrom stats glm as.formula resid fitted AIC BIC confint
anova_glm <- function(data, response_var, group_vars_vec,
                      family = gaussian(), sum_squares_type = "II",
                      plot_residuals = FALSE, full_factorial = FALSE) {
  
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
  
  # --- Create Model Formula ---
  # Always create an interaction_term column (even if only one grouping variable)
  if (length(group_vars_vec) == 1) {
    data[, interaction_term := get(group_vars_vec[1])]
    model_formula <- as.formula(paste(response_var, "~", group_vars_vec[1]))
  } else {
    if (full_factorial) {
      # Build a full factorial formula including main effects and interactions
      predictors <- paste(group_vars_vec, collapse = " * ")
      model_formula <- as.formula(paste(response_var, "~", predictors))
    } else {
      data[, interaction_term := do.call(interaction, .SD), .SDcols = group_vars_vec]
      model_formula <- as.formula(paste(response_var, "~ interaction_term"))
    }
  }
  
  # --- Fit the GLM ---
  model <- glm(model_formula, data = data, family = family)
  
  # --- Diagnostics: Residuals ---
  model_resid <- resid(model)
  fitted_vals <- model$fitted.values
  
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
  if (!full_factorial) {
    if (length(unique(data$interaction_term)) > 2) {
      posthoc_model <- tryCatch({
        if (length(group_vars_vec) == 1) {
          # For a single grouping variable, refit without an intercept so that
          # the coefficients match the factor levels (required for Tukey comparisons)
          model_posthoc <- glm(as.formula(paste(response_var, "~ interaction_term - 1")),
                               data = data, family = family)
          linfct <- mcp(interaction_term = "Tukey")
          multcomp::glht(model_posthoc, linfct = linfct, vcov = vcov(model_posthoc))
        } else {
          # For multiple grouping variables (non-full-factorial), use the existing model.
          linfct <- mcp(interaction_term = "Tukey")
          multcomp::glht(model, linfct = linfct, vcov = vcov(model))
        }
      }, error = function(e) {
        warning("Post-hoc analysis could not be performed: ", e$message)
        return(NULL)
      })
      if (!is.null(posthoc_model)) {
        posthoc <- summary(posthoc_model)
      }
    }
  }
  
  # --- Confidence Intervals ---
  conf_intervals <- as.data.frame(confint(model))
  colnames(conf_intervals) <- c("2.5 %", "97.5 %")
  
  # --- Create Boxplot of Response by Group ---
  if (length(group_vars_vec) == 1) {
    plot_box <- ggplot(data, aes_string(x = group_vars_vec[1], y = response_var,
                                        fill = group_vars_vec[1])) +
      geom_boxplot(show.legend = FALSE) +
      scale_fill_brewer(palette = "Dark2") +
      labs(title = "Boxplot of Response by Group",
           subtitle = paste("Response Variable:", response_var,
                            "| Group Variable:", group_vars_vec[1]),
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
    if (full_factorial) {
      facet_formula <- as.formula(paste("~", paste(group_vars_vec[-1], collapse = " + ")))
      plot_box <- ggplot(data, aes_string(x = group_vars_vec[1], y = response_var,
                                          fill = group_vars_vec[1])) +
        geom_boxplot(show.legend = FALSE) +
        scale_fill_brewer(palette = "Dark2") +
        labs(title = "Boxplot of Response by Group",
             subtitle = paste("Response Variable:", response_var,
                              "| Group Variables:", paste(group_vars_vec, collapse = ", ")),
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
    } else {
      plot_box <- ggplot(data, aes_string(x = "interaction_term", y = response_var,
                                          fill = "interaction_term")) +
        geom_boxplot(show.legend = FALSE) +
        scale_fill_brewer(palette = "Dark2") +
        labs(title = "Boxplot of Response by Combined Group",
             subtitle = paste("Response Variable:", response_var,
                              "| Combined Group Variables:", paste(group_vars_vec, collapse = ", ")),
             x = "Combined Group (Interaction Term)",
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
    }
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
