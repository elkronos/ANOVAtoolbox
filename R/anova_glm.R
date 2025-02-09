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
