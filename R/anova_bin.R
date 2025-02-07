#' One-Way "ANOVA" for a Binary Response via Logistic Regression
#'
#' This function performs an analysis akin to one-way ANOVA for a binary response variable by fitting a logistic regression model.
#' It conducts a likelihood ratio (deviance) test, computes effect sizes (odds ratios with confidence intervals), and creates
#' a bar plot showing the proportion of each response level within the specified groups.
#'
#' @param data A data frame containing the variables of interest.
#' @param response_var A character string specifying the binary response variable.
#' @param group_var A character string or a vector of character strings specifying the grouping variable(s).
#' @param na.rm Logical; if \code{TRUE} (default) rows with missing values in \code{response_var} or any \code{group_var} are removed.
#' @param print_plot Logical; if \code{TRUE} (default) the proportion plot is printed.
#' @param ... Additional arguments passed to \code{\link[stats]{glm}}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{deviance_anova}{A table of deviance (likelihood ratio) ANOVA results.}
#'   \item{effect_size}{A data frame of the estimated coefficients (exponentiated to show odds ratios) with confidence intervals.}
#'   \item{model_stats}{A list containing AIC, BIC, log-likelihood, and McFadden's pseudo R-squared for the fitted model.}
#'   \item{count_plot}{A \code{ggplot2} bar plot showing the proportion of the response levels by group.}
#'   \item{fitted_model}{The fitted \code{glm} object.}
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom scales percent
#' @importFrom broom tidy
#' @importFrom rlang sym syms !!! 
#' @importFrom stats anova glm logLik reformulate AIC BIC
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # Example using a single grouping variable:
#' data <- data.frame(
#'   response_var = c(0, 0, 1, 1, 1, 0, 1),
#'   group_var1   = c("A", "A", "B", "B", "C", "C", "A")
#' )
#' result <- anova_bin(data, "response_var", "group_var1")
#'
#' # Example using multiple grouping variables:
#' data$group_var2 <- rep(c("X", "Y"), length.out = nrow(data))
#' result_multi <- anova_bin(data, "response_var", c("group_var1", "group_var2"))
#' }
#'
#' @export
anova_bin <- function(data, response_var, group_var, na.rm = TRUE, print_plot = TRUE, ...) {
  
  # Ensure the pipe operator is available
  if (!exists("%>%")) {
    `%>%` <- magrittr::`%>%`
  }
  
  # --- Input checks ---
  # Check that response_var exists in data
  if (!response_var %in% names(data)) {
    stop(paste("Response variable", response_var, "not found in data."))
  }
  
  # Check that all group_var(s) exist in data
  missing_groups <- setdiff(group_var, names(data))
  if (length(missing_groups) > 0) {
    stop(paste("Grouping variable(s)", paste(missing_groups, collapse = ", "), "not found in data."))
  }
  
  # Optionally remove rows with missing values in required columns
  required_vars <- c(response_var, group_var)
  if (na.rm) {
    n_before <- nrow(data)
    # Use if_all() instead of across() to avoid deprecation warnings
    data <- data %>% dplyr::filter(if_all(all_of(required_vars), ~ !is.na(.)))
    n_after <- nrow(data)
    if(n_after < n_before) {
      message(sprintf("Removed %d rows with missing values in required columns.", n_before - n_after))
    }
  }
  
  # Check that the response variable is binary
  unique_vals <- unique(data[[response_var]])
  if (length(unique_vals) != 2) {
    stop(sprintf("The response variable '%s' must be binary. Found %d unique values.", response_var, length(unique_vals)))
  }
  
  # Ensure group variables are factors (if not already)
  data <- data %>% dplyr::mutate(across(all_of(group_var), ~ if (!is.factor(.)) as.factor(.) else .))
  
  # --- Model Fitting ---
  # Construct the logistic regression formula
  formula <- stats::reformulate(termlabels = group_var, response = response_var)
  
  # Fit the logistic regression model
  logistic_model <- glm(formula, data = data, family = binomial(), ...)
  if (!logistic_model$converged) {
    warning("The logistic regression model did not converge.")
  }
  
  # --- Model Diagnostics ---
  # Perform deviance (likelihood ratio) ANOVA
  deviance_anova <- stats::anova(logistic_model, test = "Chisq")
  
  # Compute effect sizes (odds ratios) with confidence intervals using broom
  effect_size <- broom::tidy(logistic_model, exponentiate = TRUE, conf.int = TRUE)
  
  # Calculate additional model statistics: AIC, BIC, logLik, and McFadden's pseudo R²
  null_formula <- stats::reformulate(termlabels = "1", response = response_var)
  null_model <- glm(null_formula, data = data, family = binomial())
  logLik_full <- as.numeric(stats::logLik(logistic_model))
  logLik_null <- as.numeric(stats::logLik(null_model))
  pseudo_r2 <- 1 - (logLik_full / logLik_null)
  model_stats <- list(
    AIC       = AIC(logistic_model),
    BIC       = BIC(logistic_model),
    logLik    = logLik_full,
    pseudo_r2 = pseudo_r2
  )
  
  # --- Data Preparation for Plotting ---
  # Compute counts and proportions of response levels within each group combination
  prop_data <- data %>%
    dplyr::group_by(across(all_of(c(group_var, response_var)))) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::group_by(across(all_of(group_var))) %>%
    dplyr::mutate(prop = count / sum(count)) %>%
    dplyr::ungroup()
  
  # For plotting: if more than one grouping variable is provided, combine them into one factor.
  if (length(group_var) == 1) {
    plot_data <- prop_data
    x_var <- rlang::sym(group_var)
    x_label <- group_var
  } else {
    combined_name <- "group_combined"
    plot_data <- prop_data %>%
      dplyr::mutate(!!rlang::sym(combined_name) := interaction(!!!rlang::syms(group_var), sep = " : "))
    x_var <- rlang::sym(combined_name)
    x_label <- "Combined Group"
  }
  
  # Ensure the response variable is a factor (for proper legend ordering)
  plot_data[[response_var]] <- as.factor(plot_data[[response_var]])
  
  # --- Plotting ---
  count_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = !!x_var, y = prop, fill = !!rlang::sym(response_var))) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::geom_text(ggplot2::aes(label = scales::percent(prop)),
                       position = ggplot2::position_fill(vjust = 0.5), size = 4) +
    ggplot2::labs(x = x_label,
                  y = "Proportion",
                  fill = response_var,
                  title = paste("Proportion of", response_var, "by", x_label)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title   = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.text  = ggplot2::element_text(size = 12),
      legend.position = "bottom"
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent)
  
  if (print_plot) {
    print(count_plot)
  }
  
  # --- Return Results ---
  results <- list(
    deviance_anova = deviance_anova,
    effect_size    = effect_size,
    model_stats    = model_stats,
    count_plot     = count_plot,
    fitted_model   = logistic_model
  )
  
  return(results)
}
