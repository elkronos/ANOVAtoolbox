library(ggplot2)
library(dplyr)
library(scales)

#' Perform One-way "ANOVA" for Binary Response Variable
#'
#' This function performs a one-way "ANOVA" (analysis of variance) for a binary response variable. It fits a logistic regression model with the specified response variable and group variable(s) and performs deviance ANOVA to assess the significance of the group effect. It also computes the effect size (odds ratios) and confidence intervals, and generates a bar plot showing the proportion of the response variable by group.
#'
#' @param data A data frame containing the variables of interest.
#' @param response_var A character string specifying the binary response variable.
#' @param group_var A character string or a vector of character strings specifying the grouping variable(s).
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{deviance_anova}{A table summarizing the deviance ANOVA results.}
#'   \item{effect_size}{A data frame containing the effect size (odds ratios) and confidence intervals for each group.}
#'   \item{count_plot}{A bar plot showing the proportion of the response variable by group.}
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom scales percent
#' @importFrom broom tidy
#' @importFrom rlang sym all_of
#' @importFrom stats anova
#'
#' @examples
#' \dontrun{
#' # Example 1: One-way "ANOVA" for a binary response variable
#' data <- data.frame(response_var = c(0, 0, 1, 1, 1),
#'                    group_var1 = c("A", "A", "B", "B", "C"))
#' result1 <- anova_bin(data, "response_var", "group_var1")
#' }
#'
#' @export
anova_bin <- function(data, response_var, group_var) {
  # Load required packages
  require(ggplot2)
  
  # Check if group_var is a character vector, if not convert it to a character vector
  if (!is.character(group_var)) {
    group_var <- as.character(group_var)
  }
  
  # Convert each grouping variable to a factor if not already
  for (var in group_var) {
    if (!is.factor(data[[var]])) {
      data[[var]] <- as.factor(data[[var]])
    }
  }
  
  # Perform logistic regression
  formula <- as.formula(paste(response_var, "~", paste(group_var, collapse = "+")))
  logistic_regression <- glm(formula, data = data, family = binomial())
  
  # Perform deviance ANOVA
  deviance_anova <- stats::anova(logistic_regression, test = "Chisq")
  
  # Compute effect size (odds ratios) and confidence intervals
  effect_size <- broom::tidy(logistic_regression, exponentiate = TRUE, conf.int = TRUE)
  
  # Compute proportion by group
  prop_data <- data %>%
    group_by(across(all_of(group_var)), across(all_of(response_var))) %>%
    summarize(count = n()) %>%
    group_by(across(all_of(group_var))) %>%
    mutate(prop = count / sum(count))
  
  # Generate a bar plot of the response by group
  count_plot <- ggplot(data = prop_data, aes(x = !!rlang::sym(group_var), y = prop, fill = !!rlang::sym(response_var), group = interaction(!!rlang::sym(group_var), !!rlang::sym(response_var)))) +
    geom_bar(position = "fill", stat = "identity") +
    geom_text(aes(label = scales::percent(prop), group = interaction(!!rlang::sym(group_var), !!rlang::sym(response_var))), position = position_fill(vjust = 0.5), size = 4) +
    labs(x = group_var, y = "Proportion", fill = response_var, 
         title = paste("Proportion of", response_var, "by", group_var)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "bottom"
    )
  
  print(count_plot)
  
  # Return results
  results <- list(
    deviance_anova = deviance_anova,
    effect_size = effect_size,
    count_plot = count_plot
  )
  
  return(results)
}