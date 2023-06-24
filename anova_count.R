library(data.table)
library(stats)
library(ggplot2)
library(broom)
library(rlang)
library(multcomp)

#' Perform Deviance ANOVA
#'
#' This function performs a deviance analysis of variance (ANOVA) for count data using Poisson regression.
#'
#' @param data A data frame containing the variables of interest.
#' @param response_var The name of the response variable (count data).
#' @param group_vars_list A character vector specifying the names of the grouping variables.
#'
#' @return A list containing the following elements:
#'   \item{overdispersion_check_message}{A message indicating whether overdispersion is detected.}
#'   \item{deviance_anova}{The deviance ANOVA table.}
#'   \item{posthoc}{The post-hoc test results.}
#'   \item{effect_size}{The effect size (incidence rate ratios) and confidence intervals.}
#'   \item{count_plot}{A bar plot of the response by group.}
#'
#' @import data.table
#' @import stats
#' @import ggplot2
#' @importFrom broom tidy
#' @import rlang
#' @importFrom multcomp glht mcp
#'
#' @examples
#' \dontrun{
#' # Generate an example data set for a 2x2 factorial design
#' set.seed(1)
#' data <- data.frame(
#'   group1 = factor(rep(c("A", "B"), each = 100)),
#'   group2 = factor(rep(c("X", "Y"), times = 100)),
#'   count = c(rpois(100, 5), rpois(100, 10), rpois(100, 15), rpois(100, 20))
#' )
#'
#' # Run the analysis
#' results <- anova_count(data, "count", c("group1", "group2"))
#'
#' # Check the post-hoc results
#' summary(results$posthoc)
#' }
#'
anova_count <- function(data, response_var, group_vars_list) {
  # Convert data frame to data table
  data <- as.data.table(data)
  
  # Input checks
  for (group_var in group_vars_list) {
    if (!is.factor(data[[group_var]])) {
      stop(paste(group_var, "must be a factor or convertible to a factor"))
    }
  }
  
  if (!all(data[[response_var]] >= 0) || !all(round(data[[response_var]]) == data[[response_var]])) {
    stop(paste(response_var, "must be non-negative integers (count data)"))
  }
  
  # Convert grouping variables to factors if not already
  for (group_var in group_vars_list) {
    data[[group_var]] <- as.factor(data[[group_var]])
  }
  
  # Create interaction term
  interaction_term <- do.call(interaction, data[, ..group_vars_list])
  set(data, j = "interaction_term", value = interaction_term)
  
  # Perform Poisson regression
  poisson_regression <- glm(as.formula(paste(response_var, "~", "interaction_term")), data = data, family = poisson())
  
  # Check for overdispersion (variance significantly greater than mean)
  overdispersion_check <- sum(resid(poisson_regression, type = "pearson")^2) / df.residual(poisson_regression)
  overdispersion_check_message <- ifelse(overdispersion_check > 1.5, "Warning: Overdispersion detected. Consider using Negative Binomial regression.", "No overdispersion detected.")
  
  # Perform deviance ANOVA
  deviance_anova <- anova(poisson_regression, test = "Chisq")
  
  # Perform post-hoc test
  posthoc <- glht(poisson_regression, linfct = mcp(interaction_term = "Tukey"))
  
  # Compute effect size (incidence rate ratios) and confidence intervals
  effect_size <- broom::tidy(poisson_regression, exponentiate = TRUE, conf.int = TRUE)
  
  # Generate a bar plot of the response by group
  count_plot <- ggplot(data, aes(x = interaction_term, y = !!rlang::sym(response_var), fill = interaction_term)) +
    geom_bar(stat = "identity", color = NA, show.legend = FALSE) +
    labs(x = paste(group_vars_list, collapse = " * "), y = response_var, title = paste("Sum of", response_var, "by", paste(group_vars_list, collapse = " * "))) +
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
    overdispersion_check_message = overdispersion_check_message,
    deviance_anova = deviance_anova,
    posthoc = posthoc,
    effect_size = effect_size,
    count_plot = count_plot
  )
  
  return(results)
}