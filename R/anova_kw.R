#' Perform Kruskal-Wallis Test and Post-Hoc Analysis with Boxplot Visualization
#'
#' This function conducts a Kruskal-Wallis test and performs post-hoc analysis for comparing the distribution 
#' of a numeric response variable across different groups defined by one or more categorical variables. It also 
#' generates a boxplot to visualize the distribution of the response variable by group, and (optionally) a QQ plot 
#' of model residuals.
#'
#' @param data A \code{data.frame} or \code{data.table} containing the data.
#' @param response_var A character string indicating the name of the numeric response variable.
#' @param group_vars_vec A character vector specifying the names of the group variable columns.
#' @param plot_qq Logical; if \code{TRUE}, a QQ plot of the model residuals will be generated. Default is \code{FALSE}.
#' @param posthoc_method A character string specifying the p-value adjustment method for Dunn's test. Default is \code{"bh"} (Benjamini-Hochberg).
#'
#' @return A list containing:
#' \describe{
#'   \item{assumptions}{A list with the model residuals, the Anderson-Darling normality test result, and (if requested) the QQ plot.}
#'   \item{kruskal_test}{The result of the Kruskal-Wallis test.}
#'   \item{posthoc}{The result of Dunn's test post-hoc analysis (if applicable).}
#'   \item{boxplot}{The \code{ggplot2} boxplot object.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates and prepares the data by checking that the response variable exists and is numeric, and that the grouping variables exist and are factors.
#'   \item Converts the input data to a \code{data.table} (if not already) and creates an interaction term combining the grouping variables.
#'   \item Fits an ANOVA model (\code{aov}) using the interaction term to obtain residuals, and then applies the Anderson-Darling test (from the \code{nortest} package) on these residuals.
#'   \item Optionally generates a QQ plot of the residuals.
#'   \item Performs a Kruskal-Wallis test on the response variable across the interaction groups.
#'   \item If there are at least two groups, conducts Dunn's post-hoc test (using the \code{dunn.test} package) with the specified p-value adjustment method.
#'   \item Creates a boxplot of the response variable by group. If multiple grouping variables are provided, additional groups are displayed via faceting.
#' }
#'
#' @examples
#' # Single-factor example:
#' set.seed(123)
#' group <- rep(c("Group A", "Group B", "Group C"), each = 50)
#' value <- c(rnorm(50, mean = 5, sd = 1),
#'            rnorm(50, mean = 7, sd = 1.5),
#'            rnorm(50, mean = 4, sd = 0.8))
#' data_single <- data.frame(group, value)
#' data_single$group <- as.factor(data_single$group)
#' results_single <- anova_kw(data_single, "value", c("group"), plot_qq = TRUE)
#'
#' # Two-factor example:
#' group1 <- rep(c("Group A", "Group B"), each = 50)
#' group2 <- rep(c("Group X", "Group Y"), times = 50)
#' value <- c(rnorm(50, mean = 5, sd = 1),
#'            rnorm(50, mean = 7, sd = 1.5),
#'            rnorm(50, mean = 4, sd = 0.8),
#'            rnorm(50, mean = 6, sd = 1.2))
#' data_double <- data.frame(group1, group2, value)
#' data_double$group1 <- as.factor(data_double$group1)
#' data_double$group2 <- as.factor(data_double$group2)
#' results_double <- anova_kw(data_double, "value", c("group1", "group2"), plot_qq = TRUE)
#'
#' @export
#' @import data.table
#' @import ggplot2
#' @import nortest
#' @import rlang
#' @import dunn.test
anova_kw <- function(data, response_var, group_vars_vec, plot_qq = FALSE, posthoc_method = "bh") {
  
  #### Helper Function: Validate and Prepare Data ####
  prepare_data <- function(data, response_var, group_vars_vec) {
    # Check if the response variable exists and is numeric
    if (!(response_var %in% names(data))) 
      stop("response_var not found in data")
    if (!is.numeric(data[[response_var]])) 
      stop("response_var must be numeric")
    
    # Check that all grouping variables exist
    missing_groups <- setdiff(group_vars_vec, names(data))
    if (length(missing_groups) > 0)
      stop("The following group_vars_vec elements are not found in data: ", 
           paste(missing_groups, collapse = ", "))
    
    # Ensure grouping variables are factors
    for (grp in group_vars_vec) {
      if (!is.factor(data[[grp]])) {
        data[[grp]] <- as.factor(data[[grp]])
      }
    }
    
    # Convert to data.table (if not already)
    data <- data.table::as.data.table(data)
    return(data)
  }
  
  #### Helper Function: Create Interaction Term ####
  create_interaction_term <- function(data, group_vars) {
    # Create an interaction term (as characters) to combine group factors.
    data[, interaction_term := do.call(interaction, lapply(.SD, as.character)), .SDcols = group_vars]
    return(data)
  }
  
  #### Helper Function: Compute Residuals & Anderson-Darling Test ####
  compute_residuals_and_adtest <- function(data, response_var) {
    model_formula <- as.formula(paste(response_var, "~ interaction_term"))
    fit <- aov(model_formula, data = data)
    resids <- stats::resid(fit)
    ad_test <- nortest::ad.test(resids)
    list(residuals = resids, ad_test = ad_test)
  }
  
  #### Helper Function: Create QQ Plot using ggplot2 ####
  create_qq_plot <- function(residuals) {
    qq_data <- data.frame(residuals = residuals)
    p <- ggplot2::ggplot(qq_data, ggplot2::aes(sample = residuals)) +
      ggplot2::stat_qq(color = "steelblue") +
      ggplot2::stat_qq_line(color = "red", lwd = 1) +
      ggplot2::labs(title = "QQ Plot of Residuals",
                    x = "Theoretical Quantiles",
                    y = "Sample Quantiles") +
      ggplot2::theme_minimal()
    return(p)
  }
  
  #### Helper Function: Perform Kruskal-Wallis Test ####
  perform_kruskal_test <- function(data, response_var) {
    model_formula <- as.formula(paste(response_var, "~ interaction_term"))
    test_result <- stats::kruskal.test(model_formula, data = data)
    return(test_result)
  }
  
  #### Helper Function: Perform Dunn's Post-Hoc Test ####
  perform_dunn_posthoc <- function(data, response_var, posthoc_method) {
    n_groups <- length(unique(data$interaction_term))
    if (n_groups < 2) {
      warning("Less than 2 groups found. Post-hoc analysis not performed.")
      return(NULL)
    }
    # Run Dunn's test if there are at least 2 groups.
    dunn_res <- dunn.test::dunn.test(x = data[[response_var]], 
                                     g = data$interaction_term, 
                                     method = posthoc_method, 
                                     altp = TRUE)
    return(dunn_res)
  }
  
  #### Helper Function: Create Boxplot Visualization ####
  create_boxplot <- function(data, response_var, group_vars_vec) {
    # Base plot: use the first grouping variable on the x-axis
    p <- ggplot2::ggplot(data, ggplot2::aes_string(x = group_vars_vec[1], 
                                                   y = response_var, 
                                                   fill = group_vars_vec[1])) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_fill_brewer(palette = "Dark2") +
      ggplot2::labs(title = "Boxplot of Response by Group",
                    subtitle = paste("Response Variable:", response_var, 
                                     "| Group Variables:", paste(group_vars_vec, collapse = ", ")),
                    x = group_vars_vec[1],
                    y = response_var) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 15),
                     axis.title.x = ggplot2::element_text(face = "bold", size = 14),
                     axis.title.y = ggplot2::element_text(face = "bold", size = 14),
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
                     axis.text.y = ggplot2::element_text(size = 12))
    
    # If there are additional grouping variables, create a facet variable
    if (length(group_vars_vec) > 1) {
      data[, facet_var := do.call(interaction, .SD), .SDcols = group_vars_vec[-1]]
      p <- p + ggplot2::facet_wrap(~ facet_var, scales = "free_x")
    }
    return(p)
  }
  
  #### Main Execution ####
  
  # 1. Validate and prepare the data
  data <- prepare_data(data, response_var, group_vars_vec)
  
  # 2. Create the interaction term (used in both testing and plotting)
  data <- create_interaction_term(data, group_vars_vec)
  
  # Ensure there are at least 2 unique groups
  if (length(unique(data$interaction_term)) < 2) {
    stop("Not enough groups for analysis. At least 2 unique groups are required.")
  }
  
  # 3. Compute residuals and perform the Anderson-Darling test
  resid_results <- compute_residuals_and_adtest(data, response_var)
  
  # 4. Optionally, generate and display a QQ plot of the residuals
  qqplot_obj <- NULL
  if (plot_qq) {
    qqplot_obj <- create_qq_plot(resid_results$residuals)
    print(qqplot_obj)
  }
  
  # 5. Perform the Kruskal-Wallis test
  kruskal_result <- perform_kruskal_test(data, response_var)
  
  # 6. Perform Dunn's post-hoc test (if there are at least 2 groups)
  posthoc_result <- perform_dunn_posthoc(data, response_var, posthoc_method)
  
  # 7. Generate and display the boxplot visualization
  boxplot_obj <- create_boxplot(data, response_var, group_vars_vec)
  print(boxplot_obj)
  
  # 8. Return all results as a list
  results <- list(
    assumptions = list(
      residuals = resid_results$residuals,
      ad_test = resid_results$ad_test,
      qqplot = qqplot_obj
    ),
    kruskal_test = kruskal_result,
    posthoc = posthoc_result,
    boxplot = boxplot_obj
  )
  
  return(results)
}
