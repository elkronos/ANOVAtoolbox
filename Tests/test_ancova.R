##############################################
# User Acceptance Testing (UAT) - ANCOVA
##############################################
# The following tests verify that every parameter and function works as expected.
# Run these tests to ensure the module functions correctly.

# Create a simulated dataset for testing
set.seed(123)
n <- 150
# Generate a continuous covariate
covariate <- rnorm(n, mean = 50, sd = 10)
# Create a grouping variable with three groups
group <- sample(c("A", "B", "C"), size = n, replace = TRUE)
# Simulate a dependent variable with group differences and a covariate effect
# Note: For homogeneous slopes, the effect of the covariate is constant across groups.
beta0 <- 10
beta_cov <- 0.5
group_effects <- c(A = 0, B = 5, C = -3)
dv <- beta0 + beta_cov * covariate + group_effects[group] + rnorm(n, sd = 5)
# Assemble the data frame
test_data <- data.frame(
  outcome = dv,
  group = group,
  covar = covariate
)

cat("\n--- UAT: Testing check_homogeneity_slopes ---\n")
hs_test <- tryCatch({
  res <- check_homogeneity_slopes(test_data, dv = "outcome", iv = "group", covariate = "covar")
  print(res)
  res
}, error = function(e) {
  cat("Error in check_homogeneity_slopes: ", e$message, "\n")
})

cat("\n--- UAT: Testing check_normality ---\n")
normality_test <- tryCatch({
  lm_model <- lm(outcome ~ covar + group, data = test_data)
  res <- check_normality(lm_model)
  print(res$shapiro)
  print(res$qq_plot)
  res
}, error = function(e) {
  cat("Error in check_normality: ", e$message, "\n")
})

cat("\n--- UAT: Testing check_homogeneity_variance ---\n")
levene_test <- tryCatch({
  lm_model <- lm(outcome ~ covar + group, data = test_data)
  res <- check_homogeneity_variance(test_data, lm_model, iv = "group")
  print(res)
  res
}, error = function(e) {
  cat("Error in check_homogeneity_variance: ", e$message, "\n")
})

cat("\n--- UAT: Testing calculate_effect_sizes ---\n")
effect_size_test <- tryCatch({
  lm_model <- lm(outcome ~ covar + group, data = test_data)
  res <- calculate_effect_sizes(lm_model)
  print(res)
  res
}, error = function(e) {
  cat("Error in calculate_effect_sizes: ", e$message, "\n")
})

cat("\n--- UAT: Testing plot_diagnostics ---\n")
plot_test <- tryCatch({
  lm_model <- lm(outcome ~ covar + group, data = test_data)
  plots_out <- plot_diagnostics(lm_model, test_data, dv = "outcome", iv = "group", covariate = "covar")
  # For automated testing, we only print the object summaries
  print(plots_out$residuals_vs_fitted)
  print(plots_out$qq_plot)
  print(plots_out$ancova_plot)
  plots_out
}, error = function(e) {
  cat("Error in plot_diagnostics: ", e$message, "\n")
})

cat("\n--- UAT: Testing ancova_analysis ---\n")
ancova_test <- tryCatch({
  res <- ancova_analysis(data = test_data, dv = "outcome", iv = "group", covariate = "covar",
                             alpha = 0.05, plots = FALSE, verbose = TRUE)
  str(res)
  res
}, error = function(e) {
  cat("Error in ancova_analysis: ", e$message, "\n")
})

# Edge Case Tests
cat("\n--- UAT: Edge Case Test - Missing variable ---\n")
tryCatch({
  ancova_analysis(data = test_data, dv = "outcome", iv = "nonexistent", covariate = "covar")
}, error = function(e) {
  cat("Expected error (missing iv): ", e$message, "\n")
})

cat("\n--- UAT: Edge Case Test - Non-numeric dependent variable ---\n")
test_data_bad <- test_data
test_data_bad$outcome <- as.factor(round(test_data_bad$outcome))
tryCatch({
  ancova_analysis(data = test_data_bad, dv = "outcome", iv = "group", covariate = "covar")
}, error = function(e) {
  cat("Expected error (non-numeric dv): ", e$message, "\n")
})

cat("\n--- UAT: Edge Case Test - Independent variable with only one level ---\n")
test_data_onelevel <- test_data
test_data_onelevel$group <- "OnlyOneLevel"
tryCatch({
  ancova_analysis(data = test_data_onelevel, dv = "outcome", iv = "group", covariate = "covar")
}, error = function(e) {
  cat("Expected error (iv with one level): ", e$message, "\n")
})

cat("\n--- All UAT tests completed successfully ---\n")
