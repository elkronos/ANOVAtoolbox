#' Perform Repeated Measures ANOVA with Diagnostics and Effect Size Calculation
#'
#' This function fits a repeated measures analysis of variance (ANOVA) model to data in long format,
#' using one or more within-subject (and optionally between-subject) factors. It uses \code{afex::aov_ez}
#' for model fitting, assesses model assumptions (normality and sphericity), computes effect sizes (partial eta²),
#' and generates various diagnostic plots (residuals, Q-Q, and LM diagnostic plots). Optionally, if an estimated marginal
#' means specification is provided, the function computes and plots estimated marginal means.
#'
#' @param data A \code{data.frame} containing the dataset in long format.
#' @param subject A character string specifying the name of the subject identifier column.
#' @param dv A character string specifying the name of the dependent variable column.
#' @param within A character vector specifying the names of the within-subject factor(s).
#' @param between (Optional) A character vector specifying the names of the between-subject factor(s).
#' @param factorize Logical. If \code{TRUE} (default), the function converts the subject and grouping columns to factors,
#'   and the dependent variable to numeric.
#' @param emm_specs (Optional) Specification for estimated marginal means (e.g., a factor name or formula). If provided,
#'   an estimated marginal means plot is generated.
#' @param diagnostics Logical. If \code{TRUE} (default), standard diagnostic plots of the underlying linear model are generated.
#' @param verbose Logical. If \code{TRUE} (default), progress messages are printed to the console.
#' @param ... Additional arguments passed to \code{afex::aov_ez}.
#'
#' @details
#' The function carries out the following steps:
#' \enumerate{
#'   \item \strong{Input Checks and Data Preparation:} Verifies that the specified subject, dependent variable,
#'     and factor columns exist in \code{data}. Optionally converts these columns to the proper data types.
#'   \item \strong{Model Fitting:} Fits a repeated measures ANOVA using \code{afex::aov_ez}.
#'   \item \strong{Normality Diagnostics:} Computes residuals from the underlying linear model, performs a Shapiro–Wilk
#'         test for normality, and generates a residuals vs index plot and a normal Q-Q plot.
#'   \item \strong{Sphericity Diagnostics:} If applicable (multiple within factors or one factor with more than two levels),
#'         performs Mauchly's test for sphericity using \code{car::Anova}.
#'   \item \strong{Effect Size Calculation:} Computes partial eta² for the ANOVA using the \code{effectsize} package.
#'   \item \strong{Additional Diagnostics:} Optionally produces standard diagnostic plots of the underlying linear model.
#'   \item \strong{Estimated Marginal Means (EMMs):} If \code{emm_specs} is provided, computes and plots estimated marginal means
#'         using \code{emmeans} and \code{ggplot2}.
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{anova_obj}}{The fitted repeated measures ANOVA object from \code{afex::aov_ez} (with the processed data attached).}
#'   \item{\code{normality}}{A list with the residuals and Shapiro–Wilk test results for normality.}
#'   \item{\code{sphericity}}{The results of Mauchly's test for sphericity (if applicable), otherwise \code{NULL}.}
#'   \item{\code{effect_sizes}}{A table of effect sizes (partial eta²).}
#'   \item{\code{diagnostic_plots}}{A recorded plot object (using \code{recordPlot}) for the diagnostic plots of the LM, if generated.}
#'   \item{\code{emm_plot}}{A \code{ggplot2} object displaying the estimated marginal means (if \code{emm_specs} is provided), otherwise \code{NULL}.}
#' }
#'
#' @seealso \code{\link[afex]{aov_ez}}, \code{\link[emmeans]{emmeans}}, \code{\link[car]{Anova}},
#'   \code{\link[effectsize]{eta_squared}}, \code{\link[ggplot2]{ggplot}}
#'
#' @examples
#' \dontrun{
#'   # Assume 'my_data' is a data.frame in long format with columns:
#'   # "ID" (subject identifier), "Score" (dependent variable), "Time" (within-subject factor),
#'   # and "Group" (between-subject factor).
#'
#'   results <- anova_rm(data = my_data,
#'                       subject = "ID",
#'                       dv = "Score",
#'                       within = c("Time"),
#'                       between = c("Group"),
#'                       factorize = TRUE,
#'                       emm_specs = "Time",
#'                       diagnostics = TRUE,
#'                       verbose = TRUE)
#'
#'   # Display the fitted ANOVA model
#'   print(results$anova_obj)
#'
#'   # View the results of the Shapiro–Wilk normality test
#'   print(results$normality$shapiro_test)
#'
#'   # If generated, display the estimated marginal means plot
#'   if (!is.null(results$emm_plot)) {
#'     print(results$emm_plot)
#'   }
#' }
#'
#' @export
#' @import afex
#' @import emmeans
#' @import car
#' @import effectsize
#' @import ggplot2
anova_rm <- function(data, subject, dv, within, between = NULL,
                     factorize = TRUE, emm_specs = NULL,
                     diagnostics = TRUE, verbose = TRUE, ...) {
  # --- Input Checks and Data Preparation ---
  req_cols <- c(subject, dv, within)
  if (!is.null(between)) {
    req_cols <- c(req_cols, between)
  }
  missing_cols <- setdiff(req_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing in data: ",
         paste(missing_cols, collapse = ", "))
  }
  
  if (factorize) {
    data[[subject]] <- as.factor(data[[subject]])
    data[[dv]] <- as.numeric(data[[dv]])
    data[within] <- lapply(data[within], as.factor)
    if (!is.null(between)) {
      data[between] <- lapply(data[between], as.factor)
    }
  }
  
  # --- Model Fitting ---
  if (verbose) cat("Fitting repeated measures ANOVA model...\n")
  anova_obj <- afex::aov_ez(id = subject,
                            dv = dv,
                            data = data,
                            within = within,
                            between = between,
                            ...)
  if (verbose) print(anova_obj)
  
  # Override the 'data' element in the returned object with our processed data.
  anova_obj$data <- data
  
  # --- Normality Diagnostics ---
  if (verbose) cat("Performing normality diagnostics...\n")
  if (!is.null(anova_obj$lm)) {
    resid_values <- residuals(anova_obj$lm)
    shapiro_res <- shapiro.test(resid_values)
    
    op <- par(mfrow = c(1, 2))
    plot(resid_values, main = "Residuals vs Index",
         xlab = "Observation Index", ylab = "Residuals")
    qqnorm(resid_values, main = "Normal Q-Q Plot")
    qqline(resid_values)
    par(op)
  } else {
    warning("No underlying linear model found; skipping normality diagnostics.")
    resid_values <- NULL
    shapiro_res <- NULL
  }
  
  # --- Sphericity Diagnostics ---
  check_sph <- FALSE
  if (length(within) > 1) {
    check_sph <- TRUE
  } else if (length(within) == 1) {
    if (nlevels(as.factor(data[[within]])) > 2) check_sph <- TRUE
  }
  
  if (check_sph) {
    if (verbose) cat("Performing sphericity diagnostics...\n")
    data[within] <- lapply(data[within], as.factor)
    idata <- unique(data[within])
    within_formula <- as.formula(paste("~", paste(within, collapse = "*")))
    fit <- lm(as.formula(paste(dv, "~ 1")), data = data)
    sphericity <- tryCatch({
      car::Anova(fit, idata = idata, idesign = within_formula, type = "III")
    }, error = function(e) {
      warning("Sphericity test failed: ", e$message)
      NULL
    })
  } else {
    if (verbose) cat("Sphericity test not applicable (insufficient levels in within-subject factor).\n")
    sphericity <- NULL
  }
  
  # --- Effect Size Calculation ---
  if (verbose) cat("Calculating effect sizes...\n")
  effect_sizes <- effectsize::eta_squared(anova_obj)
  
  # --- Additional Diagnostic Plots ---
  diagnostic_plots <- NULL
  if (diagnostics && !is.null(anova_obj$lm)) {
    diagnostic_plots <- tryCatch({
      if (inherits(anova_obj$lm, "lm") && !inherits(anova_obj$lm, "mlm")) {
        op <- par(mfrow = c(2, 2))
        plot(anova_obj$lm)
        par(op)
        recordPlot()  # Capture the base R diagnostic plots
      } else {
        # For non-univariate models, simply return NULL without warning.
        NULL
      }
    }, error = function(e) {
      warning("Diagnostic plot generation failed: ", e$message)
      NULL
    })
  }
  
  # --- Estimated Marginal Means Plot ---
  emm_plot <- NULL
  if (!is.null(emm_specs)) {
    if (verbose) cat("Computing and plotting estimated marginal means...\n")
    emm_obj <- emmeans::emmeans(anova_obj, specs = emm_specs)
    emm_df <- as.data.frame(emm_obj)
    # Use the first column as the x-axis variable for the plot.
    xvar <- names(emm_df)[1]
    emm_plot <- ggplot2::ggplot(emm_df, ggplot2::aes_string(x = xvar, y = "emmean", group = 1)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_line() +
      ggplot2::labs(x = xvar, y = "Estimated Marginal Mean") +
      ggplot2::theme_minimal()
    print(emm_plot)
  }
  
  # --- Return Results ---
  return(list(
    anova_obj = anova_obj,
    normality = list(residuals = resid_values, shapiro_test = shapiro_res),
    sphericity = sphericity,
    effect_sizes = effect_sizes,
    diagnostic_plots = diagnostic_plots,
    emm_plot = emm_plot
  ))
}
