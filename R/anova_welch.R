#' Welch's ANOVA with Assumption Checks, Post-hoc Tests, and Effect Sizes
#'
#' Performs Welch's ANOVA (stats::oneway.test with var.equal = FALSE) for a numeric response
#' across groups defined by one or more grouping variables. Also:
#' - Checks normality of residuals (QQ plot; Anderson-Darling if n >= 8 and {nortest} available)
#' - Computes group-wise summary statistics (mean, SD, N, SE, and CI)
#' - Pairwise Welch t-tests with multiple-comparison correction (if > 2 groups)
#' - Pairwise Cohen's d (if {effsize} available)
#'
#' NOTE: The Anderson-Darling test requires at least 8 residuals and the {nortest} package.
#'
#' @param data data.frame or data.table
#' @param response_var character; name of numeric response column
#' @param group_vars character vector; one or more grouping variables
#' @param conf_level numeric in (0,1); confidence level for CIs (default 0.95)
#' @param posthoc_adjust character; p-value adjustment method for post-hoc tests
#'   (e.g., "bonferroni", "holm", "BH"); default "bonferroni"
#'
#' @return A list with:
#' \describe{
#'   \item{assumptions}{list with \code{qq_plot} (ggplot) and \code{ad_test} (htest or NA)}
#'   \item{welch_anova}{object returned by \code{stats::oneway.test}}
#'   \item{posthoc}{object returned by \code{stats::pairwise.t.test} or NULL if \eqn{<=} 2 groups}
#'   \item{summary_stats}{data.table of group-wise mean, sd, n, se, and CI}
#'   \item{effect_sizes}{data.table of pairwise Cohen's d (or NULL if unavailable)}
#'   \item{means_plot}{ggplot bar chart of means with CIs}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' group <- rep(c("A","B","C"), each = 30)
#' value <- c(rnorm(30, 10, 2), rnorm(30, 12, 2.5), rnorm(30, 9, 1.5))
#' df <- data.frame(group, value)
#' res <- anova_welch(df, "value", "group")
#' print(res$assumptions$qq_plot)
#' print(res$welch_anova)
#' res$summary_stats
#' print(res$means_plot)
#' }
#'
#' @export
anova_welch <- function(data,
                        response_var,
                        group_vars,
                        conf_level = 0.95,
                        posthoc_adjust = "bonferroni") {
  
  ## -------- Input validation --------
  if (!is.data.frame(data)) stop("`data` must be a data.frame or data.table.")
  if (!is.character(response_var) || length(response_var) != 1)
    stop("`response_var` must be a single column name (character).")
  if (!all(group_vars %in% names(data))) {
    missing_vars <- group_vars[!group_vars %in% names(data)]
    stop("Grouping variable(s) not found in data: ", paste(missing_vars, collapse = ", "))
  }
  if (!(response_var %in% names(data)))
    stop("Response variable `", response_var, "` not found in data.")
  if (!is.numeric(data[[response_var]]))
    stop("Response variable `", response_var, "` must be numeric.")
  if (!is.numeric(conf_level) || length(conf_level) != 1 || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be a single number in (0,1).")
  
  # Work on a local copy; drop incomplete rows for involved columns
  cols_needed <- c(response_var, group_vars)
  d <- data[stats::complete.cases(data[, cols_needed, drop = FALSE]), , drop = FALSE]
  
  # Require at least 2 groups with >=1 obs each
  if (nrow(d) < 2L) stop("Not enough complete observations to proceed.")
  # Build combined grouping factor
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package {data.table} is required for this function. Please install it.")
  }
  data.table::setDT(d)
  d[, group_interaction := do.call(interaction, c(.SD, sep = " : ")), .SDcols = group_vars]
  d[, group_interaction := base::as.factor(group_interaction)]
  
  # Ensure there are at least 2 distinct groups
  n_groups <- length(base::levels(d$group_interaction))
  if (n_groups < 2L)
    stop("At least two distinct groups are required.")
  
  ## -------- Model & residuals for assumption checks --------
  # Residuals from a simple fixed-effects model (normality check)
  model_formula <- stats::as.formula(paste(response_var, "~ group_interaction"))
  lm_fit <- stats::lm(model_formula, data = d)
  resid_vec <- stats::residuals(lm_fit)
  
  # QQ plot (ggplot2)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package {ggplot2} is required for plots. Please install it.")
  }
  qq_df <- data.frame(residuals = resid_vec)
  qq_plot <- ggplot2::ggplot(qq_df, ggplot2::aes(sample = residuals)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line() +
    ggplot2::labs(title = "QQ Plot of Model Residuals",
                  x = "Theoretical Quantiles",
                  y = "Sample Quantiles") +
    ggplot2::theme_minimal()
  
  # Anderson–Darling test (optional; requires {nortest} and n >= 8)
  ad_test_result <- NA
  if (length(resid_vec) >= 8L) {
    if (requireNamespace("nortest", quietly = TRUE)) {
      ad_test_result <- tryCatch(
        nortest::ad.test(resid_vec),
        error = function(e) {
          warning("Anderson-Darling test failed: ", conditionMessage(e))
          NA
        }
      )
    } else {
      message("Package {nortest} not installed; skipping Anderson-Darling normality test.")
    }
  } else {
    message("Fewer than 8 residuals; skipping Anderson-Darling normality test.")
  }
  
  ## -------- Welch's ANOVA (correct method) --------
  welch_anova <- stats::oneway.test(
    stats::as.formula(paste(response_var, "~ group_interaction")),
    data = d,
    var.equal = FALSE
  )
  
  ## -------- Post-hoc pairwise Welch t-tests (if > 2 groups) --------
  posthoc_result <- NULL
  if (n_groups > 2L) {
    posthoc_result <- stats::pairwise.t.test(
      x = d[[response_var]],
      g = d[["group_interaction"]],
      p.adjust.method = posthoc_adjust,
      pool.sd = FALSE,           # Welch (unequal variances)
      alternative = "two.sided"
    )
  }
  
  ## -------- Group-wise summary statistics --------
  alpha <- 1 - conf_level
  crit_fun <- function(df) stats::qt(1 - alpha/2, df)
  summary_stats <- d[, .(
    mean = base::mean(get(response_var), na.rm = TRUE),
    sd   = stats::sd(get(response_var), na.rm = TRUE),
    n    = .N
  ), by = group_interaction]
  summary_stats[, se := sd / sqrt(pmax(n, 1))]
  summary_stats[, tcrit := crit_fun(pmax(n - 1, 1))]
  summary_stats[, `:=`(
    ci_low  = mean - tcrit * se,
    ci_high = mean + tcrit * se
  )]
  # Drop helper
  summary_stats[, tcrit := NULL]
  
  ## -------- Pairwise effect sizes (Cohen's d) --------
  effect_sizes <- NULL
  if (n_groups >= 2L) {
    if (requireNamespace("effsize", quietly = TRUE)) {
      combs <- utils::combn(base::levels(d$group_interaction), 2, simplify = FALSE)
      es_list <- lapply(combs, function(pair) {
        x <- d[group_interaction == pair[1], get(response_var)]
        y <- d[group_interaction == pair[2], get(response_var)]
        if (length(x) < 2L || length(y) < 2L) return(NULL)
        est <- tryCatch(
          effsize::cohen.d(x, y, hedges.correction = TRUE, na.rm = TRUE),
          error = function(e) NULL
        )
        if (is.null(est)) return(NULL)
        d_val <- unname(est$estimate)
        mag <- if (is.finite(d_val)) {
          absd <- abs(d_val)
          if (absd < 0.2) "negligible"
          else if (absd < 0.5) "small"
          else if (absd < 0.8) "medium"
          else "large"
        } else NA_character_
        data.table::data.table(
          group1 = pair[1], group2 = pair[2],
          cohen_d = d_val, magnitude = mag
        )
      })
      es_dt <- data.table::rbindlist(es_list, fill = TRUE)
      if (nrow(es_dt)) effect_sizes <- es_dt
    } else {
      message("Package {effsize} not installed; skipping Cohen's d effect sizes.")
    }
  }
  
  ## -------- Plot: Group means with 95% CI --------
  means_plot <- ggplot2::ggplot(summary_stats,
                                ggplot2::aes(x = group_interaction, y = mean, fill = group_interaction)) +
    ggplot2::geom_col(color = "black", position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high),
                           width = 0.2, position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::geom_text(ggplot2::aes(label = round(mean, 2)),
                       vjust = -0.5, position = ggplot2::position_dodge(width = 0.9), size = 3) +
    ggplot2::labs(x = "Group", y = "Mean", title = paste0("Group Means with ", round(conf_level*100), "% CIs")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  ## -------- Return --------
  list(
    assumptions   = list(qq_plot = qq_plot, ad_test = ad_test_result),
    welch_anova   = welch_anova,
    posthoc       = posthoc_result,
    summary_stats = summary_stats[],
    effect_sizes  = effect_sizes,
    means_plot    = means_plot
  )
}
