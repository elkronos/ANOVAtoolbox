###############################################################################
# anova_bin (Logistic "One-Way ANOVA") - Modular & Robust
###############################################################################

#' One-Way "ANOVA" for a Binary Response via Logistic Regression (Modular, Robust)
#'
#' Fits a logistic regression for a binary response and provides:
#' - Type I deviance ANOVA (sequential) or Type II/III via car::Anova()
#' - Optional robust (HC) SEs via sandwich + lmtest
#' - Odds ratios with Wald CIs (uses robust SEs if requested)
#' - Proportions table and stacked bar plot
#'
#' @param data data.frame with variables.
#' @param response_var character(1) name of binary response.
#' @param group_var character() grouping variable(s).
#' @param success_level optional; which response level is the modeled “success”. If NULL,
#'   uses the second factor level (glm default).
#' @param group_ref optional named list of reference levels for grouping factors,
#'   e.g. list(g1="A", g2="X").
#' @param na.rm logical; drop rows with NA in required columns (default TRUE).
#' @param print_plot logical; print the proportions plot (default TRUE).
#' @param include_intercept logical; include intercept row in OR table (default FALSE).
#' @param robust logical; use robust (HC) covariance for coef tests/OR CIs (default FALSE).
#' @param robust_type character; HC type for sandwich::vcovHC (e.g. "HC0","HC1","HC2","HC3").
#' @param anova_type "I","II","III" for Type I (deviance) or Type II/III via car::Anova (default "I").
#' @param anova_test when anova_type is "II"/"III": "LR" or "Wald" (default "LR").
#' @param ... passed to stats::glm()
#'
#' @return list with:
#'   - deviance_anova (for Type I) or car_anova (for Type II/III)
#'   - lr_comparison (null vs full LR)
#'   - coef_tests (coeftest table, robust if requested)
#'   - effect_size (ORs + Wald CIs, robust if requested)
#'   - model_stats (AIC, BIC, logLik, McFadden R2, flags)
#'   - prop_table, count_plot, fitted_model
#' @export
#' @import dplyr ggplot2
#' @importFrom stats glm binomial anova logLik qnorm reformulate AIC BIC
#' @importFrom sandwich vcovHC
#' @importFrom lmtest coeftest
#' @importFrom car Anova
#' @importFrom rlang sym syms
#' @importFrom scales percent
#' @importFrom magrittr %>%
anova_bin <- function(data,
                      response_var,
                      group_var,
                      success_level = NULL,
                      group_ref = NULL,
                      na.rm = TRUE,
                      print_plot = TRUE,
                      include_intercept = FALSE,
                      robust = FALSE,
                      robust_type = "HC0",
                      anova_type = c("I","II","III"),
                      anova_test = c("LR","Wald"),
                      ...) {
  
  if (!exists("%>%")) `%>%` <- magrittr::`%>%`
  anova_type <- match.arg(anova_type)
  anova_test <- match.arg(anova_test)
  
  # 1) Validate & prepare ------------------------------------------------------
  .assert_inputs(data, response_var, group_var)
  required_vars <- c(response_var, group_var)
  if (isTRUE(na.rm)) data <- .drop_na_required(data, required_vars)
  
  resp_out <- .prep_response(data, response_var, success_level)
  data <- resp_out$data
  success_level <- resp_out$success_level
  
  data <- .prep_groups(data, group_var, group_ref)
  
  # 2) Fit model ---------------------------------------------------------------
  model <- .fit_glm(data, response_var, group_var, ...)
  
  # 3) ANOVA tables ------------------------------------------------------------
  deviance_anova <- NULL
  car_anova <- NULL
  if (anova_type == "I") {
    deviance_anova <- stats::anova(model, test = "Chisq")
  } else {
    type_int <- switch(anova_type, "II" = 2, "III" = 3)
    car_anova <- car::Anova(model, type = type_int, test.statistic = anova_test)
  }
  
  # 4) Null vs Full LR comparison ---------------------------------------------
  lr_cmp <- .lr_null_vs_full(model, data, response_var)
  
  # 5) Coefs (robust optional) + Effect sizes ---------------------------------
  vc <- .vcov_for(model, robust = robust, robust_type = robust_type)
  coef_tab <- .coef_tests(model, vc)
  effect_size <- .effect_sizes(coef_tab, include_intercept = include_intercept,
                               conf_level = 0.95, exponentiate = TRUE)
  
  # 6) Model stats -------------------------------------------------------------
  ll_full <- as.numeric(stats::logLik(model))
  null_model <- stats::glm(stats::reformulate("1", response_var), data = data, family = stats::binomial())
  ll_null <- as.numeric(stats::logLik(null_model))
  pseudo_r2 <- 1 - (ll_full / ll_null)
  
  model_stats <- list(
    AIC           = stats::AIC(model),
    BIC           = stats::BIC(model),
    logLik        = ll_full,
    pseudo_r2     = as.numeric(pseudo_r2),
    success_level = success_level,
    robust        = robust,
    robust_type   = if (robust) robust_type else NULL
  )
  
  # 7) Proportions & Plot ------------------------------------------------------
  prop_table <- .prop_table(data, response_var, group_var)
  count_plot <- .plot_props(prop_table, response_var, group_var, success_level)
  if (isTRUE(print_plot)) print(count_plot)
  
  # 8) Return ------------------------------------------------------------------
  list(
    deviance_anova = deviance_anova,
    car_anova      = car_anova,
    lr_comparison  = lr_cmp,
    coef_tests     = coef_tab,
    effect_size    = effect_size,
    model_stats    = model_stats,
    prop_table     = prop_table,
    count_plot     = count_plot,
    fitted_model   = model
  )
}

###############################################################################
# Helpers
###############################################################################

.assert_inputs <- function(data, response_var, group_var) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.character(response_var) || length(response_var) != 1)
    stop("`response_var` must be a single character string.")
  if (missing(group_var)) stop("Please provide at least one grouping variable in `group_var`.")
  if (!is.character(group_var)) stop("`group_var` must be a character vector (or scalar).")
  if (!(response_var %in% names(data)))
    stop(sprintf("Response variable '%s' not found in `data`.", response_var))
  missing_groups <- setdiff(group_var, names(data))
  if (length(missing_groups) > 0)
    stop(paste("Grouping variable(s) not in `data`:", paste(missing_groups, collapse = ", ")))
}

.drop_na_required <- function(data, required_vars) {
  n_before <- nrow(data)
  keep <- stats::complete.cases(data[, required_vars, drop = FALSE])
  data2 <- data[keep, , drop = FALSE]
  n_after <- nrow(data2)
  if (n_after < n_before) {
    message(sprintf("Removed %d rows with missing values in required columns.", n_before - n_after))
  }
  data2
}

.prep_response <- function(data, response_var, success_level) {
  resp <- data[[response_var]]
  if (!is.factor(resp)) {
    if (is.logical(resp)) {
      resp <- factor(resp, levels = c(FALSE, TRUE))
    } else if (is.numeric(resp) && all(na.omit(unique(resp)) %in% c(0,1))) {
      resp <- factor(resp, levels = c(0,1))
    } else {
      resp <- factor(resp)
    }
  }
  lvls <- levels(resp)
  if (length(lvls) != 2) {
    stop(sprintf("The response '%s' must be binary (2 levels). Found %d levels: %s",
                 response_var, length(lvls), paste(lvls, collapse = ", ")))
  }
  if (!is.null(success_level)) {
    if (!(success_level %in% lvls)) {
      stop(sprintf("`success_level` ('%s') not found in response levels: %s",
                   success_level, paste(lvls, collapse = ", ")))
    }
    other <- setdiff(lvls, success_level)
    resp <- factor(resp, levels = c(other, success_level))
  } else {
    success_level <- levels(resp)[2]
  }
  data[[response_var]] <- resp
  list(data = data, success_level = success_level)
}

.prep_groups <- function(data, group_var, group_ref) {
  for (g in group_var) {
    if (!is.factor(data[[g]])) data[[g]] <- factor(data[[g]])
    if (!is.null(group_ref) && !is.null(group_ref[[g]])) {
      ref <- group_ref[[g]]
      if (!(ref %in% levels(data[[g]]))) {
        stop(sprintf("Reference level '%s' not found in '%s'. Levels: %s",
                     ref, g, paste(levels(data[[g]]), collapse = ", ")))
      }
      data[[g]] <- stats::relevel(data[[g]], ref = ref)
    }
  }
  data
}

.fit_glm <- function(data, response_var, group_var, ...) {
  form <- stats::reformulate(termlabels = group_var, response = response_var)
  mod <- stats::glm(form, data = data, family = stats::binomial(), ...)
  if (!isTRUE(mod$converged)) warning("The logistic regression model did not converge.")
  mod
}

.lr_null_vs_full <- function(model, data, response_var) {
  null_form <- stats::reformulate(termlabels = "1", response = response_var)
  null_model <- stats::glm(null_form, data = data, family = stats::binomial())
  ll_full <- as.numeric(stats::logLik(model))
  ll_null <- as.numeric(stats::logLik(null_model))
  df_full <- attr(stats::logLik(model), "df")
  df_null <- attr(stats::logLik(null_model), "df")
  lr_stat <- 2 * (ll_full - ll_null)
  lr_df <- df_full - df_null
  lr_p <- stats::pchisq(lr_stat, df = lr_df, lower.tail = FALSE)
  list(
    models = data.frame(
      Model = c("Null", "Full"),
      logLik = c(ll_null, ll_full),
      df = c(df_null, df_full)
    ),
    test = data.frame(
      Test = "Likelihood ratio (Full vs Null)",
      Chisq = lr_stat,
      Df = lr_df,
      `Pr(>Chisq)` = lr_p,
      check.names = FALSE
    )
  )
}

.vcov_for <- function(model, robust, robust_type) {
  if (robust) sandwich::vcovHC(model, type = robust_type) else stats::vcov(model)
}

# FIXED: resilient to odd shapes/labels from lmtest::coeftest()
.coef_tests <- function(model, vc) {
  ct <- lmtest::coeftest(model, vcov. = vc)
  ct_mat <- as.matrix(ct)
  
  term <- rownames(ct_mat)
  if (is.null(term)) term <- paste0("param_", seq_len(nrow(ct_mat)))
  
  cn_raw  <- colnames(ct_mat)
  cn_norm <- tolower(gsub("[^a-zA-Z0-9]+", "", cn_raw))
  
  get_col <- function(target) {
    opts <- switch(target,
                   estimate = c("estimate","coef","coefficients"),
                   stderr   = c("stderr","stderror","se"),
                   stat     = c("zvalue","tvalue","z","t"),
                   pvalue   = c("prgtzltz","pvalue","p","prgt|z|","prgt|t|")
    )
    idx <- match(opts, cn_norm, nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) > 0L) idx[1L] else NA_integer_
  }
  
  i_est <- get_col("estimate")
  i_se  <- get_col("stderr")
  i_st  <- get_col("stat")
  i_p   <- get_col("pvalue")
  
  if (any(is.na(c(i_est, i_se, i_st, i_p))) && ncol(ct_mat) >= 4) {
    if (is.na(i_est)) i_est <- 1L
    if (is.na(i_se))  i_se  <- 2L
    if (is.na(i_st))  i_st  <- 3L
    if (is.na(i_p))   i_p   <- 4L
  }
  
  safe_get <- function(j) {
    if (is.na(j) || j < 1L || j > ncol(ct_mat)) rep(NA_real_, nrow(ct_mat)) else ct_mat[, j]
  }
  
  out <- data.frame(
    term      = term,
    estimate  = as.numeric(safe_get(i_est)),
    std_error = as.numeric(safe_get(i_se)),
    z_value   = as.numeric(safe_get(i_st)),
    p_value   = as.numeric(safe_get(i_p)),
    row.names = NULL,
    check.names = FALSE
  )
  
  num_cols <- c("estimate","std_error","z_value","p_value")
  for (nc in num_cols) out[[nc]] <- suppressWarnings(as.numeric(out[[nc]]))
  
  out
}

.effect_sizes <- function(coef_tab, include_intercept, conf_level = 0.95, exponentiate = TRUE) {
  z <- stats::qnorm(1 - (1 - conf_level)/2)
  est <- coef_tab$estimate
  se  <- coef_tab$std_error
  lo  <- est - z * se
  hi  <- est + z * se
  
  if (exponentiate) {
    est <- exp(est); lo <- exp(lo); hi <- exp(hi)
    out <- data.frame(
      term = coef_tab$term,
      odds_ratio = est,
      conf.low   = lo,
      conf.high  = hi,
      p_value    = coef_tab$p_value,
      row.names = NULL
    )
  } else {
    out <- data.frame(
      term = coef_tab$term,
      estimate = est,
      conf.low = lo,
      conf.high = hi,
      p_value  = coef_tab$p_value,
      row.names = NULL
    )
  }
  if (!include_intercept) out <- out[out$term != "(Intercept)", , drop = FALSE]
  out
}

.prop_table <- function(data, response_var, group_var) {
  counts <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_var, response_var)))) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop")
  
  totals <- counts %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) %>%
    dplyr::summarise(group_total = sum(count), .groups = "drop")
  
  counts %>%
    dplyr::left_join(totals, by = group_var) %>%
    dplyr::mutate(prop = count / group_total)
}

.plot_props <- function(prop_table, response_var, group_var, success_level) {
  if (length(group_var) == 1) {
    x_var <- rlang::sym(group_var)
    plot_data <- prop_table
    x_label <- group_var
  } else {
    combined_name <- "group_combined"
    plot_data <- prop_table %>%
      dplyr::mutate(!!rlang::sym(combined_name) := interaction(!!!rlang::syms(group_var), sep = " : "))
    x_var <- rlang::sym(combined_name)
    x_label <- "Combined Group"
  }
  plot_data[[response_var]] <- as.factor(plot_data[[response_var]])
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = !!x_var, y = prop, fill = !!rlang::sym(response_var))
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::percent(prop)),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 4
    ) +
    ggplot2::labs(
      x = x_label,
      y = "Proportion",
      fill = response_var,
      title = sprintf("Proportion of '%s' levels by %s (success: %s)",
                      response_var, x_label, success_level)
    ) +
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
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1))
}
