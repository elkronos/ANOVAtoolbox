#' Perform ANOVA for a Generalized Linear Model with Grouping Variables
#'
#' Fits a GLM, performs Type II/III ANOVA via car::Anova(), optional Tukey/emmeans
#' post-hoc, diagnostics, and plots. Modular helpers are defined inside.
#'
#' @param data data.frame or data.table
#' @param response_var character, name of response column
#' @param group_vars_vec character vector, names of grouping variables (>=1)
#' @param family a \code{family} object for GLM; default \code{gaussian()}
#' @param sum_squares_type "II" (default) or "III"
#' @param plot_residuals logical, make residuals vs fitted plot
#' @param full_factorial logical; if TRUE and >1 group vars, use full factorial (\code{*}),
#'   else collapse groups into a single interaction factor
#' @param posthoc_method one of "none","multcomp","emmeans" (default "multcomp")
#' @param emmeans_simple_for if using emmeans & full factorial with >=2 factors,
#'   choose which factor to compare within the second factor (default = first name)
#'
#' @return list(model, model_stats, anova_test, posthoc, effect_sizes,
#'   confidence_intervals, residuals, residuals_plot, boxplot, notes)
#'
#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom car Anova
#' @importFrom multcomp glht mcp
#' @importFrom stats glm as.formula resid fitted AIC BIC confint vcov coef residuals
anova_glm <- function(
    data, response_var, group_vars_vec,
    family = gaussian(), sum_squares_type = "II",
    plot_residuals = FALSE, full_factorial = FALSE,
    posthoc_method = c("multcomp","emmeans","none"),
    emmeans_simple_for = NULL
) {
  
  # ========================= Helpers =========================
  .check_inputs <- function(dt, y, groups) {
    if (!y %in% names(dt)) stop("Response variable '", y, "' not found.")
    missing_groups <- setdiff(groups, names(dt))
    if (length(missing_groups)) {
      stop("Grouping variable(s) not found: ", paste(missing_groups, collapse = ", "))
    }
    invisible(NULL)
  }
  
  .coerce_groups_to_factor <- function(dt, groups) {
    for (v in groups) {
      if (!is.factor(dt[[v]])) {
        dt[[v]] <- as.factor(dt[[v]])
        message(sprintf("Variable '%s' was converted to a factor.", v))
      }
    }
    dt
  }
  
  .build_formula_and_data <- function(dt, y, groups, full_factorial) {
    if (length(groups) == 1) {
      form <- as.formula(paste(y, "~", groups[1]))
      list(formula = form, data = dt, used_interaction = FALSE)
    } else {
      if (isTRUE(full_factorial)) {
        preds <- paste(groups, collapse = " * ")
        form  <- as.formula(paste(y, "~", preds))
        list(formula = form, data = dt, used_interaction = FALSE)
      } else {
        # Single combined interaction factor
        dt[, interaction_term := do.call(interaction, .SD), .SDcols = groups]
        form <- as.formula(paste(y, "~ interaction_term"))
        list(formula = form, data = dt, used_interaction = TRUE)
      }
    }
  }
  
  .fit_glm <- function(formula, dt, fam) {
    stats::glm(formula, data = dt, family = fam)
  }
  
  .anova_car <- function(model, type, fam) {
    type <- toupper(type)
    if (!type %in% c("II","III")) stop("sum_squares_type must be 'II' or 'III'.")
    test_stat <- if (fam$family == "gaussian") "F" else "LR"
    
    if (type == "III") {
      old_opts <- options(contrasts = c("contr.sum", "contr.poly"))
      on.exit(options(old_opts), add = TRUE)
    }
    car::Anova(model, type = if (type == "II") 2L else 3L, test.statistic = test_stat)
  }
  
  # Compute effect sizes (R^2 from model$y for Gaussian; McFadden pseudo-R^2 for GLMs)
  .effect_sizes <- function(model, fam, ycol, dt) {
    pseudo <- if (isTRUE(is.finite(model$null.deviance)) && model$null.deviance != 0) {
      1 - (model$deviance / model$null.deviance)
    } else {
      NA_real_
    }
    
    r2_gaussian <- NA_real_
    if (fam$family == "gaussian") {
      y_used  <- tryCatch(as.numeric(model$y), error = function(e) NULL)
      res_used <- tryCatch(residuals(model, type = "response"), error = function(e) NULL)
      if (!is.null(y_used) && !is.null(res_used)) {
        n <- min(length(y_used), length(res_used))
        if (n > 1) {
          y_used <- y_used[seq_len(n)]
          res_used <- res_used[seq_len(n)]
          rss <- sum(res_used^2)
          tss <- sum((y_used - mean(y_used))^2)
          if (is.finite(tss) && tss > 0) r2_gaussian <- 1 - rss / tss
        }
      }
    }
    
    list(
      pseudo_r2_mcfadden = pseudo,
      r2_gaussian = r2_gaussian
    )
  }
  
  .overdispersion_note <- function(model, fam) {
    if (fam$family %in% c("poisson", "binomial")) {
      disp <- model$deviance / model$df.residual
      if (is.finite(disp) && disp > 1.5) {
        return(sprintf(
          "Possible overdispersion detected (dispersion ≈ %.2f). Consider quasi-%s or a different variance model (e.g., negative binomial).",
          disp, fam$family
        ))
      }
    }
    NULL
  }
  
  .residuals_plot <- function(model, make_plot) {
    if (!isTRUE(make_plot)) return(NULL)
    df <- data.frame(
      Fitted = fitted(model),
      Residuals = residuals(model, type = "deviance")
    )
    ggplot(df, aes(Fitted, Residuals)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "loess", se = FALSE, linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Deviance Residuals vs Fitted", x = "Fitted values", y = "Deviance residuals") +
      theme_minimal()
  }
  
  .boxplot <- function(dt, y, groups, used_interaction, full_factorial, fam) {
    if (length(groups) == 1) {
      return(
        ggplot(dt, aes_string(x = groups[1], y = y, fill = groups[1])) +
          geom_boxplot(show.legend = FALSE) +
          labs(title = "Boxplot of Response by Group",
               subtitle = sprintf("Response: %s | Group: %s", y, groups[1]),
               x = groups[1], y = y) +
          theme_minimal()
      )
    }
    
    if (isTRUE(full_factorial)) {
      facet_formula <- as.formula(paste("~", paste(groups[-1], collapse = " + ")))
      return(
        ggplot(dt, aes_string(x = groups[1], y = y, fill = groups[1])) +
          geom_boxplot(show.legend = FALSE) +
          labs(title = "Boxplot of Response by Group",
               subtitle = sprintf("Response: %s | Groups: %s", y, paste(groups, collapse = ", ")),
               x = groups[1], y = y) +
          theme_minimal() +
          facet_grid(facet_formula)
      )
    }
    
    if (isTRUE(used_interaction)) {
      return(
        ggplot(dt, aes(x = interaction_term, y = .data[[y]], fill = interaction_term)) +
          geom_boxplot(show.legend = FALSE) +
          labs(title = "Boxplot by Combined Group (Interaction)",
               subtitle = sprintf("Response: %s | Combined groups: %s", y, paste(groups, collapse = ", ")),
               x = "Combined group", y = y) +
          theme_minimal()
      )
    }
    
    NULL
  }
  
  .confint_safe <- function(model) {
    tryCatch({
      ci <- as.data.frame(confint(model))
      colnames(ci) <- c("2.5 %", "97.5 %")
      ci
    }, error = function(e) {
      se <- sqrt(diag(vcov(model)))
      co <- coef(model)
      ci <- cbind(co - 1.96 * se, co + 1.96 * se)
      colnames(ci) <- c("2.5 %", "97.5 %")
      as.data.frame(ci)
    })
  }
  
  .posthoc_multcomp <- function(model, dt, groups, used_interaction, y, fam) {
    # Returns summary(glht) or NULL
    if (length(groups) == 1) {
      lv <- levels(dt[[groups[1]]]); if (length(lv) < 3) return(NULL)
      m2 <- glm(as.formula(paste(y, "~", groups[1], "- 1")), data = dt, family = fam)
      # Build mcp() call dynamically: mcp(group = "Tukey")
      mcps <- list("Tukey"); names(mcps) <- groups[1]
      linf <- do.call(multcomp::mcp, mcps)
      res <- tryCatch(multcomp::glht(m2, linfct = linf), error = function(e) NULL)
      return(if (is.null(res)) NULL else summary(res))
    }
    
    if (isTRUE(used_interaction)) {
      lv <- levels(dt[["interaction_term"]]); if (length(lv) < 3) return(NULL)
      mcps <- list("Tukey"); names(mcps) <- "interaction_term"
      linf <- do.call(multcomp::mcp, mcps)
      res <- tryCatch(multcomp::glht(model, linfct = linf), error = function(e) NULL)
      return(if (is.null(res)) NULL else summary(res))
    }
    
    # Full factorial: not running omnibus Tukey over all cells here (prefer emmeans)
    NULL
  }
  
  .posthoc_emmeans <- function(model, groups, used_interaction, full_factorial, emmeans_simple_for = NULL) {
    if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)
    
    if (length(groups) == 1) {
      em <- emmeans::emmeans(model, specs = groups[1])
      return(list(
        table = summary(emmeans::contrast(em, method = "tukey")),
        emm_obj = em
      ))
    }
    
    if (isTRUE(full_factorial)) {
      # FIX: specs must be a one-sided formula (~ A | B), not "A | B"
      if (is.null(emmeans_simple_for)) emmeans_simple_for <- groups[1]
      # ensure valid target and conditioning factors
      if (!emmeans_simple_for %in% groups) emmeans_simple_for <- groups[1]
      other <- setdiff(groups, emmeans_simple_for)
      if (length(other) == 0) other <- groups[1]
      
      spec <- as.formula(paste("~", emmeans_simple_for, "|", other[1]))
      em <- emmeans::emmeans(model, specs = spec)
      return(list(
        table = summary(emmeans::contrast(em, method = "pairwise", adjust = "tukey")),
        emm_obj = em
      ))
    }
    
    if (isTRUE(used_interaction)) {
      em <- emmeans::emmeans(model, specs = "interaction_term")
      return(list(
        table = summary(emmeans::contrast(em, method = "tukey")),
        emm_obj = em
      ))
    }
    
    NULL
  }
  # ======================= End Helpers =======================
  
  posthoc_method <- match.arg(posthoc_method)
  
  dt <- data.table::as.data.table(data)
  .check_inputs(dt, response_var, group_vars_vec)
  dt <- .coerce_groups_to_factor(dt, group_vars_vec)
  
  built <- .build_formula_and_data(dt, response_var, group_vars_vec, full_factorial)
  model <- .fit_glm(built$formula, built$data, family)
  
  model_resid <- residuals(model, type = "deviance")
  
  anova_obj <- .anova_car(model, sum_squares_type, family)
  anova_df  <- as.data.frame(anova_obj)
  
  eff <- .effect_sizes(model, family, response_var, built$data)
  
  model_stats <- list(
    AIC = AIC(model),
    BIC = BIC(model),
    NullDeviance = model$null.deviance,
    ResidualDeviance = model$deviance,
    DF_null = model$df.null,
    DF_residual = model$df.residual
  )
  
  posthoc <- NULL
  if (posthoc_method == "multcomp") {
    posthoc <- .posthoc_multcomp(model, built$data, group_vars_vec, built$used_interaction, response_var, family)
  } else if (posthoc_method == "emmeans") {
    posthoc <- .posthoc_emmeans(model, group_vars_vec, built$used_interaction, full_factorial, emmeans_simple_for)
  }
  
  conf_intervals <- .confint_safe(model)
  
  residuals_plot <- .residuals_plot(model, plot_residuals)
  plot_box       <- .boxplot(built$data, response_var, group_vars_vec, built$used_interaction, full_factorial, family)
  
  notes <- c()
  od_note <- .overdispersion_note(model, family)
  if (!is.null(od_note)) notes <- c(notes, od_note)
  if (family$family != "gaussian") {
    notes <- c(notes, "For non-Gaussian families, interpret boxplots cautiously; consider marginal means plots (e.g., emmeans) on response or link scale.")
  }
  
  list(
    model = model,
    model_stats = model_stats,
    anova_test = anova_df,
    posthoc = posthoc,
    effect_sizes = eff,
    confidence_intervals = conf_intervals,
    residuals = model_resid,
    residuals_plot = residuals_plot,
    boxplot = plot_box,
    notes = if (length(notes)) notes else NULL
  )
}
