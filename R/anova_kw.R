#' Kruskal–Wallis Test with Dunn Post-hoc and Visualization
#'
#' @description
#' Conduct a Kruskal–Wallis test to compare a numeric response across groups
#' formed by one or more categorical variables. Optionally compute a
#' normality diagnostic (Anderson–Darling) on ANOVA residuals **for
#' diagnostics only** (Kruskal–Wallis does not assume normality), run
#' Dunn's post-hoc test with p-value adjustment, and return publication-ready
#' ggplot2 boxplots (with optional faceting for multi-factor designs).
#'
#' @param data A \code{data.frame} or \code{data.table} containing the input data.
#' @param response_var Character scalar. The name of the numeric response variable.
#' @param group_vars_vec Character vector. One or more grouping (categorical) variables.
#'   If multiple are supplied, an internal interaction term is created for inference
#'   and the boxplot is optionally faceted by the remaining group variables.
#' @param plot_qq Logical; if \code{TRUE}, returns a ggplot2 QQ-plot of ANOVA residuals
#'   (diagnostic only). Default \code{FALSE}.
#' @param posthoc_method Character; adjustment method for Dunn's post-hoc p-values.
#'   Passed to \code{dunn.test::dunn.test(method = ...)} and used to normalize the
#'   method for \code{stats::p.adjust()} fallback when needed. Accepted values include
#'   \code{"bh"}, \code{"fdr"}, \code{"bonferroni"}, \code{"holm"}, \code{"hochberg"},
#'   \code{"hommel"}, \code{"by"}, or \code{"none"}. Default \code{"bh"}.
#'
#' @returns A named \code{list} with:
#' \describe{
#'   \item{data_used}{\code{data.table} used in the analysis (copy of \code{data}, filtered
#'     to complete cases, with an \code{interaction_term} column and optional \code{facet_var}).}
#'   \item{assumptions}{\code{list} with \code{residuals} (from ANOVA fit on \code{interaction_term}),
#'     \code{ad_test} (Anderson–Darling on residuals), and \code{qqplot} (ggplot2 object or \code{NULL}).}
#'   \item{kruskal_test}{\code{htest} object from \code{stats::kruskal.test()}.}
#'   \item{posthoc}{\code{list} with \code{raw} (returned from \code{dunn.test}) and a
#'     \code{tidy} \code{data.frame} containing \code{comparison}, \code{Z},
#'     \code{p_unadjusted}, \code{p_adjusted}, \code{group1}, \code{group2}. May be \code{NULL}
#'     if fewer than two groups.}
#'   \item{plots}{\code{list} of ggplot2 objects: \code{boxplot} (always), \code{qqplot} (optional).}
#'   \item{n_dropped}{Integer count of rows removed due to missing values in response or group vars.}
#'   \item{notes}{Character vector of analysis notes (e.g., dropped rows, diagnostic caveats).}
#' }
#'
#' @details
#' \itemize{
#'   \item All rows with \code{NA} in \code{response_var} or any of \code{group_vars_vec} are removed
#'         prior to analysis. The number dropped is returned as \code{n_dropped}, and a note is included.
#'   \item Grouping variables are coerced to \code{factor} (if not already).
#'   \item An internal \code{interaction_term} combines \code{group_vars_vec} levels (with \code{sep=":"})
#'         and is used for KW and Dunn tests. If multiple grouping variables are provided,
#'         a \code{facet_var} (interaction of the remaining group variables) is created for faceting the boxplot.
#'   \item The ANOVA residual normality check (AD test) is for **diagnostic context** only;
#'         Kruskal–Wallis does not require normality.
#'   \item Dunn's post-hoc results are returned whether or not KW is significant (to match common workflows).
#'         You can inspect KW \code{p.value} to decide if you want to consider post-hoc inferences.
#'   \item The boxplot labels show group counts and may be reordered by median for readability.
#'         Faceting shows additional group combinations when \code{length(group_vars_vec) > 1}.
#'   \item The function does not print plots; it returns ggplot objects for the caller to \code{print()} if desired.
#' }
#'
#' @section Robustness Notes:
#' \itemize{
#'   \item The \code{dunn.test} package has varied element names across versions. This function
#'         extracts post-hoc results robustly and, if needed, computes adjusted p-values via
#'         \code{stats::p.adjust()} after normalizing the adjustment method string.
#'   \item No in-place mutation of the caller's data occurs; an internal copy is made.
#'   \item Where appropriate, \code{data.table::get()} is used within \code{data.table} verbs
#'         (rather than tidy eval pronouns) for reliability.
#' }
#'
#' @examples
#' \dontrun{
#' # Single-factor example
#' set.seed(123)
#' group <- rep(c("Group A", "Group B", "Group C"), each = 50)
#' value <- c(rnorm(50, mean = 5, sd = 1),
#'            rnorm(50, mean = 7, sd = 1.5),
#'            rnorm(50, mean = 4, sd = 0.8))
#' data_single <- data.frame(group = factor(group), value = value)
#' res1 <- anova_kw(data_single, "value", "group", plot_qq = TRUE, posthoc_method = "bh")
#' print(res1$plots$boxplot)
#' if (!is.null(res1$plots$qqplot)) print(res1$plots$qqplot)
#'
#' # Two-factor example (faceted boxplot)
#' group1 <- rep(c("Group A", "Group B"), each = 50)
#' group2 <- rep(c("Group X", "Group Y"), times = 50)
#' value2 <- c(rnorm(50, mean = 5, sd = 1),
#'             rnorm(50, mean = 7, sd = 1.5),
#'             rnorm(50, mean = 4, sd = 0.8),
#'             rnorm(50, mean = 6, sd = 1.2))
#' data_double <- data.frame(group1 = factor(group1),
#'                           group2 = factor(group2),
#'                           value = value2)
#' res2 <- anova_kw(data_double, "value", c("group1","group2"), plot_qq = FALSE, posthoc_method = "bonferroni")
#' print(res2$plots$boxplot)
#' }
#'
#' @seealso \code{\link[stats]{kruskal.test}}, \code{\link[dunn.test]{dunn.test}},
#'   \code{\link[nortest]{ad.test}}, \code{\link[ggplot2]{ggplot}}
#'
#' @export
#' @import data.table
#' @import ggplot2
#' @import nortest
#' @import dunn.test
anova_kw <- function(data, response_var, group_vars_vec, plot_qq = FALSE, posthoc_method = "bh") {
  
  #### Helper: normalize method for p.adjust() (case-insensitive) ####
  ._normalize_p_adjust_method <- function(m) {
    if (is.null(m) || length(m) != 1L) return("BH")
    m <- tolower(m)
    switch(m,
           "holm"       = "holm",
           "hochberg"   = "hochberg",
           "hommel"     = "hommel",
           "bonferroni" = "bonferroni",
           "bh"         = "BH",
           "fdr"        = "BH",   # alias
           "by"         = "BY",
           "none"       = "none",
           "BH"
    )
  }
  
  #### Helper: validate inputs ####
  .validate_inputs <- function(data, response_var, group_vars_vec) {
    if (!is.character(response_var) || length(response_var) != 1L)
      stop("response_var must be a single character string.")
    if (!(response_var %in% names(data)))
      stop("response_var not found in data.")
    if (!is.numeric(data[[response_var]]))
      stop("response_var must be numeric.")
    if (!all(group_vars_vec %in% names(data))) {
      missing_groups <- setdiff(group_vars_vec, names(data))
      stop("The following group_vars_vec elements are not in data: ",
           paste(missing_groups, collapse = ", "))
    }
    invisible(TRUE)
  }
  
  #### Helper: coerce & clean data (no side effects to caller) ####
  .prepare_data <- function(data, response_var, group_vars_vec, sep = ":") {
    DT <- data.table::as.data.table(data.table::copy(data))
    
    # Ensure grouping variables are factors
    for (g in group_vars_vec) {
      if (!is.factor(DT[[g]])) DT[[g]] <- factor(DT[[g]])
    }
    
    # NA filtering across response + groups
    keep <- stats::complete.cases(DT[, c(response_var, group_vars_vec), with = FALSE])
    n_dropped <- sum(!keep)
    DT <- DT[keep]
    
    # Create interaction term with readable labels
    DT[, interaction_term := do.call(
      interaction,
      c(lapply(.SD, as.character), list(drop = TRUE, sep = sep))
    ), .SDcols = group_vars_vec]
    DT[, interaction_term := factor(interaction_term, levels = unique(interaction_term))]
    
    # Optional facet var for additional grouping vars beyond the first
    if (length(group_vars_vec) > 1) {
      DT[, facet_var := do.call(
        interaction,
        c(lapply(.SD, as.character), list(drop = TRUE, sep = sep))
      ), .SDcols = group_vars_vec[-1]]
    }
    list(DT = DT, n_dropped = n_dropped)
  }
  
  #### Helper: diagnostics (ANOVA residuals + AD normality) ####
  .compute_residuals_and_ad <- function(DT, response_var) {
    model_formula <- stats::as.formula(paste(response_var, "~ interaction_term"))
    aov_fit <- stats::aov(model_formula, data = DT)
    resids <- stats::residuals(aov_fit)
    resids <- stats::na.omit(resids)
    ad <- nortest::ad.test(resids)
    list(residuals = resids, ad_test = ad)
  }
  
  #### Helper: QQ plot ####
  .qq_plot <- function(residuals) {
    df <- data.frame(residuals = residuals)
    ggplot2::ggplot(df, ggplot2::aes(sample = residuals)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line(linewidth = 0.8) +
      ggplot2::labs(
        title = "QQ Plot of ANOVA Residuals (diagnostic only)",
        x = "Theoretical Quantiles",
        y = "Sample Quantiles"
      ) +
      ggplot2::theme_minimal()
  }
  
  #### Helper: Kruskal–Wallis ####
  .kruskal <- function(DT, response_var) {
    ng <- DT[, data.table::uniqueN(interaction_term)]
    if (ng < 2) stop("Not enough groups for analysis. At least 2 unique groups with data are required.")
    stats::kruskal.test(stats::as.formula(paste(response_var, "~ interaction_term")), data = DT)
  }
  
  #### Helper: Dunn post-hoc + robust tidy ####
  .dunn_posthoc <- function(DT, response_var, method_user) {
    ng <- DT[, data.table::uniqueN(interaction_term)]
    if (ng < 2) {
      warning("Less than 2 groups found. Post-hoc analysis not performed.")
      return(NULL)
    }
    
    # Run dunn.test with the user's method (accepts lower-case like "bh")
    dres <- dunn.test::dunn.test(
      x = DT[[response_var]],
      g = DT[["interaction_term"]],
      method = method_user,
      altp = TRUE
    )
    
    # Robust extraction across dunn.test versions
    get_slot <- function(obj, primary, alternates = character()) {
      nm <- names(obj)
      if (primary %in% nm) return(obj[[primary]])
      for (alt in alternates) if (alt %in% nm) return(obj[[alt]])
      NULL
    }
    
    comps <- get_slot(dres, "comparisons", c("Comparison", "comparisons2"))
    Z     <- get_slot(dres, "Z",          c("z", "Zstat", "Z.statistics"))
    P     <- get_slot(dres, "P",          c("p", "P.unadj", "P.values"))
    Padj  <- get_slot(dres, "P.adjusted", c("P.adjusted2", "P.adj", "p.adjusted"))
    
    # Fallbacks
    if (is.null(comps) && !is.null(P)) {
      lv <- levels(DT$interaction_term)
      comb <- utils::combn(lv, 2, simplify = FALSE)
      comps <- vapply(comb, function(v) paste(v[1], "-", v[2]), character(1))
    }
    if (is.null(P) && !is.null(Z)) {
      P <- 2 * stats::pnorm(-abs(as.numeric(Z)))
    }
    if (!is.null(P) && is.null(Padj)) {
      method_padj <- ._normalize_p_adjust_method(method_user)
      Padj <- stats::p.adjust(as.numeric(P), method = method_padj)
    }
    
    # Coerce to vectors (may arrive as matrices)
    as_vec <- function(x) if (is.null(x)) numeric(0) else as.numeric(x)
    
    `%||%` <- function(x, y) if (is.null(x)) y else x
    
    comps <- as.character(comps %||% character(0))
    Z     <- as_vec(Z)
    P     <- as_vec(P)
    Padj  <- as_vec(Padj)
    
    # Length reconcile
    n <- length(comps)
    if (n == 0) {
      n <- max(length(Z), length(P), length(Padj))
      if (is.finite(n) && n > 0) comps <- rep_len("", n)
    }
    Z    <- if (length(Z)    == n) Z    else rep_len(NA_real_, n)
    P    <- if (length(P)    == n) P    else rep_len(NA_real_, n)
    Padj <- if (length(Padj) == n) Padj else rep_len(NA_real_, n)
    
    tidy <- data.frame(
      comparison   = comps,
      Z            = Z,
      p_unadjusted = P,
      p_adjusted   = Padj,
      stringsAsFactors = FALSE
    )
    
    if (n > 0 && any(nzchar(tidy$comparison))) {
      parts <- strsplit(tidy$comparison, " - ", fixed = TRUE)
      tidy$group1 <- vapply(parts, function(v) if (length(v) >= 1) v[1] else NA_character_, character(1))
      tidy$group2 <- vapply(parts, function(v) if (length(v) >= 2) v[2] else NA_character_, character(1))
    } else {
      tidy$group1 <- NA_character_
      tidy$group2 <- NA_character_
    }
    
    list(raw = dres, tidy = tidy)
  }
  
  #### Helper: Boxplot ####
  .boxplot <- function(DT, response_var, group_vars_vec, reorder_by_median = TRUE, show_counts = TRUE) {
    xvar <- group_vars_vec[1]
    
    if (show_counts) {
      counts <- DT[, .N, by = c(xvar)]
      DT <- merge(DT, counts, by = xvar, all.x = TRUE, sort = FALSE)
      DT[, x_label := sprintf("%s (n=%d)", get(xvar), N)]
      if (reorder_by_median) {
        med <- DT[, .(med = stats::median(get(response_var), na.rm = TRUE)), by = x_label]
        ord <- med[order(med), x_label]
        DT[, x_label := factor(x_label, levels = ord)]
      } else {
        DT[, x_label := factor(x_label, levels = unique(x_label))]
      }
      x_mapping <- "x_label"
    } else {
      if (reorder_by_median) {
        med <- DT[, .(med = stats::median(get(response_var), na.rm = TRUE)), by = xvar]
        ord <- med[order(med), get(xvar)]
        DT[, (xvar) := factor(get(xvar), levels = ord)]
      }
      x_mapping <- xvar
    }
    
    p <- ggplot2::ggplot(
      DT,
      ggplot2::aes(x = .data[[x_mapping]], y = .data[[response_var]], fill = .data[[xvar]])
    ) +
      ggplot2::geom_boxplot(outlier.alpha = 0.6) +
      ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(title = xvar)) +
      ggplot2::labs(
        title = "Distribution of Response by Group",
        subtitle = paste0("Response: ", response_var,
                          " | Grouping: ", paste(group_vars_vec, collapse = ", ")),
        x = xvar, y = response_var
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    if (length(group_vars_vec) > 1 && "facet_var" %in% names(DT)) {
      p <- p + ggplot2::facet_wrap(~ facet_var, scales = "free_x")
    }
    
    list(plot = p, data_for_plot = DT)
  }
  
  #### Main ####
  .validate_inputs(data, response_var, group_vars_vec)
  
  prep <- .prepare_data(data, response_var, group_vars_vec, sep = ":")
  DT   <- prep$DT
  n_dropped <- prep$n_dropped
  
  if (DT[, data.table::uniqueN(interaction_term)] < 2L) {
    stop("Not enough groups for analysis after filtering NAs. At least 2 groups are required.")
  }
  
  diag_res <- .compute_residuals_and_ad(DT, response_var)
  qqplot_obj <- if (plot_qq) .qq_plot(diag_res$residuals) else NULL
  
  kw <- .kruskal(DT, response_var)
  kw_sig <- !is.na(kw$p.value) && kw$p.value < 0.05
  
  dunn <- .dunn_posthoc(DT, response_var, posthoc_method)
  
  box_out <- .boxplot(DT, response_var, group_vars_vec, reorder_by_median = TRUE, show_counts = TRUE)
  
  notes <- c(
    sprintf("Dropped %d row(s) due to NA in response/group vars.", n_dropped),
    "Kruskal–Wallis is nonparametric; normality of residuals is a diagnostic only.",
    if (kw_sig) "KW is significant at alpha=0.05; examine Dunn's post-hoc."
    else "KW is not significant at alpha=0.05; post-hoc shown for completeness."
  )
  
  list(
    data_used     = DT,
    assumptions   = list(residuals = diag_res$residuals, ad_test = diag_res$ad_test, qqplot = qqplot_obj),
    kruskal_test  = kw,
    posthoc       = dunn,
    plots         = list(boxplot = box_out$plot, qqplot = qqplot_obj),
    n_dropped     = n_dropped,
    notes         = notes
  )
}
