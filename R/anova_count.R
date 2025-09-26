#' Analysis of Count Data with Poisson / Negative Binomial Models
#'
#' Fits a Poisson regression for count outcomes with one or more grouping
#' factors (including their interaction), assesses dispersion, optionally
#' refits using a Negative Binomial model, computes marginal means and
#' pairwise comparisons, and returns a plot of estimated means with
#' confidence intervals.
#'
#' @section Overview:
#' \itemize{
#'   \item Validates inputs and factor levels.
#'   \item Fits \code{glm(..., family = poisson(link = "log"))}.
#'   \item Computes the Pearson dispersion statistic.
#'   \item If dispersion exceeds a threshold and requested, refits with
#'         \code{MASS::glm.nb(...)}.
#'   \item Obtains estimated marginal means (EMMs) using \pkg{emmeans}
#'         on the response scale and pairwise comparisons with Tukey
#'         adjustment.
#'   \item Normalizes common column names from \pkg{emmeans} outputs so the
#'         returned tables consistently contain \code{response}, \code{lower.CL},
#'         and \code{upper.CL}.
#'   \item Provides a bar chart of estimated means (with intervals when available).
#' }
#'
#' @section Helpers:
#' Internal utilities used by the function include:
#' \itemize{
#'   \item Input checks and formula builder.
#'   \item Pearson dispersion computation.
#'   \item Optional HC0 variance-covariance matrix for inference.
#'   \item Name normalization for \pkg{emmeans} outputs.
#'   \item Assembly of EMM tables and pairwise contrasts (with fallbacks).
#'   \item A prediction-based fallback grid used if \pkg{emmeans} fails.
#'   \item Data preparation for plotting.
#' }
#'
#' @param data A \code{data.frame} containing the analysis variables.
#' @param response_var Character string. Name of the count response variable.
#'   Values must be non-negative integers; \code{NA}s are allowed and dropped.
#' @param group_vars Character vector. One or more grouping variables. When
#'   multiple variables are supplied, their interaction is modeled.
#' @param overdispersion_threshold Numeric. Pearson dispersion threshold used
#'   to flag overdispersion. Default \code{1.5}.
#' @param use_nb_if_overdispersed Logical. If \code{TRUE}, refit a Negative
#'   Binomial model via \code{MASS::glm.nb} when overdispersion is flagged.
#' @param vcov_type Character. Variance-covariance type: \code{"default"} (model
#'   based) or \code{"robust"} (HC0 via \pkg{sandwich}). Default \code{"default"}.
#' @param offset_var Optional character. Name of a strictly positive exposure
#'   variable to be used as a log-offset (e.g., time-at-risk).
#' @param ci_level Confidence level for intervals. Default \code{0.95}.
#' @param type_anova Character. \code{"LR"} or \code{"Wald"} for
#'   \code{car::Anova} test statistic. Default \code{"LR"}.
#' @param plot Logical. If \code{TRUE}, returns a \pkg{ggplot2} bar chart of
#'   estimated means by group combination. Default \code{TRUE}.
#' @param debug Logical. If \code{TRUE}, prints internal diagnostics. Default \code{FALSE}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{model}{The fitted object (\code{glm} or \code{negbin}).}
#'   \item{model_type}{Character, either \code{"poisson"} or \code{"negbin"}.}
#'   \item{overdispersion_statistic}{Pearson dispersion statistic.}
#'   \item{overdispersion_flagged}{Logical, whether dispersion exceeded threshold.}
#'   \item{anova_table}{ANOVA table from \code{car::Anova}, or \code{NULL}.}
#'   \item{emm_grid}{The \code{emmGrid} object (if constructed), or \code{NULL}.}
#'   \item{emmeans_table}{Data frame of estimated means with columns
#'         \code{response}, \code{lower.CL}, \code{upper.CL}.}
#'   \item{posthoc_pairs}{Data frame of Tukey pairwise comparisons with
#'         columns \code{contrast}, \code{IRR}, \code{lower.CL}, \code{upper.CL},
#'         \code{p.value}. May be empty if contrasts are not available.}
#'   \item{plot}{A \pkg{ggplot2} object when \code{plot = TRUE}, otherwise \code{NULL}.}
#'   \item{messages}{Character vector with analysis notes.}
#'   \item{n_removed_na}{Number of rows dropped due to missing values.}
#' }
#'
#' @details
#' The model is fit on complete cases across the response, grouping variables,
#' and the offset (if provided). For \pkg{emmeans}, results are requested on the
#' response scale; when column names differ across versions (e.g., \code{rate},
#' \code{asymp.LCL}, \code{asymp.UCL}), names are standardized to
#' \code{response}, \code{lower.CL}, and \code{upper.CL}. If \pkg{emmeans}
#' summaries or confidence intervals are unavailable, a prediction-based table
#' is produced from the fitted model using a full factorial grid.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 400
#' d <- data.frame(
#'   g1 = factor(rep(c("A","B"), each = n/2)),
#'   g2 = factor(rep(rep(c("X","Y"), each = n/4), 2))
#' )
#' lambda <- with(d, ifelse(g1=="A" & g2=="X", 5,
#'                   ifelse(g1=="A" & g2=="Y",10,
#'                   ifelse(g1=="B" & g2=="X",15,20))))
#' d$count <- rpois(n, lambda)
#'
#' out <- anova_count(
#'   data = d,
#'   response_var = "count",
#'   group_vars = c("g1","g2"),
#'   plot = TRUE
#' )
#' out$anova_table
#' out$emmeans_table
#' out$posthoc_pairs
#' out$plot
#' }
#'
#' @seealso \code{\link[stats]{glm}}, \code{\link[MASS]{glm.nb}},
#'   \code{\link[emmeans]{emmeans}}, \code{\link[car]{Anova}}
#'
#' @import stats
#' @import ggplot2
#' @importFrom MASS glm.nb
#' @importFrom car Anova
#' @importFrom emmeans emmeans
#' @importFrom sandwich vcovHC
#' @export
anova_count <- function(
    data,
    response_var,
    group_vars,
    overdispersion_threshold = 1.5,
    use_nb_if_overdispersed = TRUE,
    vcov_type = c("default", "robust"),
    offset_var = NULL,
    ci_level = 0.95,
    type_anova = c("LR", "Wald"),
    plot = TRUE,
    debug = FALSE
) {
  # ==== Helpers ================================================================
  .dbg <- function(...) if (isTRUE(debug)) cat("[anova_count]", sprintf(...), "\n")
  
  .is_whole_ignoring_na <- function(v) {
    if (!is.numeric(v)) return(FALSE)
    v_ok <- v[!is.na(v)]
    if (length(v_ok) == 0) return(TRUE)
    all(v_ok >= 0) && all(abs(v_ok - round(v_ok)) < .Machine$double.eps^0.5)
  }
  
  .build_formula <- function(resp, groups) {
    rhs <- paste(groups, collapse = " * ")
    stats::as.formula(paste(resp, "~", rhs))
  }
  
  .pearson_overdispersion <- function(fit) {
    pr <- stats::resid(fit, type = "pearson")
    sum(pr^2, na.rm = TRUE) / stats::df.residual(fit)
  }
  
  .vcov_matrix <- function(model, type = c("default","robust")) {
    type <- match.arg(type)
    if (type == "default") return(NULL)
    V <- NULL
    try({ V <- sandwich::vcovHC(model, type = "HC0") }, silent = TRUE)
    V
  }
  
  .normalize_ci_cols <- function(d) {
    if (!is.null(d) && is.data.frame(d)) {
      nm <- names(d)
      if ("asymp.LCL" %in% nm && !"lower.CL" %in% nm) names(d)[nm == "asymp.LCL"] <- "lower.CL"
      nm <- names(d)
      if ("asymp.UCL" %in% nm && !"upper.CL" %in% nm) names(d)[nm == "asymp.UCL"] <- "upper.CL"
      nm <- names(d)
      if ("LCL" %in% nm && !"lower.CL" %in% nm) names(d)[nm == "LCL"] <- "lower.CL"
      nm <- names(d)
      if ("UCL" %in% nm && !"upper.CL" %in% nm) names(d)[nm == "UCL"] <- "upper.CL"
    }
    d
  }
  
  .normalize_response_col <- function(d) {
    if (is.null(d) || !is.data.frame(d)) return(d)
    nm <- names(d)
    if ("response" %in% nm) return(d)
    if ("rate" %in% nm) { names(d)[nm == "rate"] <- "response"; return(d) }
    if ("emmean" %in% nm) { names(d)[nm == "emmean"] <- "response"; return(d) }
    d
  }
  
  .ensure_emmeans_ci <- function(obj, level = 0.95) {
    # summary on response scale
    summ <- try(suppressWarnings(summary(obj, infer = c(TRUE, TRUE), level = level, type = "response")),
                silent = TRUE)
    d_base <- NULL
    if (!inherits(summ, "try-error")) {
      d_base <- as.data.frame(summ)
      d_base <- .normalize_response_col(.normalize_ci_cols(d_base))
    } else {
      # fallback summaries
      d_base <- try(as.data.frame(summary(obj, type = "response")), silent = TRUE)
      if (inherits(d_base, "try-error")) d_base <- try(as.data.frame(summary(obj)), silent = TRUE)
      if (inherits(d_base, "try-error")) d_base <- NULL
      if (!is.null(d_base)) d_base <- .normalize_response_col(.normalize_ci_cols(d_base))
    }
    # add CIs via confint if missing
    if (!is.null(d_base) && !all(c("lower.CL","upper.CL") %in% names(d_base))) {
      ci <- try(as.data.frame(confint(obj, level = level, type = "response")), silent = TRUE)
      if (!inherits(ci, "try-error")) {
        ci <- .normalize_response_col(.normalize_ci_cols(ci))
        fac_cols <- intersect(names(d_base), names(ci))
        fac_cols <- setdiff(fac_cols, c("response","SE","df","t.ratio","z.ratio","p.value",
                                        "lower.CL","upper.CL","estimate","emmean"))
        if (length(fac_cols) == 0) {
          fac_cols <- setdiff(names(ci), c("response","SE","df","t.ratio","z.ratio","p.value",
                                           "lower.CL","upper.CL","estimate","emmean"))
        }
        if (length(fac_cols) > 0) {
          keep_ci <- c(fac_cols, intersect(c("lower.CL","upper.CL"), names(ci)))
          d_base <- merge(d_base, ci[, keep_ci, drop = FALSE], by = fac_cols, all.x = TRUE)
        }
      }
    }
    d_base
  }
  
  .pairs_to_df <- function(obj, level = 0.95) {
    prs <- try(emmeans::pairs(obj, adjust = "tukey"), silent = TRUE)
    if (inherits(prs, "try-error")) {
      prs <- try(emmeans::contrast(obj, "pairwise", adjust = "tukey"), silent = TRUE)
      if (inherits(prs, "try-error")) return(NULL)
    }
    # response-scale summary with CIs
    summ <- try(suppressWarnings(summary(prs, infer = c(TRUE, TRUE), level = level, type = "response")),
                silent = TRUE)
    if (!inherits(summ, "try-error")) {
      d <- as.data.frame(summ)
      if ("ratio" %in% names(d) && !"IRR" %in% names(d)) names(d)[names(d) == "ratio"] <- "IRR"
      d <- .normalize_ci_cols(d)
      if (!"IRR" %in% names(d) && "estimate" %in% names(d)) d$IRR <- exp(d$estimate)
      return(d)
    }
    # fallback: link-scale + response CIs
    d0 <- try(as.data.frame(summary(prs)), silent = TRUE)
    if (inherits(d0, "try-error")) return(NULL)
    if (!"IRR" %in% names(d0) && "estimate" %in% names(d0)) d0$IRR <- exp(d0$estimate)
    ci <- try(as.data.frame(confint(prs, level = level, type = "response")), silent = TRUE)
    if (!inherits(ci, "try-error")) {
      ci <- .normalize_ci_cols(ci)
      key <- if ("contrast" %in% intersect(names(d0), names(ci))) "contrast" else intersect(names(d0), names(ci))[1]
      if (length(key)) d0 <- merge(d0, ci[, c(key, intersect(c("lower.CL","upper.CL"), names(ci)))],
                                   by = key, all.x = TRUE)
    }
    .normalize_ci_cols(d0)
  }
  
  .emm_from_predict <- function(fit, groups, offset_var = NULL, level = 0.95) {
    levs <- lapply(groups, function(g) levels(model.frame(fit)[[g]]))
    names(levs) <- groups
    grid <- do.call(expand.grid, c(levs, stringsAsFactors = FALSE))
    newdata <- grid
    if (!is.null(offset_var)) newdata[[offset_var]] <- 1
    pr <- stats::predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
    z <- stats::qnorm(0.5 + level/2)
    response  <- exp(pr$fit)
    lower.CL  <- exp(pr$fit - z * pr$se.fit)
    upper.CL  <- exp(pr$fit + z * pr$se.fit)
    out <- data.frame(grid, response = as.numeric(response),
                      lower.CL = as.numeric(lower.CL),
                      upper.CL = as.numeric(upper.CL),
                      stringsAsFactors = FALSE)
    for (g in groups) out[[g]] <- base::droplevels(as.factor(out[[g]]))
    out
  }
  
  .coerce_plot_df <- function(d, factors, response_col = "response") {
    keep <- intersect(c(factors, response_col, "lower.CL", "upper.CL"), names(d))
    out <- as.data.frame(d[keep])
    for (f in factors) out[[f]] <- base::droplevels(as.factor(out[[f]]))
    out[[response_col]] <- as.numeric(out[[response_col]])
    if ("lower.CL" %in% names(out)) out[["lower.CL"]] <- as.numeric(out[["lower.CL"]])
    if ("upper.CL" %in% names(out)) out[["upper.CL"]] <- as.numeric(out[["upper.CL"]])
    out$interaction_label <- do.call(interaction, c(out[factors], list(drop = TRUE, sep = " : ")))
    out
  }
  
  # ==== Input validation =======================================================
  if (missing(data) || !is.data.frame(data) || nrow(data) == 0)
    stop("`data` must be a non-empty data.frame.")
  if (missing(response_var) || !is.character(response_var) || length(response_var) != 1)
    stop("`response_var` must be a single character string.")
  if (!response_var %in% names(data))
    stop(sprintf("Response variable '%s' not found.", response_var))
  if (missing(group_vars) || !is.character(group_vars) || length(group_vars) < 1)
    stop("`group_vars` must be a character vector of one or more grouping variables.")
  if (!all(group_vars %in% names(data))) {
    missing_vars <- group_vars[!group_vars %in% names(data)]
    stop(sprintf("Grouping variable(s) missing: %s", paste(missing_vars, collapse = ", ")))
  }
  
  resp_vec <- data[[response_var]]
  if (!is.numeric(resp_vec))
    stop(sprintf("`%s` must be numeric (non-negative integer counts).", response_var))
  if (!.is_whole_ignoring_na(resp_vec))
    stop(sprintf("`%s` must be non-negative integer counts (NAs allowed and will be dropped).", response_var))
  
  cols_needed <- c(response_var, group_vars)
  if (!is.null(offset_var)) cols_needed <- c(cols_needed, offset_var)
  data_local <- data[, cols_needed, drop = FALSE]
  
  for (g in group_vars) data_local[[g]] <- as.factor(data_local[[g]])
  
  complete_idx <- stats::complete.cases(data_local)
  n_removed_na <- sum(!complete_idx)
  data_local <- data_local[complete_idx, , drop = FALSE]
  if (nrow(data_local) == 0)
    stop("All rows were removed due to missing data in analyzed columns.")
  
  for (g in group_vars) {
    if (nlevels(base::droplevels(data_local[[g]])) < 2)
      stop(sprintf("Grouping variable '%s' must have at least 2 levels after removing NAs.", g))
  }
  if (!is.null(offset_var)) {
    if (!is.numeric(data_local[[offset_var]]) || any(data_local[[offset_var]] <= 0))
      stop("`offset_var` must be numeric and strictly positive (for log-offset) after removing NAs.")
  }
  
  # Sparse/zero cell note
  messages <- character(0)
  {
    levs <- lapply(group_vars, function(g) levels(data_local[[g]])); names(levs) <- group_vars
    full <- do.call(expand.grid, c(levs, stringsAsFactors = FALSE))
    obs <- stats::aggregate(rep(1, nrow(data_local)), data_local[group_vars], sum, drop = FALSE)
    names(obs)[ncol(obs)] <- "Freq"
    merged <- merge(full, obs, by = group_vars, all.x = TRUE)
    merged$Freq[is.na(merged$Freq)] <- 0L
    small <- sum(merged$Freq > 0L & merged$Freq < 5L)
    zero  <- sum(merged$Freq == 0L)
    messages <- c(messages,
                  sprintf("small counts check: %d cells < 5; zero observations check: %d cells.",
                          small, zero))
    if (zero > 0) messages <- c(messages, "Some interaction cells have zero observations; certain contrasts may be undefined or unstable.")
    if (small > 0) messages <- c(messages, "Some interaction cells have very small counts (< 5 rows); inference may be unstable.")
  }
  
  # ==== Model fitting ==========================================================
  type_anova <- match.arg(type_anova)
  vcov_type <- match.arg(vcov_type)
  
  fml <- .build_formula(response_var, group_vars)
  .dbg("Formula: %s", deparse(fml))
  
  poisson_fit <- stats::glm(
    formula = fml,
    data    = data_local,
    family  = stats::poisson(link = "log"),
    offset  = if (!is.null(offset_var)) log(data_local[[offset_var]]) else NULL
  )
  
  disp <- .pearson_overdispersion(poisson_fit)
  overdispersion_flag <- is.finite(disp) && (disp > overdispersion_threshold)
  .dbg("Dispersion: %.3f (threshold %.3f) -> flagged=%s",
       disp, overdispersion_threshold, as.character(overdispersion_flag))
  
  final_fit <- poisson_fit
  model_type <- "poisson"
  if (overdispersion_flag && isTRUE(use_nb_if_overdispersed)) {
    nb_fit <- try(MASS::glm.nb(
      formula = fml,
      data    = data_local,
      link    = log,
      offset  = if (!is.null(offset_var)) log(data_local[[offset_var]]) else NULL
    ), silent = TRUE)
    
    if (!inherits(nb_fit, "try-error")) {
      final_fit <- nb_fit
      model_type <- "negbin"
      messages <- c(messages, sprintf("Overdispersion detected (%.3f). Refit with Negative Binomial.", disp))
    } else {
      messages <- c(messages, sprintf("Overdispersion detected (%.3f), but NB refit failed; proceeding with Poisson.", disp))
    }
  } else if (overdispersion_flag) {
    messages <- c(messages, sprintf("Overdispersion detected (%.3f). Consider alternative variance estimators or a Negative Binomial model.", disp))
  }
  
  anova_table <- NULL
  try({ anova_table <- car::Anova(final_fit, test.statistic = type_anova) }, silent = TRUE)
  if (is.null(anova_table)) messages <- c(messages, "ANOVA failed; returning NULL.")
  
  vc <- .vcov_matrix(final_fit, vcov_type)
  if (!is.null(vc)) .dbg("Using HC0 variance-covariance matrix")
  
  emm_spec <- if (length(group_vars) == 1)
    stats::as.formula(paste("~", group_vars))
  else
    stats::as.formula(paste("~", paste(group_vars, collapse = ":")))
  .dbg("EMM spec: %s", deparse(emm_spec))
  
  em_grid <- NULL
  emmeans_table <- NULL
  posthoc_pairs <- NULL
  
  em_ok <- TRUE
  try({
    em_grid <- emmeans::emmeans(
      final_fit,
      specs = emm_spec,
      type  = "response",
      vcov  = vc
    )
  }, silent = TRUE)
  if (is.null(em_grid)) {
    em_ok <- FALSE
    messages <- c(messages, "emmeans failed; marginal means and post-hoc not available.")
  }
  
  if (em_ok) {
    emmeans_table <- .ensure_emmeans_ci(em_grid, level = ci_level)
    posthoc_pairs <- .pairs_to_df(em_grid, level = ci_level)
    if (isTRUE(debug) && !is.null(emmeans_table)) {
      .dbg("emmeans_table names: %s", paste(names(emmeans_table), collapse = ", "))
      print(utils::head(emmeans_table))
    }
    if (isTRUE(debug) && !is.null(posthoc_pairs)) {
      .dbg("posthoc_pairs names: %s", paste(names(posthoc_pairs), collapse = ", "))
      print(utils::head(posthoc_pairs))
    }
  }
  
  # Fallbacks when emmeans is unavailable
  if (is.null(emmeans_table)) {
    emmeans_table <- .emm_from_predict(final_fit, group_vars, offset_var = offset_var, level = ci_level)
  }
  if (is.null(posthoc_pairs)) {
    posthoc_pairs <- data.frame(contrast = character(), IRR = numeric(),
                                lower.CL = numeric(), upper.CL = numeric(),
                                p.value = numeric(), stringsAsFactors = FALSE)
  }
  
  # ==== Plot ==================================================================
  plt <- NULL
  if (plot) {
    if (is.null(emmeans_table) && !is.null(em_grid)) {
      emmeans_table <- .ensure_emmeans_ci(em_grid, level = ci_level)
    }
    if (is.null(emmeans_table)) {
      emmeans_table <- .emm_from_predict(final_fit, group_vars, offset_var = offset_var, level = ci_level)
    }
    if (!is.null(emmeans_table)) {
      emmeans_table <- .normalize_ci_cols(.normalize_response_col(emmeans_table))
      factor_cols <- group_vars[group_vars %in% names(emmeans_table)]
      if (length(factor_cols) > 0 && "response" %in% names(emmeans_table)) {
        plot_df <- .coerce_plot_df(emmeans_table, factor_cols, response_col = "response")
        has_ci <- all(c("lower.CL","upper.CL") %in% names(plot_df))
        plt <- ggplot2::ggplot(plot_df,
                               ggplot2::aes(x = interaction_label, y = .data[["response"]])) +
          ggplot2::geom_col() +
          { if (has_ci) ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data[["lower.CL"]], ymax = .data[["upper.CL"]]),
            width = 0.2
          ) } +
          ggplot2::labs(
            x = paste(factor_cols, collapse = " × "),
            y = paste0("Estimated ", response_var, " (rate)"),
            title = "Model-estimated means with confidence intervals"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
          )
      } else {
        messages <- c(messages, "Plot could not be constructed; missing expected emmeans columns.")
      }
    } else {
      messages <- c(messages, "Plot skipped: emmeans_table unavailable.")
    }
  }
  
  # ==== Return ================================================================
  list(
    model                     = final_fit,
    model_type                = model_type,
    overdispersion_statistic  = disp,
    overdispersion_flagged    = isTRUE(overdispersion_flag),
    anova_table               = anova_table,
    emm_grid                  = em_grid,
    posthoc_pairs             = .normalize_ci_cols(posthoc_pairs),
    emmeans_table             = .normalize_ci_cols(.normalize_response_col(emmeans_table)),
    plot                      = plt,
    messages                  = messages,
    n_removed_na              = n_removed_na
  )
}
