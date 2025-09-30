#' Perform Repeated Measures ANOVA with Diagnostics, Sphericity, and Effect Sizes
#'
#' Fits a repeated-measures ANOVA using \code{afex::aov_ez} and provides:
#' - Normality diagnostics (residuals + Shapiro–Wilk)
#' - Valid sphericity diagnostics (Mauchly's test) via a multivariate-wide model
#' - Effect sizes (partial eta^2)
#' - Optional LM diagnostic plots
#' - Optional estimated marginal means (EMMs) with a simple plot
#'
#' @param data data.frame in long format.
#' @param subject Character. Subject identifier column.
#' @param dv Character. Dependent variable column.
#' @param within Character vector. Within-subject factor(s).
#' @param between Character vector or NULL. Between-subject factor(s).
#' @param factorize Logical. If TRUE, coerces subject/within/between to factors and dv to numeric. Default TRUE.
#' @param emm_specs Optional. Passed to \code{emmeans::emmeans(specs=...)}.
#' @param diagnostics Logical. If TRUE, produces standard LM diagnostic plots when possible. Default TRUE.
#' @param verbose Logical. If TRUE, prints progress messages. Default TRUE.
#' @param ... Passed to \code{afex::aov_ez}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{anova_obj}{afex_aov object with processed data attached.}
#'   \item{normality}{List with residuals and Shapiro–Wilk test.}
#'   \item{mauchly}{Mauchly's test table (if applicable), otherwise NULL.}
#'   \item{effect_sizes}{Partial eta^2 table.}
#'   \item{diagnostic_plots}{Recorded base R plot of LM diagnostics (if produced).}
#'   \item{emm_plot}{ggplot2 object of EMMs (if requested), otherwise NULL.}
#' }
#'
#' @export
#' @import afex
#' @import emmeans
#' @import car
#' @import effectsize
#' @import ggplot2
#' @importFrom stats shapiro.test residuals lm complete.cases
anova_rm <- function(
    data, subject, dv, within, between = NULL,
    factorize = TRUE, emm_specs = NULL,
    diagnostics = TRUE, verbose = TRUE, ...
) {
  
  # -------------------- helpers --------------------
  .msg <- function(...) if (isTRUE(verbose)) cat(...)
  
  .check_and_prepare <- function(data, subject, dv, within, between, factorize) {
    req_cols <- c(subject, dv, within, if (!is.null(between)) between)
    missing_cols <- setdiff(req_cols, names(data))
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    d <- data
    
    if (factorize) {
      d[[subject]] <- as.factor(d[[subject]])
      dv_old <- d[[dv]]
      d[[dv]] <- suppressWarnings(as.numeric(d[[dv]]))
      na_new <- sum(is.na(d[[dv]]) & !is.na(dv_old))
      if (na_new > 0) {
        warning(na_new, " NA(s) introduced by coercing '", dv, "' to numeric.")
      }
      d[within] <- lapply(d[within], as.factor)
      if (!is.null(between)) d[between] <- lapply(d[between], as.factor)
    }
    
    d
  }
  
  # Pre-balance subjects to avoid afex missing-ID warning
  .prebalance_subjects <- function(d, subject, dv, within, between = NULL, verbose = TRUE) {
    SEP <- "__SEP__"
    lab_df <- d[within]
    lab_df[within] <- lapply(lab_df[within], function(x) droplevels(as.factor(x)))
    lab_df[[".cell."]] <- do.call(interaction, c(lab_df[within], list(drop = TRUE, sep = SEP)))
    
    d2 <- cbind(d[c(subject, dv, if (!is.null(between)) between)], lab_df[".cell."])
    names(d2)[names(d2) == ".cell."] <- "cell"
    
    wide <- stats::reshape(d2, idvar = subject, timevar = "cell", v.names = dv, direction = "wide")
    
    dv_prefix <- paste0(dv, ".")
    esc <- gsub("([][\\^$.|()?*+{}])", "\\\\\\1", dv_prefix)
    resp_cols <- grep(paste0("^", esc), names(wide), value = TRUE)
    
    keep_ids <- wide[[subject]][stats::complete.cases(wide[, resp_cols, drop = FALSE])]
    dropped  <- setdiff(unique(d[[subject]]), keep_ids)
    
    if (length(dropped) > 0L && isTRUE(verbose)) {
      cat(
        "Prebalancing: dropping ", length(dropped), " subject(s) with incomplete within cells: ",
        paste(head(as.character(dropped), 10), collapse = ", "),
        if (length(dropped) > 10) ", ..." else "", "\n",
        sep = ""
      )
    }
    
    d[d[[subject]] %in% keep_ids, , drop = FALSE]
  }
  
  .fit_afex <- function(d, subject, dv, within, between, ...) {
    dots <- list(...)
    if (is.null(dots$na.rm)) dots$na.rm <- TRUE
    do.call(
      afex::aov_ez,
      c(list(
        id      = subject,
        dv      = dv,
        data    = d,
        within  = within,
        between = between
      ), dots)
    )
  }
  
  .get_residuals <- function(anova_obj) {
    if (!is.null(anova_obj$lm)) {
      return(stats::residuals(anova_obj$lm))
    }
    if (!is.null(anova_obj$aov)) {
      return(stats::residuals(anova_obj$aov))
    }
    NULL
  }
  
  .normality <- function(resid_values) {
    if (is.null(resid_values)) {
      list(residuals = NULL, shapiro_test = NULL)
    } else {
      list(
        residuals    = resid_values,
        shapiro_test = stats::shapiro.test(resid_values)
      )
    }
  }
  
  .plot_lm_diagnostics <- function(model) {
    if (is.null(model)) return(NULL)
    if (!inherits(model, "lm") || inherits(model, "mlm")) return(NULL)
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mfrow = c(2, 2))
    graphics::plot(model)
    grDevices::recordPlot()
  }
  
  # Robust sphericity/Mauchly via MLM on wide data
  .mauchly_test <- function(d, subject, dv, within, between) {
    keep_cols <- c(subject, dv, within, if (!is.null(between)) between)
    d <- d[stats::complete.cases(d[, keep_cols, drop = FALSE]), , drop = FALSE]
    
    # Count within-cells (product of levels); need >=3
    within_levels <- lapply(within, function(w) levels(as.factor(d[[w]])))
    n_cells <- prod(vapply(within_levels, length, integer(1)))
    if (n_cells < 3) return(NULL)
    
    # Build cell label
    SEP <- "__SEP__"
    lab_df <- d[within]
    lab_df[within] <- lapply(lab_df[within], function(x) droplevels(as.factor(x)))
    lab_df[[".cell."]] <- do.call(interaction, c(lab_df[within], list(drop = TRUE, sep = SEP)))
    
    d2 <- cbind(d[c(subject, dv, if (!is.null(between)) between)], lab_df[".cell."])
    names(d2)[names(d2) == ".cell."] <- "cell"
    
    wide <- stats::reshape(
      d2, idvar = subject, timevar = "cell", v.names = dv, direction = "wide"
    )
    
    dv_prefix <- paste0(dv, ".")
    esc <- gsub("([][\\^$.|()?*+{}])", "\\\\\\1", dv_prefix)
    resp_cols <- grep(paste0("^", esc), names(wide), value = TRUE)
    if (length(resp_cols) < 3) return(NULL)
    
    keep_resp <- vapply(resp_cols, function(col) any(!is.na(wide[[col]])), logical(1))
    resp_cols <- resp_cols[keep_resp]
    if (length(resp_cols) < 3) return(NULL)
    
    Y <- as.matrix(wide[, resp_cols, drop = FALSE])
    
    # Between formula
    if (is.null(between) || length(between) == 0) {
      form <- stats::as.formula("Y ~ 1")
    } else {
      form <- stats::as.formula(paste("Y ~", paste(between, collapse = " * ")))
      for (b in between) wide[[b]] <- as.factor(wide[[b]])
    }
    
    fit_mlm <- stats::lm(form, data = wide)  # Y is found in parent frame
    
    # Build idata by reversing label construction
    cell_names <- sub(paste0("^", esc), "", resp_cols, perl = TRUE)
    split_levels <- strsplit(cell_names, split = SEP, fixed = TRUE)
    
    pieces_len <- lengths(split_levels)
    consistent <- (length(unique(pieces_len)) == 1L) && (unique(pieces_len) == length(within))
    
    if (!consistent) {
      idata <- data.frame(tmp = factor(cell_names))
      names(idata) <- within[1]
    } else {
      idata <- as.data.frame(do.call(rbind, split_levels), stringsAsFactors = TRUE)
      names(idata) <- within
      idata[] <- lapply(idata, factor)
    }
    
    idesign <- stats::as.formula(paste("~", paste(within, collapse = " * ")))
    anv <- car::Anova(fit_mlm, idata = idata, idesign = idesign, type = "III")
    
    # Mauchly's test per within effect
    tryCatch(car::Mauchly.test(anv), error = function(e) NULL)
  }
  
  .eta_sq_safe <- function(anova_obj) {
    if (!is.null(anova_obj$Anova)) {
      es <- try(effectsize::eta_squared(anova_obj$Anova, partial = TRUE), silent = TRUE)
      if (!inherits(es, "try-error")) return(es)
    }
    es <- try(effectsize::eta_squared(anova_obj, partial = TRUE), silent = TRUE)
    if (!inherits(es, "try-error")) return(es)
    warning("Effect size calculation failed; returning NULL.")
    NULL
  }
  
  .emm_plotter <- function(anova_obj, emm_specs) {
    if (is.null(emm_specs)) return(NULL)
    emm_obj <- emmeans::emmeans(anova_obj, specs = emm_specs)
    df <- as.data.frame(emm_obj)
    if (ncol(df) < 2L || !"emmean" %in% names(df)) return(NULL)
    xvar <- names(df)[1L]
    
    # Avoid aes_string(): create a temporary .x column for plotting
    df$.x <- df[[xvar]]
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .x, y = emmean, group = 1)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_line() +
      {
        if ("SE" %in% names(df)) {
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = emmean - SE, ymax = emmean + SE),
            width = 0.15
          )
        } else {
          ggplot2::geom_blank()
        }
      } +
      ggplot2::labs(x = xvar, y = "Estimated Marginal Mean") +
      ggplot2::theme_minimal()
    p
  }
  
  # -------------------- workflow --------------------
  .msg("Preparing data...\n")
  data_prep <- .check_and_prepare(data, subject, dv, within, between, factorize)
  
  # NEW: pre-balance subjects to avoid afex's missing-ID warning banner
  data_prep <- .prebalance_subjects(data_prep, subject, dv, within, between, verbose)
  
  .msg("Fitting repeated-measures ANOVA model (afex::aov_ez)...\n")
  anova_obj <- .fit_afex(data_prep, subject, dv, within, between, ...)
  
  # Attach processed data for convenience
  anova_obj$data <- data_prep
  if (isTRUE(verbose)) print(anova_obj)
  
  # Normality diagnostics
  .msg("Computing normality diagnostics...\n")
  resid_values <- .get_residuals(anova_obj)
  normality <- .normality(resid_values)
  
  # Optional LM diagnostics
  diagnostic_plots <- NULL
  if (isTRUE(diagnostics) && !is.null(anova_obj$lm)) {
    .msg("Producing LM diagnostic plots...\n")
    diagnostic_plots <- tryCatch(.plot_lm_diagnostics(anova_obj$lm),
                                 error = function(e) { warning(e$message); NULL })
  }
  
  # Sphericity / Mauchly
  .msg("Evaluating sphericity (Mauchly's test) when applicable...\n")
  mauchly <- tryCatch(
    .mauchly_test(data_prep, subject, dv, within, between),
    error = function(e) { warning("Sphericity test failed: ", e$message); NULL }
  )
  
  # Effect sizes
  .msg("Computing partial eta^2 effect sizes...\n")
  effect_sizes <- .eta_sq_safe(anova_obj)
  
  # EMM plot (optional)
  emm_plot <- NULL
  if (!is.null(emm_specs)) {
    .msg("Computing and plotting estimated marginal means...\n")
    emm_plot <- tryCatch(.emm_plotter(anova_obj, emm_specs),
                         error = function(e) { warning(e$message); NULL })
  }
  
  # Return
  list(
    anova_obj        = anova_obj,
    normality        = normality,
    mauchly          = mauchly,
    effect_sizes     = effect_sizes,
    diagnostic_plots = diagnostic_plots,
    emm_plot         = emm_plot
  )
}
