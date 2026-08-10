#' Compare Conditional Standard Error of Measurement (CSEM) Curves
#'
#' This function compares two CSEM curves (or two methods) using:
#' - Polynomial regression with interaction to assess global differences.
#' - Point‑wise effect sizes (CSEM ratio) with confidence intervals.
#' - Point‑wise hypothesis tests with multiple comparison correction (if sample sizes are provided).
#'
#' @param score Numeric vector of scores (common axis). Must be increasing.
#' @param csem.m1 Vector of CSEM for method/group 1.
#' @param csem.m2 Vector of CSEM for method/group 2.
#' @param n.m1 Vector of sample sizes for method/group 1 (optional, required for tests and CIs).
#' @param n.m2 Vector of sample sizes for method/group 2 (optional).
#' @param name.m1 Name of method/group 1 (optional; if NULL, inferred from object).
#' @param name.m2 Name of method/group 2 (optional).
#' @param poly.degree Degree of polynomial regression (default 3).
#' @param nboot Number of bootstrap replicates for global ratio CI (default 1000, only used if n provided).
#' @param seed Random seed for reproducibility (default 123).
#' @param adjust.method Method for multiple comparison correction (default "fdr"). See \code{\link{p.adjust.methods}} for options.
#' @param alpha Significance level for significant regions (default 0.05).
#' @param conf.level Confidence level for intervals (default 0.95).
#' @param plot Logical; if TRUE generates plots (default FALSE).
#' @param plot.type Type of plot: "csem" (CSEM ratio with CI and significant regions) or "var" (variance ratio with CI and significant regions). Default "csem".
#' @param plot.theme A ggplot2 theme object (optional). If NULL, \code{theme_classic()} is used.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{csem.compare}}{Data frame with score, CSEM values, difference, ratio, and confidence limits (CSEM scale).}
#'   \item{\code{var.compare}}{Data frame with variance components, difference, ratio, standard errors, z‑statistics, p‑values (raw and adjusted), significance flag, and confidence limits (variance scale).}
#'   \item{\code{global.effect}}{Data frame with global summaries: geometric mean ratios and mean differences for both CSEM and variance, with bootstrap CIs.}
#'   \item{\code{regression}}{List with ANOVA table, fit statistics (R², adjusted R²), and coefficient table (estimates, SE, t, p, CI).}
#'   \item{\code{significant.regions}}{Data frame with intervals where adjusted p‑value < \code{alpha} (if \code{n} provided).}
#'   \item{\code{alpha}}{The significance level used.}
#'   \item{\code{conf.level}}{The confidence level used.}
#' }
#'
#' @details
#' The statistical tests are based on the error variances (CSEM²) because they are additive
#' and their sampling variance has a closed form:
#' \deqn{\mathrm{Var}(\hat{\sigma}^2) \approx \frac{2(\hat{\sigma}^2)^2}{n-1}}{Var(s²) ≈ 2*(s²)²/(n-1)}
#' The standard error of the difference between two independent variances is:
#' \deqn{\mathrm{SE}(\hat{\sigma}_1^2 - \hat{\sigma}_2^2) = \sqrt{ \frac{2\hat{\sigma}_1^4}{n_1-1} + \frac{2\hat{\sigma}_2^4}{n_2-1} }}
#' Confidence intervals for the ratio of CSEM are obtained by taking the square root of the
#' variance‑ratio confidence limits (after applying the log transformation).
#'
#' @examples
#' \donttest{
#' # Simulated example without sample sizes, and non plot
#' score <- seq(0, 100, length = 50)
#' csem1 <- 1 + 0.02*score + rnorm(50, 0, 0.1)
#' csem2 <- 1.2 + 0.015*score + rnorm(50, 0, 0.1)
#' res <- curvesCSEMcomp(score, csem1, csem2, name.m1 = "Method A", name.m2 = "Method B")
#'
#' # With sample sizes (to get tests and CIs)
#' res2 <- curvesCSEMcomp(score, csem1, csem2, n.m1 = rep(30,50), n.m2 = rep(30,50),
#'                        plot = TRUE, plot.type = "csem", alpha = 0.01)
#'
#' # Global effect output
#' res2$global.effect
#'
#' # Output in score-point level
#' res2$csem.compare
#' res2$var.compare
#' }
#'
#' @importFrom stats lm anova p.adjust pnorm quantile poly confint
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_hline labs theme_classic
#' @importFrom patchwork wrap_plots
#' @export
curvesCSEMcomp <- function(score,
                              csem.m1,
                              csem.m2,
                              n.m1 = NULL,
                              n.m2 = NULL,
                              name.m1 = NULL,
                              name.m2 = NULL,
                              poly.degree = 3,
                              nboot = 1000,
                              seed = 123,
                              adjust.method = "fdr",
                              alpha = 0.05,
                              conf.level = 0.95,
                              plot = FALSE,
                              plot.type = c("csem", "var"),
                              plot.theme = NULL) {

  # --------------------------------------------------------------------------
  # 1. Capture names and validate inputs
  # --------------------------------------------------------------------------
  if (is.null(name.m1)) {
    name.m1 <- deparse(substitute(csem.m1))
    if (length(name.m1) > 1 || name.m1 %in% c("csem.m1", "NULL")) name.m1 <- "Method 1"
  }
  if (is.null(name.m2)) {
    name.m2 <- deparse(substitute(csem.m2))
    if (length(name.m2) > 1 || name.m2 %in% c("csem.m2", "NULL")) name.m2 <- "Method 2"
  }

  # Length checks
  if (length(score) != length(csem.m1) || length(score) != length(csem.m2))
    stop("'score', 'csem.m1', and 'csem.m2' must have the same length.")
  if (!is.numeric(score)) stop("'score' must be numeric.")
  if (anyNA(score) || anyNA(csem.m1) || anyNA(csem.m2))
    stop("NA values are not allowed in 'score', 'csem.m1', or 'csem.m2'.")
  if (any(!is.finite(score)) || any(!is.finite(csem.m1)) || any(!is.finite(csem.m2)))
    stop("Infinite values are not allowed.")
  if (any(csem.m1 < 0) || any(csem.m2 < 0))
    stop("CSEM values must be non‑negative.")
  if (is.unsorted(score)) {
    warning("'score' is not sorted. Data will be reordered.")
    ord <- order(score)
    score <- score[ord]
    csem.m1 <- csem.m1[ord]
    csem.m2 <- csem.m2[ord]
    if (!is.null(n.m1)) n.m1 <- n.m1[ord]
    if (!is.null(n.m2)) n.m2 <- n.m2[ord]
  }

  use_n <- !is.null(n.m1) && !is.null(n.m2)
  if (use_n) {
    if (length(n.m1) != length(score) || length(n.m2) != length(score))
      stop("'n.m1' and 'n.m2' must have the same length as 'score'.")
    if (!is.numeric(n.m1) || !is.numeric(n.m2))
      stop("'n.m1' and 'n.m2' must be numeric.")
    if (anyNA(n.m1) || anyNA(n.m2))
      stop("NA values are not allowed in 'n.m1' or 'n.m2'.")
    if (any(n.m1 < 2) || any(n.m2 < 2))
      stop("Sample sizes must be >= 2.")
    n.m1 <- as.integer(n.m1)
    n.m2 <- as.integer(n.m2)
  }

  if (!is.numeric(poly.degree) || poly.degree < 1 || poly.degree != round(poly.degree))
    stop("'poly.degree' must be a positive integer.")
  if (nboot < 1 || nboot != round(nboot))
    stop("'nboot' must be a positive integer.")
  if (!(adjust.method %in% p.adjust.methods))
    stop("'adjust.method' must be one of: ", paste(p.adjust.methods, collapse = ", "))
  if (alpha <= 0 || alpha >= 1)
    stop("'alpha' must be between 0 and 1.")
  if (conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be between 0 and 1.")

  plot.type <- match.arg(plot.type)
  z_crit <- qnorm(1 - (1 - conf.level) / 2)

  # --------------------------------------------------------------------------
  # 2. Compute error variances
  # --------------------------------------------------------------------------
  var1 <- csem.m1^2
  var2 <- csem.m2^2
  eps <- 1e-10

  # --------------------------------------------------------------------------
  # 3. Polynomial regression with interaction
  # --------------------------------------------------------------------------
  df_long <- data.frame(
    score = rep(score, 2),
    var = c(var1, var2),
    method = factor(rep(c(name.m1, name.m2), each = length(score)))
  )

  form <- as.formula(paste0("var ~ poly(score, ", poly.degree, ") * method"))
  lm_mod <- lm(form, data = df_long)
  anova_res <- anova(lm_mod)
  coef_sum <- summary(lm_mod)$coefficients
  conf_int <- confint(lm_mod, level = conf.level)
  colnames(conf_int) <- c("lwr.ci", "upr.ci")
  coef_df <- data.frame(
    estimate = coef_sum[, "Estimate"],
    s.e. = coef_sum[, "Std. Error"],
    t.value = coef_sum[, "t value"],
    p.value = coef_sum[, "Pr(>|t|)"],
    lwr.ci = conf_int[, "lwr.ci"],
    upr.ci = conf_int[, "upr.ci"]
  )

  regression <- list(
    anova = anova_res,
    fit = data.frame(
      r.squared = summary(lm_mod)$r.squared,
      adj.r.squared = summary(lm_mod)$adj.r.squared
    ),
    coefficients = coef_df
  )

  # --------------------------------------------------------------------------
  # 4. Point‑wise effect sizes (CSEM and variance)
  # --------------------------------------------------------------------------
  diff_csem <- csem.m1 - csem.m2
  ratio_csem <- csem.m1 / (csem.m2 + eps)
  diff_var <- var1 - var2
  ratio_var <- var1 / (var2 + eps)

  # --------------------------------------------------------------------------
  # 5. Statistical tests and CIs (if n provided)
  # --------------------------------------------------------------------------
  if (use_n) {
    se_var1 <- sqrt(2 * var1^2 / (n.m1 - 1))
    se_var2 <- sqrt(2 * var2^2 / (n.m2 - 1))
    se_diff_var <- sqrt(se_var1^2 + se_var2^2)

    z_var <- diff_var / (se_diff_var + eps)
    p_var <- 2 * pnorm(-abs(z_var))
    p_adj <- p.adjust(p_var, method = adjust.method)
    sig <- p_adj < alpha

    # CI for variance ratio (log scale)
    se_log_ratio_var <- sqrt(2 / (n.m1 - 1) + 2 / (n.m2 - 1))
    log_ratio_var <- log(ratio_var + eps)
    ci_lower_log_var <- log_ratio_var - z_crit * se_log_ratio_var
    ci_upper_log_var <- log_ratio_var + z_crit * se_log_ratio_var
    ci_lower_var <- exp(ci_lower_log_var)
    ci_upper_var <- exp(ci_upper_log_var)

    # CI for CSEM ratio (square root of variance ratio limits)
    ci_lower_csem <- sqrt(ci_lower_var)
    ci_upper_csem <- sqrt(ci_upper_var)

    # CI for variance difference
    ci_lower_diff_var <- diff_var - z_crit * se_diff_var
    ci_upper_diff_var <- diff_var + z_crit * se_diff_var

    # Build csem.compare
    csem_compare <- data.frame(
      score = score,
      csem.m1 = csem.m1,
      csem.m2 = csem.m2,
      diff.csem = diff_csem,
      ratio.csem = ratio_csem,
      lwr.ci.csem = ci_lower_csem,
      upr.ci.csem = ci_upper_csem
    )

    # Build var.compare (including SEs for plotting)
    var_compare <- data.frame(
      score = score,
      var.m1 = var1,
      var.m2 = var2,
      diff.var = diff_var,
      ratio.var = ratio_var,
      se.var1 = se_var1,
      se.var2 = se_var2,
      se.diff.var = se_diff_var,
      z.var = z_var,
      p.value = p_var,
      p.adj = p_adj,
      significant = sig,
      lwr.ci.var = ci_lower_var,
      upr.ci.var = ci_upper_var
    )

    # Significant regions
    sig_idx <- which(sig)
    if (length(sig_idx) > 0) {
      runs <- split(sig_idx, cumsum(c(1, diff(sig_idx) != 1)))
      significant_regions <- do.call(rbind, lapply(runs, function(r) {
        data.frame(start = score[r[1]], end = score[r[length(r)]])
      }))
    } else {
      significant_regions <- NULL
    }

  } else {
    # No n: only descriptive tables (no tests)
    csem_compare <- data.frame(
      score = score,
      csem.m1 = csem.m1,
      csem.m2 = csem.m2,
      diff.csem = diff_csem,
      ratio.csem = ratio_csem
    )

    var_compare <- data.frame(
      score = score,
      var.m1 = var1,
      var.m2 = var2,
      diff.var = diff_var,
      ratio.var = ratio_var
    )

    significant_regions <- NULL
  }

  # --------------------------------------------------------------------------
  # 6. Global effect summaries
  # --------------------------------------------------------------------------
  # Geometric mean of CSEM ratio
  valid_csem <- ratio_csem > 0 & is.finite(ratio_csem)
  if (any(valid_csem)) {
    log_ratio_csem <- log(ratio_csem[valid_csem])
    geom_csem <- exp(mean(log_ratio_csem))
    if (use_n) {
      set.seed(seed)
      boot_csem <- replicate(nboot, mean(sample(log_ratio_csem, replace = TRUE)))
      ci_csem <- quantile(boot_csem, probs = c((1 - conf.level)/2, 1 - (1 - conf.level)/2), na.rm = TRUE)
      ci_lower_csem_g <- exp(ci_csem[1])
      ci_upper_csem_g <- exp(ci_csem[2])
    } else {
      ci_lower_csem_g <- NA
      ci_upper_csem_g <- NA
    }
  } else {
    geom_csem <- NA
    ci_lower_csem_g <- NA
    ci_upper_csem_g <- NA
  }

  # Geometric mean of variance ratio
  valid_var <- ratio_var > 0 & is.finite(ratio_var)
  if (any(valid_var)) {
    log_ratio_var_all <- log(ratio_var[valid_var])
    geom_var <- exp(mean(log_ratio_var_all))
    if (use_n) {
      set.seed(seed)
      boot_var <- replicate(nboot, mean(sample(log_ratio_var_all, replace = TRUE)))
      ci_var <- quantile(boot_var, probs = c((1 - conf.level)/2, 1 - (1 - conf.level)/2), na.rm = TRUE)
      ci_lower_var_g <- exp(ci_var[1])
      ci_upper_var_g <- exp(ci_var[2])
    } else {
      ci_lower_var_g <- NA
      ci_upper_var_g <- NA
    }
  } else {
    geom_var <- NA
    ci_lower_var_g <- NA
    ci_upper_var_g <- NA
  }

  # Mean differences (CSEM and variance)
  mean_diff_csem <- mean(diff_csem, na.rm = TRUE)
  mean_diff_var <- mean(diff_var, na.rm = TRUE)

  if (use_n) {
    set.seed(seed)
    boot_diff_csem <- replicate(nboot, mean(sample(diff_csem, replace = TRUE)))
    ci_diff_csem <- quantile(boot_diff_csem, probs = c((1 - conf.level)/2, 1 - (1 - conf.level)/2), na.rm = TRUE)
    boot_diff_var <- replicate(nboot, mean(sample(diff_var, replace = TRUE)))
    ci_diff_var <- quantile(boot_diff_var, probs = c((1 - conf.level)/2, 1 - (1 - conf.level)/2), na.rm = TRUE)
  } else {
    ci_diff_csem <- c(NA, NA)
    ci_diff_var <- c(NA, NA)
  }

  global_effect <- data.frame(
    metric = c("CSEM", "Variance"),
    geom.mean.ratio = c(geom_csem, geom_var),
    lwr.ci.ratio = c(ci_lower_csem_g, ci_lower_var_g),
    upr.ci.ratio = c(ci_upper_csem_g, ci_upper_var_g),
    mean.diff = c(mean_diff_csem, mean_diff_var),
    lwr.ci.diff = c(ci_diff_csem[1], ci_diff_var[1]),
    upr.ci.diff = c(ci_diff_csem[2], ci_diff_var[2])
  )

  # --------------------------------------------------------------------------
  # 7. Plots (if requested)
  # --------------------------------------------------------------------------
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE))
      stop("Package 'ggplot2' is required for plotting.")

    # Helper function to generate ratio plot (works for both CSEM and variance)
    plot_ratio <- function(df, ratio_col, lwr_col, upr_col,
                           x_lab = "Score", y_lab = NULL,
                           title = NULL) {
      y_ratio <- df[[ratio_col]]
      y_lwr <- df[[lwr_col]]
      y_upr <- df[[upr_col]]

      p <- ggplot2::ggplot(df, ggplot2::aes(x = score)) +
        ggplot2::geom_line(ggplot2::aes(y = y_ratio), color = "purple", size = 1.2) +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
        ggplot2::geom_hline(yintercept = c(0.8, 1.2), linetype = "dotted", color = "gray50", alpha = 0.6) +
        ggplot2::labs(
          title = title,
          x = x_lab,
          y = y_lab
        ) +
        (if (is.null(plot.theme)) ggplot2::theme_classic() else plot.theme)

      if (use_n) {
        p <- p +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = y_lwr, ymax = y_upr),
                               alpha = 0.2, fill = "steelblue")
        if (!is.null(significant_regions)) {
          for (i in 1:nrow(significant_regions)) {
            p <- p +
              ggplot2::annotate("rect",
                                xmin = significant_regions$start[i],
                                xmax = significant_regions$end[i],
                                ymin = -Inf, ymax = Inf,
                                alpha = 0.15, fill = "steelblue")
          }
          # Add label for significant regions
          y_max <- max(y_ratio, na.rm = TRUE) * 0.9
          p <- p +
            ggplot2::annotate("text", x = min(df$score), y = y_max,
                              label = paste0("Significant regions (FDR < ", alpha, ")"),
                              hjust = 0, color = "steelblue", size = 3)
        }
      }
      p
    }

    if (plot.type == "csem") {
      # Plot CSEM ratio
      p <- plot_ratio(
        df = csem_compare,
        ratio_col = "ratio.csem",
        lwr_col = "lwr.ci.csem",
        upr_col = "upr.ci.csem",
        y_lab = paste0("CSEM ratio (", name.m1, " / ", name.m2, ")"),
        title = "Effect size: CSEM ratio"
      )
      print(p)
    } else { # plot.type == "var"
      # Plot variance ratio
      p <- plot_ratio(
        df = var_compare,
        ratio_col = "ratio.var",
        lwr_col = "lwr.ci.var",
        upr_col = "upr.ci.var",
        y_lab = paste0("Variance ratio (", name.m1, " / ", name.m2, ")"),
        title = "Effect size: Variance ratio"
      )
      print(p)
    }
  }

  # --------------------------------------------------------------------------
  # 8. Assemble output
  # --------------------------------------------------------------------------
  list(
    csem.compare = csem_compare,
    var.compare = var_compare,
    global.effect = global_effect,
    regression = regression,
    significant.regions = significant_regions,
    alpha = alpha,
    conf.level = conf.level,
    poly.degree = poly.degree
  ) -> result

  class(result) <- "csem_compare"
  return(result)
}


#' Print method for csem_compare objects
#'
#' @param x An object of class \code{csem_compare}.
#' @param details Logical; if \code{TRUE}, show regression coefficients.
#' @param ... Additional arguments passed to \code{print}.
#' @export
print.csem_compare <- function(x, details = FALSE, ...) {
  cat("Comparison of CSEM curves\n")
  cat("\n--- Global effects ---\n")
  print(x$global.effect, digits = 4, row.names = FALSE)

  cat("\n--- Polynomial regression ---\n")
  cat("Polynomial degree:", x$poly.degree, "\n")
  cat("R²:", round(x$regression$fit$r.squared, 4), "\n")
  cat("p‑value for method effect:",
      round(x$regression$anova["method", "Pr(>F)"], 4), "\n")
  inter_row <- grep("poly.*:method", rownames(x$regression$anova))
  if (length(inter_row) > 0) {
    cat("p‑value for interaction (different shape):",
        round(x$regression$anova[inter_row, "Pr(>F)"], 4), "\n")
  }

  if (!is.null(x$significant.regions)) {
    cat("\n--- Significant regions (FDR <", x$alpha, ") ---\n", sep = "")
    print(x$significant.regions, row.names = FALSE)
  } else {
    cat("\nNo significant regions found (or sample sizes not provided).\n")
  }

  if (details) {
    cat("\n--- Regression coefficients (", x$conf.level * 100, "% CI) ---\n", sep = "")
    print(x$regression$coefficients, digits = 4)
  } else {
    cat("\nTo see regression coefficients, use: print(x, details = TRUE)\n")
  }

  cat("\nFor full point‑wise tables, use: x$csem.compare  or  x$var.compare\n")
  invisible(x)
}
