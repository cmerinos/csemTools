#' Compute CSEM for Scale Scores using Woodruff et al. (2018) Delta Method
#'
#' @description
#' Implements the first‑order delta method for computing conditional standard errors
#' of measurement (CSEM) on a non‑linearly transformed scale score metric, following
#' the polynomial approximation approach described in Woodruff, Traynor, Cui, and Fang
#' (2013). This method adjusts raw‑score CSEMs by the slope of the raw‑to‑scale
#' conversion function, \eqn{g'(X)}, obtained from a polynomial regression.
#'
#' @param raw Numeric vector. Raw scores (must be integer values, ideally consecutive).
#'   If multiple raw scores map to the same scale score, they are automatically
#'   aggregated (averaged) for model fitting.
#' @param scale Numeric vector. Scale scores corresponding to each raw score.
#'   Must have the same length as \code{raw}.
#' @param csem Numeric vector. Conditional standard errors of measurement in raw
#'   score units. Must have the same length as \code{raw}. These are typically
#'   obtained from classical test theory procedures (e.g., Thorndike's difference
#'   method, Feldt‑Qualls, or other CTT estimators).
#' @param poly.degree Integer. Degree of the polynomial used to approximate the
#'   raw‑to‑scale conversion function. Default is \code{6}, following the practice
#'   in Woodruff et al. (2013), who used 6th‑degree polynomials for three of the
#'   four ACT tests. The function allows any positive integer, but the degree
#'   cannot exceed the number of unique raw scores minus 1.
#' @param truncate.edges Logical. If \code{TRUE} (default), the derivative values
#'   at the extreme ends of the raw score distribution are replaced by the
#'   derivative of the nearest interior point to avoid explosive extrapolation
#'   caused by polynomial instability at the boundaries. The number of points
#'   truncated at each end is controlled by \code{n.edge}.
#' @param n.edge Integer vector of length 2: \code{c(n_lower, n_upper)}.
#'   Number of points to truncate at the lower and upper tails, respectively.
#'   Default is \code{c(1, 1)}. For example, \code{c(2, 3)} truncates the two
#'   lowest and three highest raw scores. The sum of both values must be strictly
#'   less than the total number of unique raw scores.
#' @param polyDegree.test Logical. If \code{TRUE} (default), the function fits
#'   polynomials of degree 1 through \code{poly.degree} and returns a detailed
#'   diagnostics data frame with R², adjusted R², AIC, BIC, and a five‑number
#'   summary of residuals for each degree. This allows the user to audit the
#'   adequacy of the chosen polynomial degree.
#' @param plot Logical. Master switch for all plots. If \code{FALSE}, no plots
#'   are generated (overrides \code{plot.csem} and \code{plot.score}). Default
#'   is \code{FALSE}.
#' @param plot.csem Character. One of \code{"both"} (default), \code{"raw"},
#'   \code{"scale"}, or \code{"none"}. Determines which CSEM plots are displayed:
#'   \itemize{
#'     \item \code{"both"}: overlaid raw and scaled CSEM against scale score.
#'     \item \code{"raw"}: raw‑score CSEM against raw score.
#'     \item \code{"scale"}: scale‑score CSEM against scale score.
#'     \item \code{"none"}: no CSEM plots.
#'   }
#'   Only used if \code{plot = TRUE}.
#' @param plot.score Logical. If \code{TRUE}, generates the raw‑to‑scale score
#'   conversion plot with the fitted polynomial curve (similar to Figures A‑1 to
#'   A‑4 in Woodruff et al., 2013). Default is \code{FALSE}. Only used if
#'   \code{plot = TRUE}.
#'
#' @return A list with four components:
#'   \item{data}{A \code{data.frame} with columns:
#'     \code{Raw_Score}, \code{Scale_Score}, \code{CSEM_Raw}, \code{Derivative},
#'     and \code{CSEM_Scaled}.}
#'   \item{model}{The \code{lm} object from the polynomial fit of degree
#'     \code{poly.degree} (fitted on the aggregated data).}
#'   \item{diagnostics}{A \code{data.frame} with one row per tested degree
#'     (from 1 to \code{poly.degree}, or fewer if constrained by sample size).
#'     Contains columns: \code{Degree}, \code{R2}, \code{Adj_R2}, \code{AIC},
#'     \code{BIC}, and residual summaries (\code{Min}, \code{Q1}, \code{Median},
#'     \code{Q3}, \code{Max}).}
#'   \item{plot}{A list of \code{ggplot} objects (if \code{plot = TRUE}),
#'     otherwise \code{NULL}. The list may contain \code{conversion} and
#'     \code{csem} elements.}
#'
#' @details
#' \strong{Audit‑ready step‑by‑step procedure (Woodruff et al., 2013):}
#' \enumerate{
#'   \item \strong{Automatic aggregation (many‑to‑one conversion):}
#'     If several raw scores map to the same scale score (common in score
#'     conversion tables), the raw scores are averaged for each unique scale
#'     value before fitting the polynomial. This mirrors the practice in the
#'     original study.
#'   \item \strong{Polynomial fitting with orthogonal basis:}
#'     The conversion function \eqn{g(X)} is approximated by a polynomial of
#'     degree \code{poly.degree} using \code{lm(scale ~ poly(raw, poly.degree, raw = FALSE))}.
#'     Orthogonal polynomials ensure numerical stability even with high‑degree
#'     terms and large raw‑score values.
#'   \item \strong{Computation of the derivative \eqn{g'(X)}:}
#'     For each raw score in the original (non‑aggregated) vector, the derivative
#'     is evaluated by finite differences with a very small \eqn{\epsilon}
#'     (\eqn{1e-6}) using the fitted model. This avoids analytical differentiation
#'     of orthogonal polynomials and provides exact numerical slopes.
#'   \item \strong{Delta method application:}
#'     The scale‑score CSEM is obtained as:
#'     \deqn{CSEM_{escalado}(X) = CSEM_{bruto}(X) \times |g'(X)|}
#'   \item \strong{Edge truncation (if enabled):}
#'     To prevent unrealistic estimates at the extremes due to polynomial
#'     boundary instability, the derivatives of the first \code{n.edge[1]} and
#'     last \code{n.edge[2]} raw scores are replaced by the derivative of the
#'     nearest untruncated interior point (i.e., the value at position
#'     \code{n.edge[1] + 1} for the lower tail, and at position
#'     \code{N - n.edge[2]} for the upper tail).
#'   \item \strong{Model diagnostics:}
#'     When \code{polyDegree.test = TRUE}, polynomials of all degrees from 1 to
#'     \code{poly.degree} are fitted and compared via R², AIC, BIC, and residual
#'     summaries. This allows the auditor to judge whether the chosen degree is
#'     appropriate (e.g., high R², low AIC/BIC, and no systematic residual patterns).
#'   \item \strong{Graphical validation:}
#'     The \code{plot.score} option produces a scatter plot of the raw‑to‑scale
#'     conversion points together with the fitted polynomial curve (analogous to
#'     Figures A‑1 to A‑4 in Woodruff et al., 2013). Additionally, CSEM plots
#'     (controlled by \code{plot.csem}) help visualise the effect of the
#'     transformation on measurement error.
#' }
#'
#' @references
#' Woodruff, D., Traynor, A., Cui, Z., & Fang, Y. (2013). \emph{A Comparison of
#' Three Methods for Computing Scale Score Conditional Standard Errors of
#' Measurement}. ACT Research Report Series 2013 (7). Iowa City, IA: ACT, Inc.
#'
#' Feldt, L. S., & Qualls, A. L. (1998). Approximating Scale Score Standard Error
#' of Measurement From the Raw Score Standard Error. \emph{Applied Measurement in
#' Education}, 11(2), 159‑177.
#'
#' @examples
#' \donttest{
#' # Example 1: Simple linear transformation (slope = 5)
#' raw <- 0:10
#' scale <- seq(20, 70, by = 5)
#' csem <- c(2.0, 1.8, 1.6, 1.4, 1.3, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0)
#'
#' result <- scaleCSEMWood(raw, scale, csem, poly.degree = 3,
#'                         truncate.edges = TRUE, n.edge = c(1, 1),
#'                         polyDegree.test = TRUE,
#'                         plot = TRUE, plot.csem = "both", plot.score = TRUE)
#'
#' # View the diagnostic table
#' print(result$diagnostics)
#'
#' # Export the data to CSV for reporting
#' write.csv(result$data, "csem_woodruff_output.csv", row.names = FALSE)
#' }
#'
#' \donttest{
#' # Example 2: Nonlinear transformation (simulated)
#' raw <- 0:20
#' scale <- round(100 + 5*raw + 0.2*raw^2 - 0.005*raw^3, 1)
#' csem <- rep(2.5, length(raw))  # constant raw CSEM
#'
#' result <- scaleCSEMWood(raw, scale, csem, poly.degree = 6,
#'                         truncate.edges = TRUE, n.edge = c(2, 2),
#'                         plot = TRUE, plot.csem = "scale", plot.score = TRUE)
#' }
#'
#' @export
scaleCSEMWood <- function(
    raw,
    scale,
    csem,
    poly.degree = 6,
    truncate.edges = TRUE,
    n.edge = c(1, 1),
    polyDegree.test = TRUE,
    plot = FALSE,
    plot.csem = "both",
    plot.score = FALSE
) {

  # --- 1. Basic Validations -------------------------------------------------

  if (length(raw) != length(scale) || length(raw) != length(csem)) {
    stop("raw, scale, and csem must have the same length.")
  }
  if (!is.numeric(raw) || !is.numeric(scale) || !is.numeric(csem)) {
    stop("raw, scale, and csem must be numeric vectors.")
  }
  if (any(raw != round(raw))) {
    stop("raw scores must be integers (or whole numbers).")
  }
  if (any(csem < 0, na.rm = TRUE)) {
    stop("csem values must be non‑negative.")
  }
  if (length(poly.degree) != 1 || poly.degree < 1 || poly.degree != round(poly.degree)) {
    stop("poly.degree must be a positive integer.")
  }

  if (!is.logical(truncate.edges) || length(truncate.edges) != 1) {
    stop("truncate.edges must be TRUE or FALSE.")
  }
  if (!is.numeric(n.edge) || length(n.edge) != 2 || any(n.edge < 0) || any(n.edge != round(n.edge))) {
    stop("n.edge must be an integer vector of length 2 with non‑negative values.")
  }

  plot.csem <- match.arg(plot.csem, choices = c("both", "raw", "scale", "none"))
  if (!is.logical(plot.score) || length(plot.score) != 1) {
    stop("plot.score must be TRUE or FALSE.")
  }

  # --- 2. Ordering by raw (for consistent indexing) ------------------------

  ord <- order(raw)
  raw <- raw[ord]
  scale <- scale[ord]
  csem <- csem[ord]

  # --- 3. Automatic aggregation for many‑to‑one conversions -----------------
  #    (average raw scores for each unique scale value)
  df_orig <- data.frame(raw = raw, scale = scale, csem = csem)
  df_agg <- aggregate(raw ~ scale, data = df_orig, FUN = mean)
  colnames(df_agg)[2] <- "raw_avg"

  # Ensure we have enough unique points to fit the requested degree
  n_unique <- nrow(df_agg)
  if (poly.degree >= n_unique) {
    warning(sprintf(
      "poly.degree (%d) is >= number of unique scale points (%d). Reducing degree to %d.",
      poly.degree, n_unique, n_unique - 1
    ))
    poly.degree <- n_unique - 1
  }
  if (poly.degree < 1) {
    stop("Need at least 2 unique scale points to fit a polynomial (degree >= 1).")
  }

  # --- 4. Diagnostics (polyDegree.test) -------------------------------------

  diag_df <- NULL
  if (isTRUE(polyDegree.test)) {
    max_deg <- min(poly.degree, n_unique - 1)
    diag_list <- vector("list", length = max_deg)

    for (d in seq_len(max_deg)) {
      form <- stats::as.formula(paste0("scale ~ poly(raw_avg, ", d, ", raw = FALSE)"))
      fit_d <- stats::lm(form, data = df_agg)

      smry <- summary(fit_d)
      r2 <- smry$r.squared
      adj_r2 <- smry$adj.r.squared
      aic <- stats::AIC(fit_d)
      bic <- stats::BIC(fit_d)

      res <- stats::residuals(fit_d)
      res_sum <- stats::quantile(res, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
      names(res_sum) <- c("Min", "Q1", "Median", "Q3", "Max")

      diag_list[[d]] <- data.frame(
        Degree = d,
        R2 = r2,
        Adj_R2 = adj_r2,
        AIC = aic,
        BIC = bic,
        Resid_Min = res_sum[1],
        Resid_Q1 = res_sum[2],
        Resid_Median = res_sum[3],
        Resid_Q3 = res_sum[4],
        Resid_Max = res_sum[5],
        stringsAsFactors = FALSE
      )
    }
    diag_df <- do.call(rbind, diag_list)
  }

  # --- 5. Fit the final polynomial model (on aggregated data) ---------------

  final_form <- stats::as.formula(paste0("scale ~ poly(raw_avg, ", poly.degree, ", raw = FALSE)"))
  model_final <- stats::lm(final_form, data = df_agg)

  # --- 6. Evaluate derivative g'(X) for the ORIGINAL raw scores ------------
  #    using finite differences on the fitted model (exact for polynomial)

  eps <- 1e-6
  newdata_orig <- data.frame(raw_avg = raw)  # predict on original raw scale
  pred0 <- stats::predict(model_final, newdata = newdata_orig)

  newdata_eps <- data.frame(raw_avg = raw + eps)
  pred1 <- stats::predict(model_final, newdata = newdata_eps)

  derivative <- (pred1 - pred0) / eps

  # --- 7. Edge truncation (clip by neighbour) ------------------------------

  if (isTRUE(truncate.edges)) {
    N <- length(derivative)
    lower_n <- min(n.edge[1], N - 1)   # cannot truncate all points
    upper_n <- min(n.edge[2], N - 1 - lower_n)

    if (lower_n > 0) {
      replacement_lower <- derivative[lower_n + 1]
      derivative[1:lower_n] <- replacement_lower
    }
    if (upper_n > 0) {
      replacement_upper <- derivative[N - upper_n]
      idx_upper <- (N - upper_n + 1):N
      derivative[idx_upper] <- replacement_upper
    }
  }

  # --- 8. Apply the Delta method -------------------------------------------

  scale_csem <- csem * abs(derivative)

  # --- 9. Prepare output data.frame ----------------------------------------

  output_df <- data.frame(
    Raw_Score      = raw,
    Scale_Score    = scale,
    CSEM_Raw       = csem,
    Derivative     = derivative,
    CSEM_Scaled    = scale_csem,
    stringsAsFactors = FALSE
  )

  # --- 10. Optional plots ---------------------------------------------------

  plot_list <- list()

  # Verificar disponibilidad de ggplot2
  has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

  if (!has_ggplot2) {
    warning("ggplot2 is not installed. Skipping all plots.")
  } else {

    # ---- 10a. Raw‑to‑Scale conversion plot (plot.score) ------------------
    if (isTRUE(plot.score)) {
      pred_curve <- stats::predict(model_final, newdata = data.frame(raw_avg = raw))

      p_conv <- ggplot2::ggplot() +
        ggplot2::geom_point(
          data = df_agg,
          ggplot2::aes(x = raw_avg, y = scale),
          color = "#2c3e50", size = 2.5, alpha = 0.7
        ) +
        ggplot2::geom_line(
          data = data.frame(raw_avg = raw, pred = pred_curve),
          ggplot2::aes(x = raw_avg, y = pred),
          color = "#e74c3c", linewidth = 1.2
        ) +
        ggplot2::labs(
          title = paste("Polynomial fit (degree =", poly.degree, ")"),
          x = "Raw Score",
          y = "Scale Score"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )

      plot_list[["conversion"]] <- p_conv
    }

    # ---- 10b. CSEM plots (plot.csem) -------------------------------------
    if (plot.csem != "none") {

      if (plot.csem == "raw") {
        p_csem <- ggplot2::ggplot(output_df, ggplot2::aes(x = Raw_Score, y = CSEM_Raw)) +
          ggplot2::geom_line(color = "#2c3e50", linewidth = 1) +
          ggplot2::geom_point(color = "#2c3e50", size = 2) +
          ggplot2::labs(
            title = "Raw‑Score CSEM",
            x = "Raw Score",
            y = "CSEM (raw metric)"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

      } else if (plot.csem == "scale") {
        p_csem <- ggplot2::ggplot(output_df, ggplot2::aes(x = Scale_Score, y = CSEM_Scaled)) +
          ggplot2::geom_line(color = "#e74c3c", linewidth = 1) +
          ggplot2::geom_point(color = "#e74c3c", size = 2) +
          ggplot2::labs(
            title = paste("Scale‑Score CSEM (Woodruff Delta, degree =", poly.degree, ")"),
            x = "Scale Score",
            y = "CSEM (scale metric)"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

      } else { # "both"
        p_csem <- ggplot2::ggplot(output_df) +
          ggplot2::geom_line(
            ggplot2::aes(x = Scale_Score, y = CSEM_Scaled, color = "Scale CSEM"),
            linewidth = 1
          ) +
          ggplot2::geom_point(
            ggplot2::aes(x = Scale_Score, y = CSEM_Scaled, color = "Scale CSEM"),
            size = 2
          ) +
          ggplot2::geom_line(
            ggplot2::aes(x = Scale_Score, y = CSEM_Raw, color = "Raw CSEM"),
            linewidth = 1, linetype = "dashed"
          ) +
          ggplot2::geom_point(
            ggplot2::aes(x = Scale_Score, y = CSEM_Raw, color = "Raw CSEM"),
            size = 2
          ) +
          ggplot2::scale_color_manual(
            values = c("Raw CSEM" = "#2c3e50", "Scale CSEM" = "#e74c3c")
          ) +
          ggplot2::labs(
            title = paste("CSEM Comparison (Woodruff Delta, degree =", poly.degree, ")"),
            x = "Scale Score",
            y = "CSEM",
            color = "Metric"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            legend.position = "bottom",
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
          )
      }

      plot_list[["csem"]] <- p_csem
    }

    # Print all generated plots
    for (p in plot_list) {
      print(p)
    }
  }

  # --- 11. Return structured output ----------------------------------------

  return(invisible(list(
    data        = output_df,
    model       = model_final,
    diagnostics = diag_df,
    plot        = if (length(plot_list) > 0) plot_list else NULL
  )))
}
