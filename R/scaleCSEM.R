#' Compute Conditional Standard Error of Measurement (CSEM) for Non-Linear Scale Scores
#'
#' @description
#' Applies the formal methodologies proposed by Feldt and Qualls (1998) to translate
#' raw score CSEMs into non-linear scale score metrics. Supports both the calculus-based
#' Polynomial Method (using monotonic splines) and the interval-based Approximation Method.
#'
#' @param raw Numeric vector. Raw scores (must be integer values, ideally consecutive).
#' @param scale Numeric vector. Scale scores corresponding to each raw score.
#' @param csem Numeric vector. Conditional standard errors of measurement in raw score units.
#' @param method Character. "approx" (interval method) or "polym" (monotonic spline method).
#' @param C Integer. Interval width for "approx". If NULL, uses round(1.5 * mean(csem)).
#' @param plot Logical. If TRUE, generates a plot.
#' @param plot.what Character. "both", "raw", or "scale" to choose what to display.
#'
#' @return A data.frame with columns: raw, scale, csem, slope, scale_csem.
#'
#' @details
#' \strong{Important technical notes:}
#' \enumerate{
#'   \item \strong{Use of SCAM in the "polym" method:} This method uses the \code{scam} package
#'         (Shape Constrained Additive Models) with a monotonic increasing P-spline basis
#'         (\code{bs = "mpi"}). This choice mathematically enforces a strictly non-negative
#'         first derivative across the entire raw score range, matching the psychometric
#'         requirement stated by Feldt & Qualls (1998, p. 163). Unconstrained splines
#'         (\code{smooth.spline}) or high-degree polynomials can produce negative slopes
#'         at the boundaries, leading to invalid negative scale score standard errors.
#'         The \code{scam} approach avoids arbitrary post-hoc adjustments.
#'         \emph{Note:} The \code{scam} package requires \code{mgcv}, which will be installed
#'         automatically when you install \code{scam}.
#'   \item \strong{Completeness of the conversion table (method "approx"):} This method
#'         directly looks up the scale scores associated with \eqn{X_0 \pm C} in the provided
#'         vectors. Therefore, it is mandatory that the raw vector contains all integer values
#'         required to cover \eqn{X_0 \pm C} for every row within the observed range of raw scores.
#'         Ideally, the table should include all integer raw scores from the minimum to the maximum
#'         observed. If any required raw score is missing, the function stops with an informative error.
#'         No interpolation is performed.
#' }
#'
#' \strong{Dependencies:}
#' The \code{"polym"} method requires the \code{scam} package (which in turn depends on \code{mgcv}).
#' Please ensure both are installed:
#'   \code{install.packages(c("scam", "mgcv"))}.
#'   If these packages are not available, use \code{method = "approx"} instead.
#'
#' \strong{Range of raw scores:}
#'   The function does not assume that raw scores start at 0. It uses the minimum and maximum values
#'   present in the \code{raw} vector as the natural boundaries of the score scale. This makes it
#'   suitable for Likert-type scales (e.g., summing 5 items each scored 1-5 gives a minimum of 5).
#'   Internally, the intervals \eqn{X_0 \pm C} are truncated to the observed range \eqn{[min(raw), max(raw)]}.
#'
#' @references
#' Feldt, L. S., & Qualls, A. L. (1998). Approximating Scale Score Standard Error of
#' Measurement From the Raw Score Standard Error. Applied Measurement in Education, 11(2), 159-177.
#'
#' @examples
#'
#' \donttest{
#' # Example with linear transformation (slope = 5)
#' raw <- 0:10
#' scale <- seq(20, 70, by = 5)
#' csem <- c(2.0, 1.8, 1.6, 1.4, 1.3, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0)
#'
#' scaleCSEM(raw, scale, csem, method = "approx", plot = TRUE)
#' }
#'
#' \donttest{
#'  # Simulate a nonlinear scale (example: square root)
#' raw <- 0:10
#' scale <- round(20 + 30 * sqrt(raw/10))  # no lineal
#' csem <- c(2.0, 1.8, 1.6, 1.4, 1.3, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0)
#'
#' scaleCSEM(raw, scale, csem, method = "approx", plot = TRUE)
#'}
#'
#' \donttest{
#' # Full workflow
#'
#' # Loading data
#' data("bfi")
#'
#' # Choosing variables
#' data.bfi <- bfi[, c("N1", "N2", "N3", "N4", "N5", "gender", "age")]
#'
#' # Clean for missing values
#' data.bfi.nmiss <- data.bfi[complete.cases(data.bfi), ]
#'
#' # CSEM with bootstrapping
#' output.boots1 <- csemBoots(data = data.bfi.nmiss[,1:5], ci = F,
#' conf.level = .90,
#' full.range = T,
#' score.range = c(5, 30),
#' smooth = T, B = 2000)
#'
#' # Score sum
#' uno <- table(rowSums(data.bfi.nmiss[,1:5]))
#'
#' # t Scores, linear transformation
#' dos <- table(round(psych::rescale(x = rowSums(data.bfi.nmiss[,1:5]), mean = 50, sd = 10)))
#'
#' # merge ot dataframe
#' dframe <- cbind.data.frame(score = as.data.frame(uno)$Var1,
#' tScore = as.data.frame(dos)[1])
#'
#' colnames(dframe)[2] <- "tscore"
#'
#' scaleCSEM(raw = output.boots1$CSEM$score,
#' scale = as.numeric(dframe$tscore),
#' csem = output.boots1$CSEM$CSEM.smooth,
#' method = "polym", C = 5,plot = T, plot.what = "both")
#' }
#'
#' @export
scaleCSEM <- function(raw, scale, csem,
                      method = c("approx", "polym"),
                      C = NULL,
                      plot = FALSE,
                      plot.what = "both") {

  # --- 1. Basic Validations ---
  method <- match.arg(method)

  if (length(raw) != length(scale) || length(raw) != length(csem)) {
    stop("raw, scale, and csem must have the same length.")
  }
  if (!is.numeric(raw) || !is.numeric(scale) || !is.numeric(csem)) {
    stop("raw, scale, and csem must be numeric vectors.")
  }
  if (any(raw != round(raw))) {
    stop("raw scores must be integers (or whole numbers).")
  }

  # Sort by raw (important for consistency)
  ord <- order(raw)
  raw <- raw[ord]
  scale <- scale[ord]
  csem <- csem[ord]

  # Observed range (we do not assume it starts at 0)
  raw_min <- min(raw)
  raw_max <- max(raw)

  # --- 2. Set C to the "approx" method ---
  if (method == "approx") {
    if (is.null(C)) {
      C <- round(1.5 * mean(csem, na.rm = TRUE))
      C <- max(C, 1)   # minimum 1
    } else {
      if (!is.numeric(C) || length(C) != 1 || C < 1 || C != round(C)) {
        stop("C must be a positive integer (or NULL for automatic calculation).")
      }
      C <- as.integer(C)
    }

    # Calculate the truncated L and U statistics based on the observed range
    L_vals <- pmax(raw - C, raw_min)
    U_vals <- pmin(raw + C, raw_max)

    # Verify that all required values are present in raw
    needed <- unique(c(L_vals, U_vals))
    missing <- setdiff(needed, raw)
    if (length(missing) > 0) {
      missing_str <- paste(sort(missing), collapse = ", ")
      stop(sprintf(
        "The raw vector is missing the following values required for the intervals: %s.
        Please provide a complete conversion table (all integer raw scores from %d to %d).",
        missing_str, raw_min, raw_max
      ))
    }

    # Calculate slopes and scale_csem
    idx_L <- match(L_vals, raw)
    idx_U <- match(U_vals, raw)
    scale_L <- scale[idx_L]
    scale_U <- scale[idx_U]
    denom <- U_vals - L_vals   # actual interval (not always 2*C at the ends)
    slope <- (scale_U - scale_L) / denom
    scale_csem <- csem * slope

    output <- data.frame(raw = raw, scale = scale, csem = csem,
                         slope = slope, scale_csem = scale_csem)
  }

  # --- 3. "polym" method (monotonic spline) ---
  if (method == "polym") {
    if (!requireNamespace("scam", quietly = TRUE)) {
      stop("Package 'scam' is required for method 'polym'. Please install it (this will also install 'mgcv').")
    }

    df <- data.frame(raw = raw, scale = scale)
    scam_formula <- as.formula("scale ~ s(raw, bs = 'mpi')")

    scam_model <- tryCatch(
      scam::scam(scam_formula, data = df),
      error = function(e) {
        stop("scam fitting failed. Possibly too few data points or non-monotonic relationship. Try method = 'approx'.\n",
             "Original error: ", e$message)
      }
    )

    # Numerical derivative
    eps <- 1e-5
    pred0 <- predict(scam_model, newdata = df)
    df_eps <- df
    df_eps$raw <- df_eps$raw + eps
    pred1 <- predict(scam_model, newdata = df_eps)
    slope <- (pred1 - pred0) / eps
    slope <- pmax(slope, 0)   # security (should not be negative)
    scale_csem <- csem * slope

    output <- data.frame(raw = raw, scale = scale, csem = csem,
                         slope = slope, scale_csem = scale_csem)
  }

  # --- 4. Optional display ---
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 not installed. Skipping plot.")
    } else {

      method_label <- if (method == "polym") {
        "Monotonic Spline Method"
      } else {
        paste0("Approximation Method [C = ", C, "]")
      }

      if (plot.what == "raw") {
        p <- ggplot2::ggplot(output, ggplot2::aes(x = raw, y = csem)) +
          ggplot2::geom_line(color = "#2c3e50", linewidth = 1) +
          ggplot2::geom_point(color = "#2c3e50", size = 2) +
          ggplot2::labs(title = "CSEM: Raw Score Metric",
                        x = "Raw Score", y = "Raw CSEM") +
          ggplot2::theme_minimal()
        print(p)

      } else if (plot.what == "scale") {
        p <- ggplot2::ggplot(output, ggplot2::aes(x = scale, y = scale_csem)) +
          ggplot2::geom_line(color = "#e74c3c", linewidth = 1) +
          ggplot2::geom_point(color = "#e74c3c", size = 2) +
          ggplot2::labs(title = paste("CSEM:", method_label),
                        x = "Scale Score", y = "Scale CSEM") +
          ggplot2::theme_minimal()
        print(p)

      } else {
        p <- ggplot2::ggplot(output) +
          ggplot2::geom_line(ggplot2::aes(x = scale, y = scale_csem, color = "Scale CSEM"),
                             linewidth = 1) +
          ggplot2::geom_point(ggplot2::aes(x = scale, y = scale_csem, color = "Scale CSEM"),
                              size = 2) +
          ggplot2::geom_line(ggplot2::aes(x = scale, y = csem, color = "Raw CSEM"),
                             linewidth = 1, linetype = "dashed") +
          ggplot2::geom_point(ggplot2::aes(x = scale, y = csem, color = "Raw CSEM"),
                              size = 2) +
          ggplot2::scale_color_manual(values = c("Raw CSEM" = "#2c3e50",
                                                 "Scale CSEM" = "#e74c3c")) +
          ggplot2::labs(title = paste("CSEM Comparison -", method_label),
                        x = "Scale Score", y = "CSEM",
                        color = "Metric") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "bottom")
        print(p)
      }
    }
  }

  # --- 5. Output ---
  return(output)
}
