#' Compute Conditional Standard Error of Measurement (CSEM) for Non-Linear Scale Scores
#'
#' @description
#' Applies the formal methodologies proposed by Feldt and Qualls (1998) to translate
#' raw score CSEMs into non-linear scale score metrics. Supports both the calculus-based
#' Polynomial Method (evaluating exact local derivatives under monotonicity constraints)
#' and the interval-based Approximation Method (employing an empirical symmetric window of 1.5 * Mean CSEM).
#'
#' @details
#' \enumerate{
#'   \item \strong{Polynomial Method (method = "polym"):} Fits a Shape Constrained Additive Model
#'   (SCAM) with a monotonic increasing P-spline (\code{bs = "mpi"}) to the raw-to-scale conversion.
#'   This guarantees a strictly non-negative first derivative across the entire score range,
#'   avoiding arbitrary post-hoc absolute value adjustments. The derivative is evaluated
#'   numerically at each raw score point and multiplied by the raw CSEM.
#'
#'   \item \strong{Approximation Method (method = "approx"):} Uses a symmetric raw score window
#'   of width \eqn{2C} around each raw score \eqn{X_0}. The slope is approximated as
#'   \eqn{(Scale(X_0+C) - Scale(X_0-C)) / (2C)}, with boundary corrections when the window
#'   exceeds the valid raw score range (0 to k). The scale CSEM is then \eqn{slope * rawCSEM}.
#'   \strong{Important:} This method requires that the conversion table contains all raw scores
#'   needed for the intervals \eqn{X_0 \pm C}. No interpolation is performed; the function will
#'   stop with an error if any required raw score is missing.
#' }
#'
#' @param data Optional data frame containing the raw-to-scale conversion table and CSEMs.
#'        If provided, \code{raw.col}, \code{scale.col}, and \code{rawcsem.col} must be specified.
#' @param raw.col Character string. Name of the column containing raw scores (if \code{data} is used).
#' @param scale.col Character string. Name of the column containing transformed scale scores.
#' @param rawcsem.col Character string. Name of the column containing raw score CSEMs.
#' @param raw Numeric vector. Alternative to \code{data}: raw scores.
#' @param scale Numeric vector. Alternative to \code{data}: scale scores.
#' @param csem Numeric vector. Alternative to \code{data}: raw CSEMs.
#' @param method Character string. Either \code{"approx"} (default) or \code{"polym"}.
#' @param C Integer. Window width for the approximation method. If \code{NULL} (default),
#'        computed as \code{round(1.5 * mean(csem))}, with a minimum of 1.
#' @param plot Logical. If \code{TRUE}, generates a plot of the CSEMs. Default \code{FALSE}.
#' @param plot.what Character string. One of \code{"both"} (overlay raw and scale CSEMs),
#'        \code{"scale"} (only scale CSEM), or \code{"raw"} (only raw CSEM). Default \code{"both"}.
#'
#' @return A data frame (sorted by raw score) with the original columns plus:
#' \item{slope}{The estimated slope (derivative or interval slope) at each raw score.}
#' \item{scale_csem}{The CSEM in the scale score metric.}
#'
#' @references
#' Feldt, L. S., & Qualls, A. L. (1998). Approximating Scale Score Standard Error of
#' Measurement From the Raw Score Standard Error. \emph{Applied Measurement in Education},
#' 11(2), 159-177. \doi{10.1177/0013164499591001}
#'
#' @examples
#' # Example with a complete conversion table (raw scores 0 to 20)
#' set.seed(123)
#' raw_all <- 0:20
#' scale_all <- 50 + 10 * scale(raw_all)  # just an example
#' csem_all <- 2 + 0.1 * abs(raw_all - 10) # example pattern
#' df_full <- data.frame(raw = raw_all, scale = scale_all, csem = csem_all)
#'
#' # Approximation method
#' res_approx <- scaleCSEM(data = df_full, raw.col = "raw", scale.col = "scale",
#'                         rawcsem.col = "csem", method = "approx", C = 3)
#' head(res_approx)
#'
#' # Polynomial method (requires scam package)
#' if (requireNamespace("scam", quietly = TRUE)) {
#'   res_polym <- scaleCSEM(data = df_full, raw.col = "raw", scale.col = "scale",
#'                          rawcsem.col = "csem", method = "polym", plot = TRUE)
#'   head(res_polym)
#' }
#'
#' # Using direct vectors instead of data frame
#' res_vec <- scaleCSEM(raw = raw_all, scale = scale_all, csem = csem_all,
#'                      method = "approx", C = 3)
#'
#' @export
scaleCSEM <- function(data = NULL,
                      raw.col = NULL, scale.col = NULL, rawcsem.col = NULL,
                      raw = NULL, scale = NULL, csem = NULL,
                      method = c("approx", "polym"),
                      C = NULL,
                      plot = FALSE, plot.what = "both") {

  # ---- 1. Input validation and data extraction ----
  method <- match.arg(method)

  # Check for required packages
  if (method == "polym" && !requireNamespace("scam", quietly = TRUE)) {
    stop("Method 'polym' requires the 'scam' package. Please install it.")
  }
  if (plot && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Plotting requires the 'ggplot2' package. Please install it.")
  }

  # Extract vectors from either data frame or direct arguments
  if (!is.null(data)) {
    # Using data frame interface
    if (is.null(raw.col) || is.null(scale.col) || is.null(rawcsem.col)) {
      stop("When 'data' is provided, 'raw.col', 'scale.col', and 'rawcsem.col' must be specified.")
    }
    if (!all(c(raw.col, scale.col, rawcsem.col) %in% names(data))) {
      stop("One or more specified columns not found in 'data'.")
    }
    # Order by raw score and extract
    data <- data[order(data[[raw.col]]), ]
    raw_vec <- data[[raw.col]]
    scale_vec <- data[[scale.col]]
    csem_vec <- data[[rawcsem.col]]
    # Keep the original data for output
    original_data <- data
  } else {
    # Using direct vectors
    if (is.null(raw) || is.null(scale) || is.null(csem)) {
      stop("Either 'data' with column names, or 'raw', 'scale', and 'csem' vectors must be provided.")
    }
    if (length(raw) != length(scale) || length(raw) != length(csem)) {
      stop("'raw', 'scale', and 'csem' must have the same length.")
    }
    # Order by raw score
    ord <- order(raw)
    raw_vec <- raw[ord]
    scale_vec <- scale[ord]
    csem_vec <- csem[ord]
    # Build a data frame for output (preserve names)
    original_data <- data.frame(raw = raw_vec, scale = scale_vec, csem = csem_vec,
                                stringsAsFactors = FALSE)
    # Set column names for later use
    raw.col <- "raw"
    scale.col <- "scale"
    rawcsem.col <- "csem"
  }

  # Check for duplicates in raw_vec
  if (any(duplicated(raw_vec))) {
    stop("Raw scores must be unique. Please remove duplicates.")
  }

  k <- max(raw_vec)  # maximum raw score (test length)
  n <- length(raw_vec)

  # Ensure raw_vec is integer (or at least numeric)
  if (!is.numeric(raw_vec)) stop("Raw scores must be numeric.")

  # ---- 2. Compute slopes and scale CSEM ----

  if (method == "polym") {
    # ---- 2a. Polynomial Method (Monotonic Spline) ----
    # Fit a monotonic increasing P-spline
    # Use scam with bs = "mpi" (monotonic increasing P-spline)
    scam_form <- as.formula(paste(scale.col, "~ s(", raw.col, ", bs = 'mpi')"))
    # We need to build a temporary data frame for scam
    temp_df <- data.frame(x = raw_vec, y = scale_vec)
    names(temp_df) <- c(raw.col, scale.col)
    scam_model <- scam::scam(scam_form, data = temp_df)

    # Numerical derivative via finite differences
    eps <- 1e-5
    x_plus <- raw_vec + eps
    pred_x <- predict(scam_model, newdata = setNames(data.frame(raw_vec), raw.col))
    pred_x_plus <- predict(scam_model, newdata = setNames(data.frame(x_plus), raw.col))
    slope <- (pred_x_plus - pred_x) / eps

    # Ensure slope is non-negative (should be, but safeguard)
    slope <- pmax(slope, 0)

    # Compute scale CSEM
    scale_csem <- csem_vec * slope

  } else {
    # ---- 2b. Approximation Method ----
    # Determine C if not provided
    if (is.null(C)) {
      C <- round(1.5 * mean(csem_vec, na.rm = TRUE))
      C <- max(C, 1)   # at least 1
    } else {
      C <- as.integer(C)
      if (C < 1) stop("C must be a positive integer.")
    }

    # Pre-allocate vectors
    slope <- numeric(n)
    scale_csem <- numeric(n)

    # For each raw score, compute L and U, look up scale values
    # We'll use a named vector for fast lookup (scale by raw)
    scale_lookup <- setNames(scale_vec, raw_vec)

    for (i in seq_len(n)) {
      X0 <- raw_vec[i]
      L <- max(X0 - C, 0)
      U <- min(X0 + C, k)

      # Check that L and U exist in raw_vec
      if (!(L %in% raw_vec)) {
        stop(sprintf("At raw score = %g, the lower bound L = %g is not present in the raw score table. Please provide a complete table (all raw scores from 0 to %g).", X0, L, k))
      }
      if (!(U %in% raw_vec)) {
        stop(sprintf("At raw score = %g, the upper bound U = %g is not present in the raw score table. Please provide a complete table (all raw scores from 0 to %g).", X0, U, k))
      }

      scale_L <- scale_lookup[as.character(L)]
      scale_U <- scale_lookup[as.character(U)]

      # Compute slope according to boundary cases
      denom <- U - L   # this is 2C for interior, but may differ at boundaries
      # (Note: if X0-C < 0, denom = X0 + C; if X0+C > k, denom = k - X0 + C)
      # But U-L already gives the correct denominator by definition.
      if (denom == 0) {
        stop("Denominator for slope calculation is zero. Check C and raw score table.")
      }
      slope[i] <- (scale_U - scale_L) / denom
      # Slope should be non-negative (scale is monotonic)
      if (slope[i] < 0) {
        warning(sprintf("Negative slope detected at raw score = %g. Setting to zero.", X0))
        slope[i] <- 0
      }
      scale_csem[i] <- csem_vec[i] * slope[i]
    }
  }

  # ---- 3. Build output data frame ----
  # Start with original_data (already sorted by raw)
  output_df <- original_data
  output_df$slope <- slope
  output_df$scale_csem <- scale_csem

  # Restore column names if they were not the defaults
  # We'll rename to match user's original column names if data was provided
  if (!is.null(data)) {
    # Keep original column names; output_df already has them
    # But we added 'slope' and 'scale_csem'
    # If data was provided, we should preserve its column names.
    # Actually original_data was data sorted, so column names are as user provided.
    # So we just add the new columns.
    # We'll keep as is.
  } else {
    # If vectors were provided, we already named columns as raw, scale, csem
    # So output_df has columns: raw, scale, csem, slope, scale_csem
    # That's fine.
  }

  # ---- 4. Visualization ----
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 not available. Plot skipped.")
    } else {
      library(ggplot2)  # for aes_string, etc.

      # Determine labels based on method
      method_label <- if (method == "polym") {
        "Monotonic Spline Method"
      } else {
        paste0("Approximation Method [C = ", C, "]")
      }

      # Determine which plot type
      if (plot.what == "raw") {
        p <- ggplot(output_df, aes_string(x = raw.col, y = rawcsem.col)) +
          geom_line(color = "#2c3e50", size = 1) +
          geom_point(color = "#2c3e50", size = 2) +
          labs(title = "Conditional Standard Error: Raw Metric",
               x = "Raw Score", y = "Raw CSEM") +
          theme_minimal()
        print(p)
      } else if (plot.what == "scale") {
        p <- ggplot(output_df, aes_string(x = scale.col, y = "scale_csem")) +
          geom_line(color = "#e74c3c", size = 1) +
          geom_point(color = "#e74c3c", size = 2) +
          labs(title = paste("CSEM (Scale Metric):", method_label),
               x = "Scale Score", y = "Scale CSEM") +
          theme_minimal()
        print(p)
      } else {  # "both"
        # Overlay raw and scale CSEMs, both mapped to scale score on x-axis
        # We'll use the scale score as x, and plot both raw CSEM and scale CSEM vs scale score
        p <- ggplot(output_df) +
          geom_line(aes_string(x = scale.col, y = "scale_csem",
                               color = paste0("Scale CSEM (", method_label, ")"))) +
          geom_point(aes_string(x = scale.col, y = "scale_csem",
                                color = paste0("Scale CSEM (", method_label, ")"))) +
          geom_line(aes_string(x = scale.col, y = rawcsem.col,
                               color = "Raw CSEM"), linetype = "dashed") +
          geom_point(aes_string(x = scale.col, y = rawcsem.col,
                                color = "Raw CSEM")) +
          scale_color_manual(values = c("Raw CSEM" = "#2c3e50",
                                        paste0("Scale CSEM (", method_label, ")") = "#e74c3c")) +
          labs(title = "Overlay of Raw and Scale CSEMs",
               x = "Scale Score", y = "CSEM",
               color = "Metric") +
          theme_minimal() +
          theme(legend.position = "bottom")
        print(p)
      }
    }
  }

  # ---- 5. Return ----
  return(output_df)
}
