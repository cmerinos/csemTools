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
#'   \item \strong{Completeness of the conversion table (method "approx"):} This method
#'         directly looks up the scale scores associated with \eqn{X_0 \pm C} in the provided
#'         vectors. Therefore, it is mandatory that the raw vector contains all integer values
#'         required to cover \eqn{X_0 \pm C} for every row. Ideally, the table should include
#'         all integer raw scores from the minimum to the maximum observed in the data.
#'         If any required raw score is missing, the function stops with an informative error.
#'         No interpolation is performed.
#' }
#'
#' @references
#' Feldt, L. S., & Qualls, A. L. (1998). Approximating Scale Score Standard Error of
#' Measurement From the Raw Score Standard Error. Applied Measurement in Education, 11(2), 159-177.
#'
#' @examples
#' raw <- 0:10
#' scale <- seq(20, 70, by = 5)
#' csem <- c(2.0, 1.8, 1.6, 1.4, 1.3, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0)
#' result <- scaleCSEM(raw, scale, csem, method = "approx", plot = TRUE)
#' head(result)
#'
#' @export
scaleCSEM <- function(raw, scale, csem,
                      method = c("approx", "polym"),
                      C = NULL,
                      plot = FALSE,
                      plot.what = "both") {

  # --- 1. Validaciones básicas ---
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

  # Ordenar por raw (importante para consistencia)
  ord <- order(raw)
  raw <- raw[ord]
  scale <- scale[ord]
  csem <- csem[ord]

  # Determinar rango real de raw (no asumir que empieza en 0)
  min_raw <- min(raw)
  max_raw <- max(raw)

  # --- 2. Configurar C para método "approx" ---
  if (method == "approx") {
    if (is.null(C)) {
      C <- round(1.5 * mean(csem, na.rm = TRUE))
      C <- max(C, 1)   # mínimo 1
    } else {
      if (!is.numeric(C) || length(C) != 1 || C < 1 || C != round(C)) {
        stop("C must be a positive integer (or NULL for automatic calculation).")
      }
      C <- as.integer(C)
    }

    # Validar que para cada raw, los valores X0±C existan en el vector raw
    L_vals <- pmax(raw - C, min_raw)
    U_vals <- pmin(raw + C, max_raw)
    needed <- unique(c(L_vals, U_vals))
    missing <- setdiff(needed, raw)
    if (length(missing) > 0) {
      missing_str <- paste(sort(missing), collapse = ", ")
      stop(sprintf(
        "The raw vector is missing the following values required for the intervals: %s.
        Please provide a complete conversion table for all integer raw scores from %d to %d.",
        missing_str, min_raw, max_raw
      ))
    }

    # Calcular slopes y scale_csem de forma vectorizada
    idx_L <- match(L_vals, raw)
    idx_U <- match(U_vals, raw)
    scale_L <- scale[idx_L]
    scale_U <- scale[idx_U]
    denom <- U_vals - L_vals   # 2*C generalmente, pero ajustado en bordes
    slope <- (scale_U - scale_L) / denom
    scale_csem <- csem * slope

    output <- data.frame(raw = raw, scale = scale, csem = csem,
                         slope = slope, scale_csem = scale_csem)
  }

  # --- 3. Método "polym" (monotonic spline) ---
  if (method == "polym") {
    if (!requireNamespace("scam", quietly = TRUE)) {
      stop("Package 'scam' is required for method 'polym'. Please install it.")
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

    # Derivada numérica usando diferencias finitas
    eps <- 1e-5
    pred0 <- predict(scam_model, newdata = df)
    df_eps <- df
    df_eps$raw <- df_eps$raw + eps
    pred1 <- predict(scam_model, newdata = df_eps)
    slope <- (pred1 - pred0) / eps
    slope <- pmax(slope, 0)   # seguridad
    scale_csem <- csem * slope

    output <- data.frame(raw = raw, scale = scale, csem = csem,
                         slope = slope, scale_csem = scale_csem)
  }

  # --- 4. Visualización (opcional) ---
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 not installed. Skipping plot.")
    } else {
      library(ggplot2)

      method_label <- if (method == "polym") {
        "Monotonic Spline Method"
      } else {
        paste0("Approximation Method [C = ", C, "]")
      }

      if (plot.what == "raw") {
        p <- ggplot(output, aes(x = raw, y = csem)) +
          geom_line(color = "#2c3e50", linewidth = 1) +
          geom_point(color = "#2c3e50", size = 2) +
          labs(title = "CSEM: Raw Score Metric",
               x = "Raw Score", y = "Raw CSEM") +
          theme_minimal()
        print(p)
      } else if (plot.what == "scale") {
        p <- ggplot(output, aes(x = scale, y = scale_csem)) +
          geom_line(color = "#e74c3c", linewidth = 1) +
          geom_point(color = "#e74c3c", size = 2) +
          labs(title = paste("CSEM:", method_label),
               x = "Scale Score", y = "Scale CSEM") +
          theme_minimal()
        print(p)
      } else {  # "both" por defecto
        p <- ggplot(output) +
          geom_line(aes(x = scale, y = scale_csem, color = "Scale CSEM"),
                    linewidth = 1) +
          geom_point(aes(x = scale, y = scale_csem, color = "Scale CSEM"),
                     size = 2) +
          geom_line(aes(x = scale, y = csem, color = "Raw CSEM"),
                    linewidth = 1, linetype = "dashed") +
          geom_point(aes(x = scale, y = csem, color = "Raw CSEM"),
                     size = 2) +
          scale_color_manual(values = c("Raw CSEM" = "#2c3e50",
                                        "Scale CSEM" = "#e74c3c")) +
          labs(title = paste("CSEM Comparison -", method_label),
               x = "Scale Score", y = "CSEM",
               color = "Metric") +
          theme_minimal() +
          theme(legend.position = "bottom")
        print(p)
      }
    }
  }

  # --- 5. Retornar resultado ---
  return(output)
}
