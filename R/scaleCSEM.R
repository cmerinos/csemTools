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
#'   \item \strong{Monotone Spline Smoothing (method = "polym"):} Instead of unconstrained high-degree
#'   polynomials which suffer from boundary oscillations (Runge's phenomenon), this function fits a
#'   Shape Constrained Additive Model (SCAM) forcing a strictly monotonic increasing relationship. This maps
#'   to the authors' premise that the conversion function derivative must never be negative, avoiding
#'   arbitrary post-hoc absolute value adjustments.
#'   \item \strong{Exact Derivative Evaluation (method = "polym"):} Computes the point-specific instantaneous
#'   conversion slope by extracting the numerical first derivative of the fitted monotonic smoothing spline.
#'   \item \strong{Empirical Window Width (method = "approx"):} Following page 165 of Feldt & Qualls (1998),
#'   if the parameter \code{C} is omitted, it is determined dynamically as \eqn{C = 1.5 \times \text{mean}(\sigma_{E(X)})},
#'   rounded to the nearest integer. This aligns with the psychometric properties of the target test
#'   rather than relying on historic test-length heuristics. Boundary corrections follow Equations 5 and 6 exactly.
#' }
#'
#' @param data A data frame containing the raw-to-scale score conversion table and baseline errors.
#' @param raw.col Character string. The name of the column containing raw scores.
#' @param scale.col Character string. The name of the column containing transformed scale scores.
#' @param rawcsem.col Character string. The name of the column containing pre-computed raw score CSEMs.
#' @param method Character string. The psychometric approach to estimate conversion slopes.
#' Options are "approx" (Interval Approximation Method) or "polym" (Monotonic Smoothing Method). Default is "approx".
#' @param C Integer. The interval constant for the raw score window (X0 +/- C) when method = "approx". If NULL, it automatically computes \code{round(1.5 * mean(rawcsem))}.
#' @param plot Logical. If TRUE, generates a graphical visualization of the conditional errors. Default is FALSE.
#' @param plot.what Character string. Layout options: "both" (overlayed metrics), "scale" (scale CSEM only), or "raw" (raw CSEM only). Default is "both".
#'
#' @return A data frame sorted by raw score with appended columns:
#' \item{slope}{The exact conversion derivative or local interval slope computed via specified method.}
#' \item{scale_csem}{The approximated Conditional Standard Error of Measurement mapped onto the scale score metric.}
#'
#' @references
#' Feldt, L. S., & Qualls, A. L. (1998). Approximating Scale Score Standard Error of
#' Measurement From the Raw Score Standard Error. \emph{Applied Measurement in Education},
#' 11(2), 159-177. https://doi.org/10.1177/0013164499591001
#'
#' @export
scaleCSEM <- function(data, raw.col, scale.col, rawcsem.col,
                      method = c("approx", "polym"), C = NULL,
                      plot = FALSE, plot.what = "both") {

  # 1. Dependency and Argument Validation
  method <- match.arg(method)
  required_packages <- c("dplyr", "ggplot2", "scam")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(paste("Please install required packages:", paste(missing_packages, collapse = ", ")))
  }

  library(dplyr)
  library(ggplot2)
  library(scam)

  if (!all(c(raw.col, scale.col, rawcsem.col) %in% names(data))) {
    stop("One or more specified columns do not exist in the data frame.")
  }

  # Sort baseline data by raw score internally
  working_data <- data %>% arrange(!!sym(raw.col))
  raw_vec <- working_data[[raw.col]]
  scale_vec <- working_data[[scale.col]]
  rawcsem_vec <- working_data[[rawcsem.col]]
  k <- max(raw_vec)

  # 2. Compute Slopes Based on Refined Psychometric Methods
  if (method == "polym") {
    # --- MONOTONIC SMOOTHING METHOD (Feldt & Qualls, 1998, p. 162-164) ---
    # Fit a shape-constrained additive model forcing a monotonically increasing shape (bs = "mpi")
    scam_model <- scam::scam(scale_vec ~ s(raw_vec, bs = "mpi"))

    # Calculate numerical first derivative via finite differences of the smooth monotonic function
    eps <- 1e-5
    pred_x <- predict(scam_model, newdata = data.frame(raw_vec = raw_vec))
    pred_x_eps <- predict(scam_model, newdata = data.frame(raw_vec = raw_vec + eps))
    slopes_vec <- (pred_x_eps - pred_x) / eps

    output_data <- working_data %>%
      mutate(
        slope = slopes_vec,
        scale_csem = !!sym(rawcsem.col) * slope
      )

  } else {
    # --- APPROXIMATION METHOD (Feldt & Qualls, 1998, p. 164-166) ---
    # Dynamically define C based on Feldt & Qualls' true 1.5 * Mean CSEM recommendation (p. 165)
    if (is.null(C)) {
      mean_raw_csem <- mean(rawcsem_vec, na.rm = TRUE)
      C <- round(1.5 * mean_raw_csem)
      C <- max(C, 1) # Enforce a minimum window size of 1
    }

    output_data <- working_data %>%
      rowwise() %>%
      mutate(
        .X0 = !!sym(raw.col),
        .U_raw = min(.X0 + C, k),
        .L_raw = max(.X0 - C, 0),
        .GE_U = scale_vec[which(raw_vec == .U_raw)],
        .GE_L = scale_vec[which(raw_vec == .L_raw)],

        # Boundary adjusted interval slope parameters via Equations 4, 5, and 6
        slope = case_when(
          (.X0 + C) > k  ~ abs(.GE_U - .GE_L) / (k - .X0 + C), # Equation 5 (Upper tail)
          (.X0 - C) < 0  ~ abs(.GE_U - .GE_L) / (.X0 + C),     # Equation 6 (Lower tail)
          TRUE           ~ abs(.GE_U - .GE_L) / (2 * C)        # Equation 4 (Standard)
        ),
        scale_csem = !!sym(rawcsem.col) * slope
      ) %>%
      ungroup() %>%
      select(-starts_with("."))
  }

  # 3. Automated Visualization Module
  if (plot) {
    method_label <- if(method == "polym") {
      "Monotonic Spline Method"
    } else {
      paste0("Approximation Method [C = ", C, "]")
    }

    if (plot.what == "raw") {
      p <- ggplot(output_data, aes(x = !!sym(raw.col), y = !!sym(rawcsem.col))) +
        geom_line(color = "#2c3e50", linewidth = 1) + geom_point(color = "#2c3e50", size = 2) +
        labs(title = "Conditional Standard Error of Measurement: Raw Metric", x = "Raw Score", y = "Raw CSEM") +
        theme_minimal()
      print(p)
    } else if (plot.what == "scale") {
      p <- ggplot(output_data, aes(x = !!sym(scale.col), y = scale_csem)) +
        geom_line(color = "#e74c3c", linewidth = 1) + geom_point(color = "#e74c3c", size = 2) +
        labs(title = paste("CSEM:", method_label), x = "Scale Score", y = "Scale CSEM") +
        theme_minimal()
      print(p)
    } else if (plot.what == "both") {
      raw_vals <- output_data[[raw.col]]
      scale_vals <- output_data[[scale.col]]

      p <- ggplot(output_data) +
        geom_line(aes(x = !!sym(scale.col), y = scale_csem, color = paste("Scale Metric (", method_label, ")")), linewidth = 1) +
        geom_point(aes(x = !!sym(scale.col), y = scale_csem, color = paste("Scale Metric (", method_label, ")")), size = 2) +
        geom_line(aes(x = !!sym(scale.col), y = !!sym(rawcsem.col), color = "Raw Score Metric"), linewidth = 1, linetype = "dashed") +
        geom_point(aes(x = !!sym(scale.col), y = !!sym(rawcsem.col), color = "Raw Score Metric"), size = 2) +
        scale_x_continuous(
          name = "Scale Score Metric (Lower Axis)",
          sec_axis = sec_axis(~ ., name = "Equivalent Raw Score (Upper Axis)",
                              breaks = scale_vals[seq(1, length(scale_vals), length.out = 8)],
                              labels = round(raw_vals[seq(1, length(raw_vals), length.out = 8)], 1))
        ) +
        scale_color_manual(values = setNames(c("#2c3e50", "#e74c3c"), c("Raw Score Metric", paste("Scale Metric (", method_label, ")")))) +
        labs(title = "Overlayed Conditional Standard Error of Measurement (CSEM) Functions",
             subtitle = "Upper horizontal axis indicates exact structural raw score alignments.",
             y = "Error Magnitude (CSEM)", color = "Measurement Metric") +
        theme_minimal() + theme(legend.position = "bottom")
      print(p)
    }
  }

  return(output_data)
}
