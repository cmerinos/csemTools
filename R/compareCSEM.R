#' Compare Conditional SEM with Global SEM
#'
#' @description
#' Compares a conditional standard error of measurement (CSEM) function against the
#' global (constant) SEM. The global SEM is computed from the observed score standard
#' deviation and a reliability coefficient (e.g., alpha or omega). Optionally,
#' confidence intervals for the global SEM can be obtained from user supplied
#' reliability confidence bounds. The function produces a data frame with ratios
#' and, if requested, plots of the CSEM vs. global SEM and/or the ratio.
#'
#' @param data A data frame containing at least the columns specified in
#'   \code{raw.score} and \code{CSEM}. Optionally, \code{lwr.ci.CSEM} and \code{upr.ci.CSEM}
#'   can be provided for confidence bands. Ignored if \code{score} and \code{csem} are provided.
#' @param raw.score Character. Name of the column with observed scores (default = "raw.score").
#' @param CSEM Character. Name of the column with conditional SEM values (default = "CSEM").
#' @param lwr.ci.CSEM Character. Optional name of the column with lower confidence limits for CSEM.
#' @param upr.ci.CSEM Character. Optional name of the column with upper confidence limits for CSEM.
#' @param score Numeric vector of observed scores. If provided, overrides \code{data}.
#' @param csem Numeric vector of conditional SEM values. If provided, overrides \code{data}.
#' @param lwr.ci Numeric vector of lower confidence limits for CSEM (optional).
#' @param upr.ci Numeric vector of upper confidence limits for CSEM (optional).
#' @param cutoff Numeric vector of score values where vertical lines are added to the plots.
#'   Useful for highlighting cut scores or quantiles. Default = NULL.
#' @param sd.score Numeric. Standard deviation of the observed scores (calculated externally from \code{raw.score}).
#' @param reliability Numeric. Reliability coefficient (e.g., alpha, omega, Gilmer-Feldt) used to compute global SEM.
#' @param reliab.lwrci Numeric. Optional lower confidence bound for \code{reliability}. If provided, must be used with \code{reliab.uprci}.
#' @param reliab.uprci Numeric. Optional upper confidence bound for \code{reliability}.
#' @param plot Character. Type of plot to produce: \code{"csem"} (CSEM vs. global SEM),
#'   \code{"ratio"} (ratio CSEM/global SEM), \code{"all"} (both plots), or \code{"none"} (no plots). Default = \code{"none"}.
#' @param digits Numeric. Number of decimal places for rounding output numeric columns (default = 3).
#' @param conf.level Numeric. Confidence level for intervals (default = 0.95). Used only for labelling.
#' @param ... Additional arguments passed to \code{ggplot2::geom_line()} and \code{ggplot2::geom_ribbon()}.
#'
#' @details
#' The method for comparing conditional SEM against global SEM is inspired by the
#' reliability representativeness approach described in McNeish and Dumas (2025),
#' adapted here for standard errors of measurement rather than reliability.
#'
#' The global SEM is defined as:
#' \deqn{SEM_{global} = SD_{observed} \times \sqrt{1 - reliability}}
#'
#' If reliability confidence bounds (\code{reliab.lwrci}, \code{reliab.uprci}) are supplied,
#' they are transformed to SEM bounds using the inverse relationship:
#' \deqn{SEM\_lower = SD_{observed} \times \sqrt{1 - reliability\_upper}}
#' \deqn{SEM\_upper = SD_{observed} \times \sqrt{1 - reliability\_lower}}
#' because lower reliability implies larger measurement error. These bounds are used to
#' draw a confidence band around the global SEM line.
#'
#' The ratio is computed as \code{CSEM / SEM_global}. Values above 1 indicate that the
#' conditional SEM is larger (worse precision) than the global average; values below 1
#' indicate better precision.
#'
#' @return
#' \itemize{
#'   \item \code{data}: A data frame with columns \code{raw.score}, \code{CSEM}, optionally
#'     \code{lwr.ci.CSEM}/\code{upr.ci.CSEM}, and \code{ratio} (rounded to \code{digits}).
#'   \item \code{global_sem}: A data frame with two columns: \code{value} (label) and \code{estimate}
#'     (numeric). Contains \code{sem.global} and, if provided, \code{lwr.ci} and \code{upr.ci}.
#'   \item \code{plot}: A ggplot2 object (or list of two objects if \code{plot = "all"}); \code{NULL} if \code{plot = "none"}.
#' }
#'
#' @references
#' McNeish, D., Dumas, D. (2025). Reliability representativeness: How well does coefficient alpha
#' summarize reliability across the score distribution? \emph{Behavior Research Methods}, 57, 93.
#' \doi{10.3758/s13428-025-02611-8}
#'
#' @importFrom stats qt
#' @importFrom rlang .data
#'
#' @examples
#' \donttest{
#' # Example using data frame
#' df <- data.frame(
#'   raw.score = 10:40,
#'   CSEM = 2.3 + 0.005 * (10:40 - 25)^2,
#'   lwr.ci.CSEM = 2.1 + 0.005 * (10:40 - 25)^2,
#'   upr.ci.CSEM = 2.5 + 0.005 * (10:40 - 25)^2
#' )
#' result <- compareCSEM(data = df,
#'                       sd.score = 6.0,
#'                       reliability = 0.85,
#'                       reliab.lwrci = 0.82,
#'                       reliab.uprci = 0.88,
#'                       plot = "all",
#'                       cutoff = c(15, 30))
#'
#' # Example using vectors directly
#' scores <- 10:40
#' csem_vals <- 2.3 + 0.005 * (scores - 25)^2
#' compareCSEM(data = NULL, score = scores, csem = csem_vals,
#'             sd.score = 6.0, reliability = 0.85,
#'             plot = "csem")
#' }
#'
#' @export
compareCSEM <- function(data,
                        raw.score = "raw.score",
                        CSEM = "CSEM",
                        lwr.ci.CSEM = NULL,
                        upr.ci.CSEM = NULL,
                        score = NULL,
                        csem = NULL,
                        lwr.ci = NULL,
                        upr.ci = NULL,
                        cutoff = NULL,
                        sd.score,
                        reliability,
                        reliab.lwrci = NULL,
                        reliab.uprci = NULL,
                        plot = c("none", "all", "csem", "ratio"),
                        digits = 3,
                        conf.level = 0.95,
                        ...) {

  # --- Match plot argument ---
  plot <- match.arg(plot)

  # --- Determine input source: vectors or data frame --------------------
  # If both score and csem are provided, they take precedence over data.
  if (!is.null(score) && !is.null(csem)) {
    if (!is.data.frame(data) && !missing(data)) {
      message("Both 'score'/'csem' and 'data' provided. Using vectors 'score' and 'csem'; 'data' ignored.")
    }
    # Validate vectors
    if (!is.numeric(score) || !is.numeric(csem))
      stop("'score' and 'csem' must be numeric vectors.")
    if (length(score) != length(csem))
      stop("'score' and 'csem' must have the same length.")
    if (!is.null(lwr.ci) && length(lwr.ci) != length(score))
      stop("'lwr.ci' (if provided) must have the same length as 'score'.")
    if (!is.null(upr.ci) && length(upr.ci) != length(score))
      stop("'upr.ci' (if provided) must have the same length as 'score'.")

    scores <- score
    csem_vals <- csem
    csem_lwr <- lwr.ci
    csem_upr <- upr.ci
    has_cisem <- !is.null(csem_lwr) && !is.null(csem_upr)

  } else {
    # Use data frame
    if (missing(data) || !is.data.frame(data))
      stop("If 'score' and 'csem' are not provided, 'data' must be a data frame.")
    required_cols <- c(raw.score, CSEM)
    missing_cols <- required_cols[!required_cols %in% names(data)]
    if (length(missing_cols) > 0) {
      stop("Missing required columns in 'data': ", paste(missing_cols, collapse = ", "))
    }
    scores <- data[[raw.score]]
    csem_vals <- data[[CSEM]]

    # Optional confidence interval columns from data
    has_cisem <- !is.null(lwr.ci.CSEM) && !is.null(upr.ci.CSEM)
    if (has_cisem) {
      if (!lwr.ci.CSEM %in% names(data) || !upr.ci.CSEM %in% names(data)) {
        stop("Columns specified in 'lwr.ci.CSEM' and/or 'upr.ci.CSEM' not found in data.")
      }
      csem_lwr <- data[[lwr.ci.CSEM]]
      csem_upr <- data[[upr.ci.CSEM]]
    } else {
      csem_lwr <- NULL
      csem_upr <- NULL
    }
  }

  # --- Validate numeric inputs -----------------------------------------
  if (!is.numeric(scores) || !is.numeric(csem_vals))
    stop("'score' and 'csem' (or columns from 'data') must be numeric vectors.")

  if (!is.numeric(sd.score) || length(sd.score) != 1 || sd.score <= 0)
    stop("'sd.score' must be a single positive number.")

  if (!is.numeric(reliability) || length(reliability) != 1 || reliability <= 0 || reliability >= 1)
    stop("'reliability' must be a single number between 0 and 1 (exclusive).")

  if (!is.null(reliab.lwrci) && !is.null(reliab.uprci)) {
    if (!is.numeric(reliab.lwrci) || !is.numeric(reliab.uprci) ||
        reliab.lwrci >= reliab.uprci || reliab.lwrci <= 0 || reliab.uprci >= 1)
      stop("'reliab.lwrci' and 'reliab.uprci' must be numbers with 0 < lwr < upr < 1.")
  } else if (xor(is.null(reliab.lwrci), is.null(reliab.uprci))) {
    stop("Both 'reliab.lwrci' and 'reliab.uprci' must be provided together, or both NULL.")
  }

  if (has_cisem && (is.null(csem_lwr) || is.null(csem_upr)))
    stop("If confidence intervals for CSEM are requested, both lower and upper bounds must be provided.")

  if (!is.null(cutoff) && !is.numeric(cutoff))
    stop("'cutoff' must be a numeric vector.")

  if (!is.numeric(digits) || length(digits) != 1 || digits < 0)
    stop("'digits' must be a non-negative integer.")

  # --- Compute global SEM and confidence interval ----------------------
  sem_global <- sd.score * sqrt(1 - reliability)

  if (!is.null(reliab.lwrci)) {
    sem_global_lwr <- sd.score * sqrt(1 - reliab.uprci)
    sem_global_upr <- sd.score * sqrt(1 - reliab.lwrci)
    has_ci_global <- TRUE
  } else {
    sem_global_lwr <- sem_global_upr <- NA_real_
    has_ci_global <- FALSE
  }

  # --- Build output data frame -----------------------------------------
  out_df <- data.frame(
    raw.score = scores,
    CSEM = csem_vals
  )
  if (has_cisem) {
    out_df$lwr.ci.CSEM <- csem_lwr
    out_df$upr.ci.CSEM <- csem_upr
  }
  out_df$ratio <- csem_vals / sem_global

  # Round numeric columns
  numeric_cols <- sapply(out_df, is.numeric)
  out_df[numeric_cols] <- lapply(out_df[numeric_cols], round, digits = digits)

  # --- Global SEM info as data frame -----------------------------------
  global_sem_df <- data.frame(
    value = "sem.global",
    estimate = round(sem_global, digits)
  )
  if (has_ci_global) {
    global_sem_df <- rbind(
      global_sem_df,
      data.frame(value = "lwr.ci", estimate = round(sem_global_lwr, digits)),
      data.frame(value = "upr.ci", estimate = round(sem_global_upr, digits))
    )
  }

  # --- Prepare plots --------------------------------------------------
  plot_list <- list()
  has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

  if (plot != "none" && !has_ggplot2) {
    warning("Package 'ggplot2' is not installed. Plots will be omitted.")
    plot <- "none"
  }

  if (plot != "none" && has_ggplot2) {
    # Use original (unrounded) values for plotting
    plot_df <- data.frame(
      raw.score = scores,
      CSEM = csem_vals,
      ratio = csem_vals / sem_global
    )
    if (has_cisem) {
      plot_df$lwr.ci.CSEM <- csem_lwr
      plot_df$upr.ci.CSEM <- csem_upr
    }

    # ---- Plot 1: CSEM vs Global SEM ----
    if (plot %in% c("csem", "all")) {
      p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = raw.score))
      if (has_cisem) {
        p1 <- p1 + ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr.ci.CSEM,
                                                     ymax = upr.ci.CSEM),
                                        fill = "grey70", alpha = 0.5)
      }
      p1 <- p1 +
        ggplot2::geom_line(ggplot2::aes(y = CSEM, colour = "Conditional SEM"), size = 1.2) +
        ggplot2::geom_hline(yintercept = sem_global, linetype = "dashed", colour = "red", size = 1) +
        ggplot2::labs(x = "Observed Score", y = "Standard Error of Measurement",
                      title = "Conditional vs. Global SEM",
                      colour = "Line") +
        ggplot2::theme_classic()

      # Add global SEM confidence band (rectangular band)
      if (has_ci_global) {
        p1 <- p1 + ggplot2::annotate("rect",
                                     xmin = min(scores), xmax = max(scores),
                                     ymin = sem_global_lwr, ymax = sem_global_upr,
                                     fill = "red", alpha = 0.2)
      }

      # Add cutoff vertical lines if provided
      if (!is.null(cutoff)) {
        p1 <- p1 + ggplot2::geom_vline(xintercept = cutoff,
                                       linetype = "dashed", color = "gray50", size = 0.5)
      }

      plot_list$csem_plot <- p1
    }

    # ---- Plot 2: Ratio CSEM / Global SEM ----
    if (plot %in% c("ratio", "all")) {
      p2 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = raw.score, y = .data$ratio)) +
        ggplot2::geom_line(size = 1.2, colour = "steelblue") +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "darkred", size = 0.8) +
        ggplot2::labs(x = "Observed Score", y = "Ratio (CSEM / Global SEM)",
                      title = "Relative Precision: CSEM vs. Global SEM") +
        ggplot2::theme_classic()

      # Add ratio confidence band if global SEM CI exists
      if (has_ci_global) {
        ratio_lwr <- csem_vals / sem_global_upr
        ratio_upr <- csem_vals / sem_global_lwr
        p2 <- p2 + ggplot2::geom_ribbon(ggplot2::aes(ymin = ratio_lwr, ymax = ratio_upr),
                                        fill = "steelblue", alpha = 0.2)
      }

      # Add cutoff vertical lines if provided
      if (!is.null(cutoff)) {
        p2 <- p2 + ggplot2::geom_vline(xintercept = cutoff,
                                       linetype = "dashed", color = "gray50", size = 0.5)
      }

      plot_list$ratio_plot <- p2
    }
  }

  # --- Output ---------------------------------------------------------
  result <- list(data = out_df,
                 global_sem = global_sem_df,
                 plot = NULL)
  if (plot == "all" && has_ggplot2) {
    result$plot <- plot_list
  } else if (plot %in% c("csem", "ratio") && has_ggplot2) {
    result$plot <- if (plot == "csem") plot_list$csem_plot else plot_list$ratio_plot
  }
  return(result)
}
