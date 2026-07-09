#' Plot Conditional Standard Error of Measurement (CSEM) Curve
#'
#' @description
#' Creates a ggplot2 plot of the CSEM as a function of test scores, optionally
#' including confidence intervals or confidence band widths. Can be used with
#' any CSEM function that returns a data frame with columns `Score`, `CSEM`, and
#' optionally `lwr.ci`, `upr.ci`.
#'
#' @param data A data frame containing at least columns `Score` and `CSEM`.
#'   Optionally `lwr.ci` and `upr.ci` for confidence intervals.
#' @param x Character. Name of the column to use on x-axis (default = "Score").
#' @param y Character. Name of the column to use on y-axis for CSEM (default = "CSEM").
#' @param plot.type Character: `"CSEM"` (only the CSEM line/points),
#'   `"CI"` (CSEM with error bars/ribbon), or `"band"` (width of confidence band).
#' @param title Character. Plot title. If NULL, a default title is generated.
#' @param xlab Character. X-axis label. If NULL, uses `x`.
#' @param ylab Character. Y-axis label. If NULL, auto-generated.
#' @param color.line Color for the line (default = "black").
#' @param color.points Color for points (default = "darkred").
#' @param color.band Color for confidence band/ribbon (default = "lightblue").
#' @param line.type Line type (default = "solid").
#' @param point.size Size of points (default = 2).
#' @param save.path Optional file path to save plot as PNG. If NULL, plot is not saved.
#' @param width Numeric. Width of saved plot in inches (default = 8).
#' @param height Numeric. Height of saved plot in inches (default = 6).
#'
#' @return A ggplot2 object (invisibly). The plot is also printed.
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \donttest{
#' ## Load data
#' library(EFA.dimensions)
#' data("data_RSE")
#'
#' ## Recode negative items
#' data_RSE[c("Q3", "Q5", "Q8", "Q9", "Q10")] <- 5 - data_RSE[c("Q3", "Q5", "Q8", "Q9", "Q10")]
#'
#' ## Split in two halves
#' RSE.namesHalf <- checkSplit(data = data_RSE, method = "difficulty")
#'
#' ## Items in the halves?
#' RSE.namesHalf$half1
#' RSE.namesHalf$half2
#'
#' # Mollenkopst-Feldt method
#' mfres <- csemMF(RSE.namesHalf$half1,
#' RSE.namesHalf$half1,
#' degree = 2, ci = TRUE)
#'
#' plotCSEM(data = mfres$CSEM,
#' x = "score",
#' y = "CSEM.smooth", plot.type = "CSEM")
#' }
plotCSEM <- function(data,
                     x = "Score",
                     y = "CSEM",
                     plot.type = c("CSEM", "CI", "band"),
                     title = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     color.line = "black",
                     color.points = "darkred",
                     color.band = "lightblue",
                     line.type = "solid",
                     point.size = 2,
                     save.path = NULL,
                     width = 8,
                     height = 6) {

  plot.type <- match.arg(plot.type)

  # --- Check required columns ---
  if (!x %in% names(data)) stop(paste("Column", x, "not found in data."))
  if (!y %in% names(data)) stop(paste("Column", y, "not found in data."))
  if (plot.type == "CI" && (!"lwr.ci" %in% names(data) || !"upr.ci" %in% names(data))) {
    stop("For plot.type = 'CI', data must have columns 'lwr.ci' and 'upr.ci'.")
  }
  if (plot.type == "band" && !"CI_band" %in% names(data)) {
    # For band we need CI_band (width). If not present, we can compute as upr.ci - lwr.ci
    if ("lwr.ci" %in% names(data) && "upr.ci" %in% names(data)) {
      data$CI_band <- data$upr.ci - data$lwr.ci
    } else {
      stop("For plot.type = 'band', data must have columns 'lwr.ci' and 'upr.ci', or a 'CI_band' column.")
    }
  }

  # --- Default labels ---
  if (is.null(title)) {
    title <- switch(plot.type,
                    CSEM = "Conditional Standard Error of Measurement",
                    CI   = "CSEM with Confidence Intervals",
                    band = "Width of Confidence Bands")
  }
  if (is.null(xlab)) xlab <- x
  if (is.null(ylab)) {
    ylab <- switch(plot.type,
                   CSEM = "CSEM",
                   CI   = "Score",
                   band = "CI Band Width")
  }

  # --- Build plot ---
  p <- ggplot2::ggplot(data, ggplot2::aes_string(x = x, y = y)) +
    ggplot2::theme_minimal()

  if (plot.type == "CSEM") {
    p <- p + ggplot2::geom_line(color = color.line, linetype = line.type, size = 1) +
      ggplot2::geom_point(color = color.points, size = point.size)

  } else if (plot.type == "CI") {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "lwr.ci", ymax = "upr.ci"),
                                  fill = color.band, alpha = 0.4) +
      ggplot2::geom_line(color = color.line, linetype = line.type, size = 1) +
      ggplot2::geom_point(color = color.points, size = point.size)

  } else if (plot.type == "band") {
    p <- p + ggplot2::geom_line(ggplot2::aes_string(y = "CI_band"),
                                color = color.line, linetype = line.type, size = 1) +
      ggplot2::geom_point(ggplot2::aes_string(y = "CI_band"),
                          color = color.points, size = point.size) +
      ggplot2::ylab(ylab)
  }

  p <- p + ggplot2::labs(title = title, x = xlab, y = ylab)

  # --- Save if requested ---
  if (!is.null(save.path)) {
    ggplot2::ggsave(filename = save.path, plot = p, width = width, height = height, dpi = 300)
    message("Plot saved to: ", save.path)
  }

  # --- Print and return invisibly ---
  print(p)
  invisible(p)
}
