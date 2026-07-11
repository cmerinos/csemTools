#' Plot Conditional Standard Error of Measurement (CSEM) Curve
#'
#' @description
#' Creates a ggplot2 plot of the CSEM as a function of test scores. Optionally,
#' includes confidence intervals (ribbon around the CSEM curve) or the width
#' of the confidence band.
#'
#' @param score Numeric vector of scores (x-axis).
#' @param csem Numeric vector of CSEM values (y-axis). Must be same length as `score`.
#' @param lwr.ci Numeric vector of lower confidence limits for CSEM. Required if
#'   `plot.type = "CI"` or `"band"`. Must be same length as `score`.
#' @param upr.ci Numeric vector of upper confidence limits for CSEM. Required if
#'   `plot.type = "CI"` or `"band"`. Must be same length as `score`.
#' @param plot.type Character: `"CSEM"` (only the CSEM line/points),
#'   `"CI"` (CSEM with confidence ribbon), or `"band"` (width of confidence band).
#' @param title Character. Plot title. If NULL, a default title is generated.
#' @param xlab Character. X-axis label. If NULL, defaults to "Score".
#' @param ylab Character. Y-axis label. If NULL, auto-generated based on `plot.type`.
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
#' @examples
#' \donttest{
#' # Simulated data
#' scores <- 0:20
#' csem_vals <- 1 + 0.5 * abs(scores - 10) / 10
#' lwr <- csem_vals - 0.2
#' upr <- csem_vals + 0.2
#'
#' # Only CSEM curve
#' plotCSEM(score = scores, csem = csem_vals, plot.type = "CSEM")
#'
#' # With confidence interval ribbon
#' plotCSEM(score = scores, csem = csem_vals,
#'          lwr.ci = lwr, upr.ci = upr, plot.type = "CI")
#'
#' # Width of confidence band
#' plotCSEM(score = scores, csem = csem_vals,
#'          lwr.ci = lwr, upr.ci = upr, plot.type = "band")
#' }
#'
#' @export
plotCSEM <- function(score,
                     csem,
                     lwr.ci = NULL,
                     upr.ci = NULL,
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

  # --- Validate inputs ---
  if (!is.numeric(score) || !is.numeric(csem))
    stop("'score' and 'csem' must be numeric vectors.")
  if (length(score) != length(csem))
    stop("'score' and 'csem' must have the same length.")

  if (plot.type %in% c("CI", "band")) {
    if (is.null(lwr.ci) || is.null(upr.ci))
      stop("For plot.type = 'CI' or 'band', 'lwr.ci' and 'upr.ci' must be provided.")
    if (!is.numeric(lwr.ci) || !is.numeric(upr.ci))
      stop("'lwr.ci' and 'upr.ci' must be numeric vectors.")
    if (length(lwr.ci) != length(score) || length(upr.ci) != length(score))
      stop("'lwr.ci' and 'upr.ci' must have the same length as 'score'.")
  }

  # --- Prepare data frame for ggplot ---
  df <- data.frame(score = score, csem = csem)
  if (!is.null(lwr.ci)) df$lwr.ci <- lwr.ci
  if (!is.null(upr.ci)) df$upr.ci <- upr.ci

  if (plot.type == "band") {
    df$CI_band <- df$upr.ci - df$lwr.ci
  }

  # --- Default labels ---
  if (is.null(title)) {
    title <- switch(plot.type,
                    CSEM = "Conditional Standard Error of Measurement",
                    CI   = "CSEM with Confidence Intervals",
                    band = "Width of Confidence Bands")
  }
  if (is.null(xlab)) xlab <- "Score"
  if (is.null(ylab)) {
    ylab <- switch(plot.type,
                   CSEM = "CSEM",
                   CI   = "CSEM",
                   band = "CI Band Width")
  }

  # --- Build plot ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = score, y = csem)) +
    ggplot2::theme_minimal()

  if (plot.type == "CSEM") {
    p <- p + ggplot2::geom_line(color = color.line, linetype = line.type, size = 1) +
      ggplot2::geom_point(color = color.points, size = point.size)

  } else if (plot.type == "CI") {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr.ci, ymax = upr.ci),
                                  fill = color.band, alpha = 0.4) +
      ggplot2::geom_line(color = color.line, linetype = line.type, size = 1) +
      ggplot2::geom_point(color = color.points, size = point.size)

  } else if (plot.type == "band") {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = CI_band),
                                color = color.line, linetype = line.type, size = 1) +
      ggplot2::geom_point(ggplot2::aes(y = CI_band),
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
