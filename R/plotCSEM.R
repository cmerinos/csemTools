#' Plot Conditional Standard Error of Measurement (CSEM) Curve
#'
#' @description
#' Creates a ggplot2 plot of the CSEM as a function of test scores, or a plot
#' of the observed score with confidence interval band for the true score.
#'
#' @param score Numeric vector of scores (x-axis).
#' @param csem Numeric vector of CSEM values (y-axis for `"csem"` plot).
#'   For `"csemscore"`, this vector is used as the observed score (line diagonal).
#' @param lwr.ci Numeric vector of lower confidence limits (required for `"csemscore"`).
#' @param upr.ci Numeric vector of upper confidence limits (required for `"csemscore"`).
#' @param plot.type Character: `"csem"` (CSEM curve) or `"csemscore"` (confidence band).
#' @param cutoff Numeric vector of score values where vertical lines are added.
#'   Useful for highlighting cut scores or quantiles. Default = NULL.
#' @param title Character. Plot title. If NULL, a default title is generated.
#' @param xlab Character. X-axis label. If NULL, defaults to "Score".
#' @param ylab Character. Y-axis label. If NULL, auto-generated.
#' @param color.line Color for the line (default = "black").
#' @param color.points Color for points (default = "darkred").
#' @param color.band Color for confidence band/ribbon (default = "lightblue").
#' @param color.cutoff Color for cutoff lines (default = "gray50").
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
#' # Example with simulated data
#' scores <- 0:20
#' csem_vals <- 1 + 0.5 * abs(scores - 10) / 10
#'
#' # CSEM curve
#' plotCSEM(score = scores, csem = csem_vals, plot.type = "csem")
#'
#' # Confidence band for true score (needs lwr.ci and upr.ci)
#' lwr <- scores - 1.96 * csem_vals
#' upr <- scores + 1.96 * csem_vals
#' plotCSEM(score = scores, csem = scores, lwr.ci = lwr, upr.ci = upr,
#'          plot.type = "csemscore", cutoff = c(5, 10, 15))
#' }
#'
#' @export
plotCSEM <- function(score,
                     csem,
                     lwr.ci = NULL,
                     upr.ci = NULL,
                     plot.type = c("csem", "csemscore"),
                     cutoff = NULL,
                     title = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     color.line = "black",
                     color.points = "darkred",
                     color.band = "lightblue",
                     color.cutoff = "gray50",
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

  if (plot.type == "csemscore") {
    if (is.null(lwr.ci) || is.null(upr.ci))
      stop("For plot.type = 'csemscore', 'lwr.ci' and 'upr.ci' must be provided.")
    if (!is.numeric(lwr.ci) || !is.numeric(upr.ci))
      stop("'lwr.ci' and 'upr.ci' must be numeric vectors.")
    if (length(lwr.ci) != length(score) || length(upr.ci) != length(score))
      stop("'lwr.ci' and 'upr.ci' must have the same length as 'score'.")
  }

  if (!is.null(cutoff) && !is.numeric(cutoff))
    stop("'cutoff' must be a numeric vector.")

  # --- Prepare data frame ---
  df <- data.frame(score = score, csem = csem)
  if (!is.null(lwr.ci)) df$lwr.ci <- lwr.ci
  if (!is.null(upr.ci)) df$upr.ci <- upr.ci

  # --- Default labels ---
  if (is.null(title)) {
    title <- switch(plot.type,
                    csem = "Conditional Standard Error of Measurement",
                    csemscore = "Observed Score with Confidence Interval")
  }
  if (is.null(xlab)) xlab <- "Score"
  if (is.null(ylab)) {
    ylab <- switch(plot.type,
                   csem = "CSEM",
                   csemscore = "Score")
  }

  # --- Build plot with theme_classic ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = score)) +
    ggplot2::theme_classic()

  if (plot.type == "csem") {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = csem),
                                color = color.line, linetype = line.type, size = 1) +
      ggplot2::geom_point(ggplot2::aes(y = csem),
                          color = color.points, size = point.size)

  } else if (plot.type == "csemscore") {
    # Band for true score
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr.ci, ymax = upr.ci),
                                  fill = color.band, alpha = 0.4) +
      # Diagonal line (observed score = true score)
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                           color = "gray50", size = 0.5) +
      # Optionally, add observed score points (if desired)
      ggplot2::geom_point(ggplot2::aes(y = score),
                          color = color.points, size = point.size, alpha = 0.5)
  }

  # --- Add cutoff lines if provided ---
  if (!is.null(cutoff)) {
    p <- p + ggplot2::geom_vline(xintercept = cutoff,
                                 color = color.cutoff,
                                 linetype = "dashed",
                                 size = 0.5)
  }

  p <- p + ggplot2::labs(title = title, x = xlab, y = ylab)

  # --- Save if requested ---
  if (!is.null(save.path)) {
    ggplot2::ggsave(filename = save.path, plot = p, width = width, height = height, dpi = 300)
    message("Plot saved to: ", save.path)
  }

  print(p)
  invisible(p)
}
