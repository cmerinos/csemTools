#' Plot Conditional Standard Error of Measurement (CSEM) Curve
#'
#' @description
#' Creates a ggplot2 plot of the CSEM as a function of test scores, optionally
#' including confidence intervals (bands) for the true score or the width of the
#' confidence band. The function is designed to work with any CSEM function that
#' returns a data frame containing at least columns for score and CSEM, and
#' optionally `lwr.ci` and `upr.ci`.
#'
#' Three plot types are available:
#' \itemize{
#'   \item \code{"CSEM"}: plots the CSEM curve (line + points).
#'   \item \code{"CI"}: plots the observed score with a confidence band around it
#'         (based on `lwr.ci` and `upr.ci`) and optionally overlays the CSEM curve
#'         on a secondary axis (scaled by `csem.scale`).
#'   \item \code{"band"}: plots the width of the confidence band (i.e., `upr.ci - lwr.ci`).
#' }
#'
#' @param data A data frame containing at least the columns specified by `x` and `y`.
#'   For `plot.type = "CI"`, the data must also contain columns `lwr.ci` and `upr.ci`.
#' @param x Character. Name of the column to use on the x-axis (default = `"Score"`).
#' @param y Character. Name of the column to use on the y-axis for the CSEM values
#'   (default = `"CSEM"`).
#' @param plot.type Character. One of `"CSEM"`, `"CI"`, or `"band"`. Determines the
#'   type of plot to generate. Default = `"CSEM"`.
#' @param title Character. Plot title. If `NULL`, a default title is generated
#'   based on `plot.type`.
#' @param xlab Character. Label for the x-axis. If `NULL`, the value of `x` is used.
#' @param ylab Character. Label for the y-axis. If `NULL`, a default label is
#'   generated based on `plot.type`.
#' @param color.line Color for the main line (default = `"black"`).
#' @param color.points Color for the points (default = `"darkred"`).
#' @param color.band Color for the confidence band/ribbon (default = `"lightblue"`).
#' @param line.type Line type for the main line (default = `"solid"`).
#' @param point.size Size of points (default = `2`).
#' @param csem.scale Numeric scaling factor for the CSEM curve when `plot.type = "CI"`.
#'   The CSEM values are multiplied by this factor to make them visible on the
#'   secondary axis. Set to `0` to suppress the CSEM curve entirely. Default = `1`.
#' @param save.path Optional file path (e.g., `"plot.png"`) to save the plot as a
#'   PNG image. If `NULL`, the plot is not saved. Default = `NULL`.
#' @param width Numeric. Width of the saved plot in inches (default = `8`).
#' @param height Numeric. Height of the saved plot in inches (default = `6`).
#'
#' @details
#' The function uses `ggplot2` and returns a ggplot object invisibly. The plot
#' is also printed automatically.
#'
#' When `plot.type = "CI"`, the confidence band is drawn around the observed score
#' (identity line) using `lwr.ci` and `upr.ci`. The CSEM curve is optionally
#' overlaid on a secondary y‑axis, scaled by `csem.scale`. This allows the user to
#' visualize both the precision of the observed score and the conditional error
#' in a single plot.
#'
#' @return A `ggplot` object (invisibly). The plot is also printed.
#'
#' @importFrom ggplot2 ggplot aes_string geom_line geom_point geom_ribbon geom_abline scale_y_continuous sec_axis theme_minimal labs ggsave
#'
#' @examples
#' \donttest{
#' # Example using csemThorndike output
#' # (Assuming half1, half2 and csemThorndike output exist)
#' # plotCSEM(data = res$CSEM, x = "score", y = "CSEM.smooth")
#' }
#'
#' @export
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
                     csem.scale = 1,
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
                    CI   = "Observed Score with Confidence Interval",
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
  p <- ggplot2::ggplot(data, ggplot2::aes_string(x = x)) +
    ggplot2::theme_minimal()

  if (plot.type == "CSEM") {
    p <- p + ggplot2::geom_line(ggplot2::aes_string(y = y), color = color.line, linetype = line.type, size = 1) +
      ggplot2::geom_point(ggplot2::aes_string(y = y), color = color.points, size = point.size)

  } else if (plot.type == "CI") {
    # --- CI para el puntaje verdadero: banda alrededor de la diagonal ---
    p <- p + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "lwr.ci", ymax = "upr.ci"),
                                  fill = color.band, alpha = 0.4) +
      # Línea de identidad
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", size = 0.5) +
      # Puntos del puntaje observado
      ggplot2::geom_point(ggplot2::aes_string(y = x), color = "gray30", size = 1, alpha = 0.5)

    # Añadir CSEM en eje secundario si csem.scale != 0
    if (csem.scale != 0) {
      data$CSEM_scaled <- data[[y]] * csem.scale
      p <- p + ggplot2::geom_line(ggplot2::aes_string(y = "CSEM_scaled"), color = color.line, linetype = line.type, size = 1) +
        ggplot2::geom_point(ggplot2::aes_string(y = "CSEM_scaled"), color = color.points, size = point.size) +
        ggplot2::scale_y_continuous(
          name = "Score",
          sec.axis = ggplot2::sec_axis(~ . / csem.scale, name = "CSEM")
        )
    } else {
      p <- p + ggplot2::scale_y_continuous(name = "Score")
    }
    p <- p + ggplot2::labs(title = title, x = xlab)

  } else if (plot.type == "band") {
    p <- p + ggplot2::geom_line(ggplot2::aes_string(y = "CI_band"),
                                color = color.line, linetype = line.type, size = 1) +
      ggplot2::geom_point(ggplot2::aes_string(y = "CI_band"),
                          color = color.points, size = point.size) +
      ggplot2::ylab(ylab) +
      ggplot2::labs(title = title, x = xlab)
  }

  # --- Add labels (for non-CI types, already added) ---
  if (plot.type != "CI") {
    p <- p + ggplot2::labs(y = ylab)
  }

  # --- Save if requested ---
  if (!is.null(save.path)) {
    ggplot2::ggsave(filename = save.path, plot = p, width = width, height = height, dpi = 300)
    message("Plot saved to: ", save.path)
  }

  # --- Print and return invisibly ---
  print(p)
  invisible(p)
}
