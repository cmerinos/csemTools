#' Mollenkopf-Feldt method for Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Implements the Mollenkopf-Feldt method, which applies polynomial regression
#' to the individual error variance estimates from Feldt & Qualls (1996).
#' This provides smoothed CSEM values along the score scale.
#'
#' @inheritParams csemFeldtQualls
#' @param degree Integer. Degree of the polynomial (default 2).
#' @param full.range Logical. If `TRUE`, evaluates the smoothed CSEM for every
#'   integer score from `score.range[1]` to `score.range[2]` (requires `score.range`).
#' @param ci Logical. If `TRUE`, compute confidence intervals for the true score.
#' @param conf.level Numeric. Confidence level (default 0.95).
#' @param digits Integer. Rounding for output.
#' @param score.range Numeric vector length 2. Required if `full.range = TRUE`.
#' @param na.rm Logical. Remove rows with missing values.
#'
#' @return A list with elements:
#' \item{CSEM}{data.frame with columns `score`, `n`, `CSEM.smooth`, and
#'   if `ci = TRUE`, `lwr.ci` and `upr.ci`.}
#' \item{...}{No `binned.CSEM` (can be added later).}
#'
#' @references
#' Mollenkopf, W. G. (1949). Variation of the standard error of measurement.
#'   *Psychometrika*, 14(3), 189–229.
#' Feldt, L. S., & Qualls, A. L. (1996). Estimation of measurement error variance
#'   at specific score levels. *Journal of Educational Measurement*, 33(2), 141–156.
#'
#' @export
csemMF <- function(data,
                   n.parts = NULL,
                   part_items = NULL,
                   min.items.per.part = 2,
                   degree = 2,
                   full.range = FALSE,
                   ci = FALSE,
                   conf.level = 0.95,
                   digits = 3,
                   score.range = NULL,
                   na.rm = TRUE) {

  # --- Validaciones iniciales ---
  if (!is.data.frame(data) && !is.matrix(data))
    stop("`data` must be a data frame or matrix.")
  data <- as.data.frame(data)
  if (na.rm) data <- stats::na.omit(data)
  if (anyNA(data)) stop("Missing values present. Set na.rm = TRUE to remove them.")

  J <- ncol(data)
  n_persons <- nrow(data)
  if (n_persons < 2) stop("At least 2 persons required.")
  if (J < 2) stop("At least 2 items required.")

  # --- Construcción de partes (cálculo de Y_i) ---
  if (!is.null(part_items)) {
    if (!is.list(part_items)) stop("`part_items` must be a list.")
    n.parts <- length(part_items)
    part_scores <- matrix(NA, nrow = n_persons, ncol = n.parts)
    items_per_part <- numeric(n.parts)
    for (j in 1:n.parts) {
      cols <- part_items[[j]]
      if (is.character(cols)) {
        if (!all(cols %in% colnames(data)))
          stop("Some column names in part_items[[", j, "]] not found.")
        cols <- which(colnames(data) %in% cols)
      } else if (is.numeric(cols)) {
        if (any(cols < 1 | cols > J))
          stop("Column indices out of range in part_items[[", j, "]].")
      } else {
        stop("part_items[[", j, "]] must be character or integer vector.")
      }
      part_scores[, j] <- rowSums(data[, cols, drop = FALSE], na.rm = TRUE)
      items_per_part[j] <- length(cols)
    }
    if (any(items_per_part < min.items.per.part))
      warning("Some parts have fewer than ", min.items.per.part, " items.")
    if (length(unique(items_per_part)) > 1)
      warning("Parts have unequal length; results may be biased.")
    k <- J / n.parts
    d <- J / k   # = n.parts
  } else {
    if (is.null(n.parts)) n.parts <- J
    if (n.parts < 2) stop("n.parts must be at least 2.")
    if (n.parts > J) stop("n.parts cannot exceed number of items.")
    idx_split <- split(1:J, cut(1:J, breaks = n.parts, labels = FALSE))
    part_scores <- matrix(NA, nrow = n_persons, ncol = n.parts)
    items_per_part <- numeric(n.parts)
    for (j in 1:n.parts) {
      cols <- idx_split[[j]]
      part_scores[, j] <- rowSums(data[, cols, drop = FALSE], na.rm = TRUE)
      items_per_part[j] <- length(cols)
    }
    if (any(items_per_part < min.items.per.part))
      warning("Some parts have fewer than ", min.items.per.part, " items.")
    k <- J / n.parts
    d <- n.parts
  }

  # --- Calcular Y_i (estimación individual de varianza de error) ---
  Xij <- part_scores
  barX_i <- rowMeans(Xij, na.rm = TRUE)
  barX_j <- colMeans(Xij, na.rm = TRUE)
  M <- mean(barX_j)
  dev <- sweep(Xij, 1, barX_i, "-")
  dev <- sweep(dev, 2, barX_j - M, "-")
  SS_i <- rowSums(dev^2, na.rm = TRUE)
  var_adj_i <- SS_i / (n.parts - 1)
  Y_i <- d * var_adj_i

  total <- rowSums(data)   # puntaje total original

  # --- Rango para truncar CI (si no se da score.range, usar rango observado) ---
  if (!is.null(score.range)) {
    if (!is.numeric(score.range) || length(score.range) != 2)
      stop("score.range must be a numeric vector of length 2.")
    score_min_teo <- score.range[1]
    score_max_teo <- score.range[2]
  } else {
    score_min_teo <- min(total, na.rm = TRUE)
    score_max_teo <- max(total, na.rm = TRUE)
  }

  # --- Regresión polinómica de Y_i sobre total ---
  df_persons <- data.frame(total = total, Y = Y_i)
  # Eliminar posibles NA (no debería haber, pero por seguridad)
  df_persons <- df_persons[complete.cases(df_persons), ]
  if (nrow(df_persons) < degree + 1)
    stop("Not enough persons to fit polynomial of degree ", degree)

  fit <- lm(Y ~ poly(total, degree, raw = TRUE), data = df_persons)

  # --- Puntajes para predicción ---
  if (full.range) {
    if (is.null(score.range))
      stop("full.range = TRUE requires 'score.range'.")
    pred_scores <- seq(score.range[1], score.range[2], by = 1)
  } else {
    pred_scores <- sort(unique(total))
  }

  # Predicción de varianza y CSEM
  pred_var <- predict(fit, newdata = data.frame(total = pred_scores))
  pred_var <- pmax(pred_var, 0)
  csem_smooth <- sqrt(pred_var)

  # Frecuencias reales (número de personas con cada puntaje)
  n_vals <- sapply(pred_scores, function(s) sum(total == s))

  result_df <- data.frame(score = pred_scores, n = n_vals,
                          CSEM.smooth = round(csem_smooth, digits),
                          stringsAsFactors = FALSE)

  # --- Intervalos de confianza ---
  if (ci) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    lwr <- result_df$score - z * result_df$CSEM.smooth
    upr <- result_df$score + z * result_df$CSEM.smooth
    lwr <- pmax(lwr, score_min_teo)
    upr <- pmin(upr, score_max_teo)
    result_df$lwr.ci <- round(lwr, digits)
    result_df$upr.ci <- round(upr, digits)
  }

  # --- Salida ---
  out <- list(CSEM = result_df)
  # (Opcional: se podría añadir $binned.CSEM más adelante)
  return(out)
}
