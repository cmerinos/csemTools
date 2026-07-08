#' Mollenkopf-Feldt method for Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Implements the Mollenkopf-Feldt method, which applies polynomial regression
#' to the individual error variance estimates from Feldt & Qualls (1996).
#' This provides smoothed CSEM values along the score scale.
#'
#' @param data A data frame or matrix with item responses (subjects in rows,
#'   items in columns).
#' @param n.parts Integer. Number of parts into which the test will be split
#'   (by column order, as balanced as possible). Ignored if `part_items` is provided.
#'   Default is `NULL`, which sets `n.parts = ncol(data)` (each item as a part).
#' @param part_items Optional list. Each element is a character vector of column names
#'   or an integer vector of column indices defining the items in that part.
#'   If provided, `n.parts` is ignored.
#' @param min.items.per.part Integer. Minimum number of items per part (default 2).
#'   A warning is issued if any part has fewer items.
#' @param degree Integer. Degree of the polynomial (default 2).
#' @param bin.score Integer. Number of quantile groups (e.g., 5 for quintiles).
#'   If `NULL` (default), no binning is performed. If provided, the function
#'   returns a data frame `binned.CSEM` with average CSEM per quantile group.
#' @param full.range Logical. If `TRUE`, evaluates the smoothed CSEM for every
#'   integer score from `score.range[1]` to `score.range[2]` (requires `score.range`).
#' @param ci Logical. If `TRUE`, compute confidence intervals for the true score.
#' @param conf.level Numeric. Confidence level (default 0.95).
#' @param digits Integer. Rounding for output.
#' @param score.range Numeric vector length 2. Required if `full.range = TRUE`.
#'   Also used to truncate confidence intervals (if `ci = TRUE`).
#' @param na.rm Logical. Remove rows with missing values.
#'
#' @return A list with elements:
#' \item{CSEM}{data.frame with columns `score`, `n`, `CSEM.smooth`, and
#'   if `ci = TRUE`, `lwr.ci` and `upr.ci`.}
#' \item{binned.CSEM}{(if `bin.score` is provided) data.frame with quantile groups:
#'   `group`, `range`, `n`, `mean_score`, `CSEM.mean`, and intervals if `ci = TRUE`.}
#'
#' @references
#' Mollenkopf, W. G. (1949). Variation of the standard error of measurement.
#'   *Psychometrika*, 14(3), 189–229.
#' Feldt, L. S., & Qualls, A. L. (1996). Estimation of measurement error variance
#'   at specific score levels. *Journal of Educational Measurement*, 33(2), 141–156.
#'
#' @examples
#' # Examples will be added by the user.
#'
#' @export
csemMF <- function(data,
                   n.parts = NULL,
                   part_items = NULL,
                   min.items.per.part = 2,
                   degree = 2,
                   bin.score = NULL,
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

  # --- Binning (si bin.score es un entero) ---
  binned_df <- NULL
  if (!is.null(bin.score)) {
    if (!is.numeric(bin.score) || length(bin.score) != 1 || bin.score < 2)
      stop("bin.score must be an integer >= 2.")
    # Tomar solo filas con n > 0 (observadas)
    obs_df <- result_df[result_df$n > 0, c("score", "CSEM.smooth")]
    if (nrow(obs_df) == 0) stop("No observed scores for binning.")
    # Expandir a puntajes por persona (ponderado por frecuencia)
    person_scores <- rep(obs_df$score, times = result_df$n[result_df$n > 0])
    q <- stats::quantile(person_scores, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
    q <- unique(q)
    groups <- cut(person_scores, breaks = q, include.lowest = TRUE, right = TRUE)
    group_levels <- levels(groups)
    bin_list <- list()
    for (i in seq_along(group_levels)) {
      idx_in_group <- which(groups == group_levels[i])
      scores_in_group <- unique(person_scores[idx_in_group])
      sub_df <- obs_df[obs_df$score %in% scores_in_group, , drop = FALSE]
      if (nrow(sub_df) == 0) next
      csem_mean <- mean(sub_df$CSEM.smooth, na.rm = TRUE)
      range_str <- paste0(min(scores_in_group), "-", max(scores_in_group))
      n_persons <- length(idx_in_group)
      mean_score <- mean(person_scores[idx_in_group])
      bin_list[[i]] <- data.frame(
        group = i,
        range = range_str,
        n = n_persons,
        mean_score = mean_score,
        CSEM.mean = csem_mean,
        stringsAsFactors = FALSE
      )
    }
    binned_df <- do.call(rbind, bin_list)
    if (ci && !is.null(binned_df)) {
      z <- stats::qnorm(1 - (1 - conf.level) / 2)
      lwr_bin <- binned_df$mean_score - z * binned_df$CSEM.mean
      upr_bin <- binned_df$mean_score + z * binned_df$CSEM.mean
      lwr_bin <- pmax(lwr_bin, score_min_teo)
      upr_bin <- pmin(upr_bin, score_max_teo)
      binned_df$lwr.ci <- round(lwr_bin, digits)
      binned_df$upr.ci <- round(upr_bin, digits)
    }
  }

  out <- list(CSEM = result_df)
  if (!is.null(binned_df)) out$binned.CSEM <- binned_df
  return(out)
}
