#' Feldt & Qualls (1996) method for Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Implements the improved method by Feldt & Qualls (1996) for estimating conditional
#' measurement error variance. The test is partitioned into `n.parts` (or as specified
#' by `part_items`). For each person, an unbiased estimate of the error variance
#' is computed from the part-test scores. These individual estimates are then
#' averaged within each total score group to obtain the CSEM.
#'
#' This function does not apply smoothing; for the smoothed version (polynomial
#' regression) see `csemMF`.
#'
#' @param data A data frame or matrix with item responses (subjects in rows,
#'   items in columns). Items can be dichotomous or polytomous.
#' @param n.parts Integer. Number of parts into which the test will be split
#'   (by column order, as balanced as possible). Ignored if `part_items` is provided.
#'   Default is `NULL`, which with `part_items = NULL` sets `n.parts = ncol(data)`
#'   (each item as a part).
#' @param part_items Optional list. Each element is a character vector of column names
#'   or an integer vector of column indices defining the items in that part.
#'   If provided, `n.parts` is ignored.
#' @param min.items.per.part Integer. Minimum number of items per part (default 2).
#'   A warning is issued if any part has fewer items.
#' @param bin.score Integer. Number of quantile groups (e.g., 5 for quintiles).
#'   If `NULL` (default), CSEM is reported for each observed score (with `n >= 2`).
#' @param full.range Logical. If `TRUE`, report all integer scores from
#'   `score.range[1]` to `score.range[2]` (requires `score.range`).
#'   Missing values are filled with the last valid observation (carried forward).
#'   Default `FALSE`.
#' @param ci Logical. If `TRUE`, compute confidence intervals for the true score.
#'   Default `FALSE`.
#' @param conf.level Numeric. Confidence level (default 0.95).
#' @param digits Integer. Rounding for CSEM and confidence limits. Default 3.
#' @param score.range Numeric vector of length 2 (min, max). Required for `full.range = TRUE`.
#'   Also used to truncate confidence intervals (if `ci = TRUE`). If `NULL`, the
#'   observed range is used for truncation.
#' @param na.rm Logical. If `TRUE` (default), removes rows with any missing values.
#'
#' @return A list with elements:
#' \item{CSEM}{data.frame with columns `score`, `n`, and `CSEM` (raw).
#'   If `ci = TRUE`, also `lwr.ci` and `upr.ci` (truncated to possible score range).}
#' \item{binned.CSEM}{(if `bin.score` is integer) data.frame with quantile groups:
#'   `group`, `range`, `n`, `mean_score`, `CSEM.mean`, and intervals if `ci = TRUE`.}
#'
#' @references
#' Feldt, L. S., & Qualls, A. L. (1996). Estimation of measurement error variance
#' at specific score levels. *Journal of Educational Measurement*, 33(2), 141–156.
#'
#' @export
csemFeldtQualls <- function(data,
                            n.parts = NULL,
                            part_items = NULL,
                            min.items.per.part = 2,
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

  n_persons <- nrow(data)
  J <- ncol(data)                 # número total de ítems
  if (n_persons < 2) stop("At least 2 persons required.")
  if (J < 2) stop("At least 2 items required.")

  # --- Construcción de partes (part_test scores) ---
  if (!is.null(part_items)) {
    # Partición manual
    if (!is.list(part_items)) stop("`part_items` must be a list.")
    n.parts <- length(part_items)
    part_scores <- matrix(NA, nrow = n_persons, ncol = n.parts)
    colnames(part_scores) <- paste0("Part", 1:n.parts)
    items_per_part <- numeric(n.parts)
    for (j in 1:n.parts) {
      cols <- part_items[[j]]
      if (is.character(cols)) {
        if (!all(cols %in% colnames(data)))
          stop("Some column names in part_items[[", j, "]] not found in data.")
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
    if (any(items_per_part < min.items.per.part)) {
      warning("Some parts have fewer than ", min.items.per.part,
              " items. Estimates may be unstable.")
    }
    # factor d = J / k, donde k es el número de items por parte (no constante)
    # Para la fórmula (13), se necesita k (items por parte) y luego d = J/k.
    # Como las partes pueden tener distinto k, usaremos el promedio? La teoría
    # asume partes iguales. Para ser fieles, advertimos si no son iguales.
    if (length(unique(items_per_part)) > 1) {
      warning("Parts have unequal number of items. The method assumes equal length; results may be biased.")
    }
    # Usamos k = J / n.parts (no el real) porque la fórmula original espera partes iguales.
    # Pero para una implementación robusta, usaremos el k medio redondeado? Mejor seguir
    # la recomendación: si las partes son desiguales, el usuario debe reordenar.
    # Por simplicidad, usaremos k = J / n.parts (puede no ser entero) y d = n.parts.
    # Esto es lo que hace Feldt & Qualls cuando omiten items para balancear.
    k <- J / n.parts
    d <- J / k   # = n.parts
  } else {
    # Partición automática
    if (is.null(n.parts)) n.parts <- J   # por defecto, cada ítem es una parte
    if (n.parts < 2) stop("n.parts must be at least 2.")
    if (n.parts > J) stop("n.parts cannot exceed number of items.")
    # Crear partes balanceadas por orden de columnas
    idx_split <- split(1:J, cut(1:J, breaks = n.parts, labels = FALSE))
    part_scores <- matrix(NA, nrow = n_persons, ncol = n.parts)
    items_per_part <- numeric(n.parts)
    for (j in 1:n.parts) {
      cols <- idx_split[[j]]
      part_scores[, j] <- rowSums(data[, cols, drop = FALSE], na.rm = TRUE)
      items_per_part[j] <- length(cols)
    }
    if (any(items_per_part < min.items.per.part)) {
      warning("Some parts have fewer than ", min.items.per.part,
              " items. Estimates may be unstable.")
    }
    k <- J / n.parts   # puede no ser entero
    d <- J / k         # = n.parts (exactamente)
  }

  # --- Calcular Y_i para cada persona (ecuación 13, Feldt & Qualls 1996) ---
  # Y_i = d * [ sum_j ( (X_ij - barX_i) - (barX_j - M) )^2 / (n.parts - 1) ]
  # donde barX_i = media de las partes de la persona i,
  # barX_j = media de la parte j (sobre personas),
  # M = media de los barX_j.

  Xij <- part_scores
  barX_i <- rowMeans(Xij, na.rm = TRUE)           # media por persona
  barX_j <- colMeans(Xij, na.rm = TRUE)           # media por parte
  M <- mean(barX_j)                               # media de las medias de partes

  # Matriz de desviaciones: (X_ij - barX_i) - (barX_j - M)
  dev <- sweep(Xij, 1, barX_i, "-")               # X_ij - barX_i
  dev <- sweep(dev, 2, barX_j - M, "-")           # (X_ij - barX_i) - (barX_j - M)

  # Suma de cuadrados por persona
  SS_i <- rowSums(dev^2, na.rm = TRUE)
  var_adj_i <- SS_i / (n.parts - 1)               # varianza ajustada
  Y_i <- d * var_adj_i                            # estimación de error varianza para el test completo

  # --- CSEM por puntaje total (raw) ---
  total <- rowSums(data)   # puntaje total original (suma de ítems, no de partes)
  unique_scores <- sort(unique(total))
  raw_list <- list()
  for (s in unique_scores) {
    idx <- which(total == s)
    n_s <- length(idx)
    if (n_s >= 2) {
      csem_raw <- sqrt(mean(Y_i[idx]))      # promedio de Y_i, luego raíz
    } else {
      csem_raw <- NA_real_
    }
    raw_list[[length(raw_list)+1]] <- data.frame(score = s, n = n_s, CSEM = csem_raw,
                                                 stringsAsFactors = FALSE)
  }
  raw_df <- do.call(rbind, raw_list)
  raw_df <- raw_df[!is.na(raw_df$CSEM), , drop = FALSE]
  raw_df$CSEM <- round(raw_df$CSEM, digits)

  # --- Definir rango de puntajes (para truncar CI y full.range) ---
  if (!is.null(score.range)) {
    if (!is.numeric(score.range) || length(score.range) != 2)
      stop("score.range must be a numeric vector of length 2.")
    score_min_teo <- score.range[1]
    score_max_teo <- score.range[2]
  } else {
    score_min_teo <- min(total, na.rm = TRUE)
    score_max_teo <- max(total, na.rm = TRUE)
  }

  # --- Función auxiliar: last observation carried forward (locf) ---
  na_locf <- function(x) {
    idx <- !is.na(x)
    if (sum(idx) == 0) return(x)
    last_val <- x[idx][1]
    for (i in seq_along(x)) {
      if (!is.na(x[i])) last_val <- x[i]
      else x[i] <- last_val
    }
    return(x)
  }

  # --- Si full.range = TRUE, expandir a todos los enteros ---
  if (full.range) {
    if (is.null(score.range))
      stop("full.range = TRUE requires 'score.range'.")
    all_scores <- seq(score.range[1], score.range[2], by = 1)
    full_n <- sapply(all_scores, function(s) sum(total == s))
    csem_map <- setNames(raw_df$CSEM, raw_df$score)
    csem_full <- csem_map[as.character(all_scores)]
    csem_full[is.na(csem_full)] <- NA_real_
    result_df <- data.frame(score = all_scores, n = full_n, CSEM = csem_full,
                            stringsAsFactors = FALSE)
    # Rellenar NA hacia abajo
    result_df$CSEM <- na_locf(result_df$CSEM)
  } else {
    result_df <- raw_df
  }

  # --- Intervalos de confianza ---
  if (ci) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    lwr <- result_df$score - z * result_df$CSEM
    upr <- result_df$score + z * result_df$CSEM
    lwr <- pmax(lwr, score_min_teo)
    upr <- pmin(upr, score_max_teo)
    result_df$lwr.ci <- round(lwr, digits)
    result_df$upr.ci <- round(upr, digits)
  }

  # --- Agrupación por cuantiles (bin.score) ---
  binned_df <- NULL
  if (!is.null(bin.score)) {
    # Usamos raw_df (CSEM por cada score observado)
    q <- stats::quantile(total, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
    q <- unique(q)
    groups <- cut(total, breaks = q, include.lowest = TRUE, right = TRUE)
    group_levels <- levels(groups)
    bin_list <- list()
    for (i in seq_along(group_levels)) {
      idx_in_group <- which(groups == group_levels[i])
      scores_in_group <- unique(total[idx_in_group])
      sub_df <- raw_df[raw_df$score %in% scores_in_group, , drop = FALSE]
      if (nrow(sub_df) == 0) next
      csem_mean <- mean(sub_df$CSEM, na.rm = TRUE)
      range_str <- paste0(min(scores_in_group), "-", max(scores_in_group))
      n_persons <- length(idx_in_group)
      mean_score <- mean(total[idx_in_group])
      bin_list[[i]] <- data.frame(group = i, range = range_str, n = n_persons,
                                  mean_score = mean_score, CSEM.mean = csem_mean,
                                  stringsAsFactors = FALSE)
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
