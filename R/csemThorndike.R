#' Thorndike's method for Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Simplified version without manual merging of small groups.
#' Use `smooth = TRUE` for stable estimates across the whole score range.
#' When `smooth = FALSE` and `full.range = TRUE`, missing values (NA) are
#' filled with the last valid observation (carried forward) for better readability.
#'
#' @param half1 data.frame/matrix with first half items.
#' @param half2 data.frame/matrix with second half items.
#' @param bin.score integer (number of quantile groups) or NULL (individual scores).
#' @param smooth logical. If TRUE, applies polynomial smoothing.
#' @param degree integer. Polynomial degree (used if smooth=TRUE).
#' @param full.range logical. If TRUE, report all integer scores from
#'   `score.range[1]` to `score.range[2]` (requires `score.range`).
#' @param ci logical. If TRUE, compute confidence intervals for true score.
#' @param conf.level numeric. Confidence level (default 0.95).
#' @param digits integer. Rounding for output.
#' @param score.range numeric vector of length 2 (min, max). Required for full.range=TRUE.
#'   Also used to truncate confidence intervals (if ci=TRUE).
#' @param na.rm logical. Not used (kept for compatibility).
#'
#' @return A list with elements:
#'   \item{CSEM}{data.frame with columns `score`, `n`,
#'     and either `CSEM` (raw) or `CSEM.smooth` (smoothed). If ci=TRUE, also
#'     `lwr.ci` and `upr.ci` (truncated to possible score range).}
#'   \item{binned.CSEM}{(if bin.score is integer) data.frame with quantile groups.}
#'
#' @references
#'Thorndike, R. L. (1951). Reliability. In E. F. Lindquist (Ed.), Educational measurement
#'(pp. 560–620). American Council on Education.
#'
#'Lee, W., & Harris, D. J. (2025). Reliability in educational measurement. In L. L. Cook & M. J. Pitoniak (Eds.),
#'Educational measurement (5th ed., pp. 277–381). Oxford University Press. \doi{10.1093/oso/9780197654965.003.0005}
#'
#'
#'
#' \donttest{
#' ## Load data
#' library(EFA.dimensions)
#' data("data_RSE")
#'
#' ## Recode negative items
#' data_RSE[c("Q3", "Q5", "Q8", "Q9", "Q10")] <- 5 - data_RSE[c("Q3", "Q5", "Q8", "Q9", "Q10")]
#'
#' ## Choosing split by difficulty criteria
#' RSE.namesHalf <- checkSplit(data = data_RSE, method = "difficulty")
#'
#' RSE.namesHalf$half1
#' RSE.namesHalf$half2
#'
#' #' # Thorndike csem, basic ouput
#' csemThorndike(half1 = data_RSE[, RSE.namesHalf$half1],
#' half2 = data_RSE[, RSE.namesHalf$half2],
#' smooth = F,
#' ci = F)
#'
#' # Thorndike csem, smoothing and binned score
#' csemThorndike(half1 = data_RSE[, RSE.namesHalf$half1],
#' half2 = data_RSE[, RSE.namesHalf$half2],
#' smooth = T,
#' degree = 2,
#' ci = F,
#' bin.score = 5)
#'
#' # Thorndike csem, smoothing, binned score, and confidence interval
#' csemThorndike(half1 = data_RSE[, RSE.namesHalf$half1],
#' half2 = data_RSE[, RSE.namesHalf$half2],
#' smooth = T,
#' degree = 2,
#' ci = T,
#' conf.level = .90,
#' bin.score = 5)
#' }
#'
#' @export
csemThorndike <- function(half1, half2,
                          bin.score = NULL,
                          smooth = FALSE,
                          degree = 2,
                          full.range = FALSE,
                          ci = FALSE,
                          conf.level = 0.95,
                          digits = 3,
                          score.range = NULL,
                          na.rm = TRUE) {

  # --- Validaciones ---
  if (missing(half1) || missing(half2))
    stop("Both 'half1' and 'half2' must be provided.")
  half1 <- as.data.frame(half1)
  half2 <- as.data.frame(half2)
  if (nrow(half1) != nrow(half2))
    stop("'half1' and 'half2' must have the same number of rows.")

  total1 <- rowSums(half1, na.rm = TRUE)
  total2 <- rowSums(half2, na.rm = TRUE)
  total <- total1 + total2
  diff <- total1 - total2

  # --- Definir rango de puntajes (para truncar CI) ---
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

  # --- Modo sin suavizamiento (smooth = FALSE) ---
  if (!smooth) {
    # Obtener CSEM crudo para cada puntaje único con n>=2
    unique_scores <- sort(unique(total))
    raw_list <- list()
    for (s in unique_scores) {
      idx <- which(total == s)
      n_s <- length(idx)
      if (n_s >= 2) {
        csem_raw <- sd(diff[idx])
      } else {
        csem_raw <- NA_real_
      }
      raw_list[[length(raw_list)+1]] <- data.frame(score = s, n = n_s, CSEM = csem_raw,
                                                   stringsAsFactors = FALSE)
    }
    raw_df <- do.call(rbind, raw_list)
    raw_df <- raw_df[!is.na(raw_df$CSEM), , drop = FALSE]
    raw_df$CSEM <- round(raw_df$CSEM, digits)

    if (full.range) {
      if (is.null(score.range))
        stop("full.range = TRUE requires 'score.range'.")
      all_scores <- seq(score.range[1], score.range[2], by = 1)
      full_n <- sapply(all_scores, function(s) sum(total == s))
      csem_map <- setNames(raw_df$CSEM, raw_df$score)
      csem_full <- csem_map[as.character(all_scores)]
      csem_full[is.na(csem_full)] <- NA_real_
      raw_df <- data.frame(score = all_scores, n = full_n, CSEM = csem_full,
                           stringsAsFactors = FALSE)
      # Rellenar NA hacia abajo para mejorar presentación
      raw_df$CSEM <- na_locf(raw_df$CSEM)
    }

    # Intervalos de confianza (si se solicitan)
    if (ci) {
      z <- stats::qnorm(1 - (1 - conf.level) / 2)
      lwr <- raw_df$score - z * raw_df$CSEM
      upr <- raw_df$score + z * raw_df$CSEM
      # Truncar al rango de puntajes (teórico u observado)
      lwr <- pmax(lwr, score_min_teo)
      upr <- pmin(upr, score_max_teo)
      raw_df$lwr.ci <- round(lwr, digits)
      raw_df$upr.ci <- round(upr, digits)
    }

    result_df <- raw_df
  } else {
    # --- Modo con suavizamiento (smooth = TRUE) ---
    unique_scores <- sort(unique(total))
    raw_list <- list()
    for (s in unique_scores) {
      idx <- which(total == s)
      n_s <- length(idx)
      if (n_s >= 2) {
        csem_raw <- sd(diff[idx])
      } else {
        csem_raw <- NA_real_
      }
      raw_list[[length(raw_list)+1]] <- data.frame(score = s, n = n_s, CSEM.raw = csem_raw,
                                                   stringsAsFactors = FALSE)
    }
    raw_df <- do.call(rbind, raw_list)
    raw_df <- raw_df[!is.na(raw_df$CSEM.raw), , drop = FALSE]

    if (nrow(raw_df) < degree + 1)
      stop("Not enough valid scores to fit polynomial of degree ", degree)

    fit <- lm(CSEM.raw^2 ~ poly(score, degree, raw = TRUE), data = raw_df)

    if (full.range) {
      if (is.null(score.range))
        stop("full.range = TRUE requires 'score.range'.")
      all_scores <- seq(score.range[1], score.range[2], by = 1)
      pred_var <- predict(fit, newdata = data.frame(score = all_scores))
      pred_var <- pmax(pred_var, 0)
      csem_smooth <- sqrt(pred_var)
      full_n <- sapply(all_scores, function(s) sum(total == s))
      result_df <- data.frame(score = all_scores, n = full_n,
                              CSEM.smooth = round(csem_smooth, digits),
                              stringsAsFactors = FALSE)
    } else {
      pred_var <- predict(fit, newdata = data.frame(score = raw_df$score))
      pred_var <- pmax(pred_var, 0)
      csem_smooth <- sqrt(pred_var)
      result_df <- data.frame(score = raw_df$score, n = raw_df$n,
                              CSEM.smooth = round(csem_smooth, digits),
                              stringsAsFactors = FALSE)
    }

    if (ci) {
      z <- stats::qnorm(1 - (1 - conf.level) / 2)
      lwr <- result_df$score - z * result_df$CSEM.smooth
      upr <- result_df$score + z * result_df$CSEM.smooth
      lwr <- pmax(lwr, score_min_teo)
      upr <- pmin(upr, score_max_teo)
      result_df$lwr.ci <- round(lwr, digits)
      result_df$upr.ci <- round(upr, digits)
    }
  }

  # --- Agrupación por cuantiles (bin.score) ---
  binned_df <- NULL
  if (!is.null(bin.score)) {
    # Necesitamos los valores de CSEM (raw o smooth) para cada score único observado
    if (smooth) {
      # Para suavizado, usamos predicciones del modelo para los scores observados
      scores_obs <- sort(unique(total))
      pred_obs <- predict(fit, newdata = data.frame(score = scores_obs))
      pred_obs <- pmax(pred_obs, 0)
      csem_obs <- sqrt(pred_obs)
      temp_df <- data.frame(score = scores_obs, CSEM = csem_obs)
    } else {
      # Usamos los datos crudos (raw_df) que ya tienen CSEM
      if (full.range) {
        # Tomamos solo los scores observados (con n>0) de la tabla expandida
        temp_df <- result_df[result_df$n > 0 & !is.na(result_df$CSEM), c("score", "CSEM")]
      } else {
        temp_df <- result_df[, c("score", "CSEM")]
      }
    }
    # Cuantiles sobre las personas
    q <- stats::quantile(total, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
    q <- unique(q)
    groups <- cut(total, breaks = q, include.lowest = TRUE, right = TRUE)
    group_levels <- levels(groups)
    bin_list <- list()
    for (i in seq_along(group_levels)) {
      idx_in_group <- which(groups == group_levels[i])
      scores_in_group <- unique(total[idx_in_group])
      sub_df <- temp_df[temp_df$score %in% scores_in_group, , drop = FALSE]
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
