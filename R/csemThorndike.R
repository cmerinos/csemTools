#' Thorndike's method for Conditional Standard Error of Measurement (CSEM)
#'
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
#'
#' @return A list with elements:
#'   \item{CSEM}{data.frame with columns `score`, `n` (number of subjects),
#'     and either `CSEM` (raw) or `CSEM.smooth` (smoothed).}
#'   \item{binned.CSEM}{(if bin.score is integer) data.frame with quantile groups.}
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
                          score.range = NULL) {

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

  # --- Función auxiliar: last observation carried forward (locf) ---
  na_locf <- function(x) {
    idx <- !is.na(x)
    if (sum(idx) == 0) return(x)
    last_val <- x[idx][1]
    for (i in seq_along(x)) {
      if (!is.na(x[i])) {
        last_val <- x[i]
      } else {
        x[i] <- last_val
      }
    }
    return(x)
  }

  # --- Modo sin suavizamiento (smooth = FALSE) ---
  if (!smooth) {
    # Obtener CSEM crudo solo para puntajes con n>=2
    get_raw_csem <- function() {
      unique_scores <- sort(unique(total))
      res <- data.frame(score = unique_scores,
                        n = as.integer(table(total)[as.character(unique_scores)]))
      res$CSEM <- sapply(unique_scores, function(s) {
        idx <- which(total == s)
        if (length(idx) >= 2) sd(diff[idx]) else NA_real_
      })
      res <- res[!is.na(res$CSEM), , drop = FALSE]
      res$CSEM <- round(res$CSEM, digits)
      return(res)
    }

    raw_df <- get_raw_csem()

    # Expandir a rango completo si se solicita
    if (full.range) {
      if (is.null(score.range))
        stop("full.range = TRUE requires 'score.range' (e.g., c(0,36)).")
      all_scores <- seq(score.range[1], score.range[2], by = 1)
      full_n <- sapply(all_scores, function(s) sum(total == s))
      csem_map <- setNames(raw_df$CSEM, raw_df$score)
      csem_full <- csem_map[as.character(all_scores)]
      csem_full[is.na(csem_full)] <- NA_real_
      raw_df <- data.frame(score = all_scores, n = full_n, CSEM = csem_full,
                           stringsAsFactors = FALSE)
    }

    # Intervalos de confianza (si se piden)
    if (ci) {
      z <- stats::qnorm(1 - (1 - conf.level) / 2)
      lwr <- raw_df$score - z * raw_df$CSEM
      upr <- raw_df$score + z * raw_df$CSEM
      raw_df$lwr.ci <- round(lwr, digits)
      raw_df$upr.ci <- round(upr, digits)
    }

    # Rellenar NAs con el último valor válido (mejora de presentación)
    if (full.range) {
      raw_df$CSEM <- na_locf(raw_df$CSEM)
      if (ci) {
        raw_df$lwr.ci <- na_locf(raw_df$lwr.ci)
        raw_df$upr.ci <- na_locf(raw_df$upr.ci)
      }
    }

    # Agrupación por cuantiles (bin.score)
    binned_df <- NULL
    if (!is.null(bin.score)) {
      # Usamos solo las filas con CSEM no NA para promediar
      valid <- raw_df[!is.na(raw_df$CSEM), ]
      if (nrow(valid) == 0) stop("No valid scores for binning.")
      q <- stats::quantile(total, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
      q <- unique(q)
      groups <- cut(total, breaks = q, include.lowest = TRUE, right = TRUE)
      group_levels <- levels(groups)
      bin_list <- list()
      for (i in seq_along(group_levels)) {
        idx_in_group <- which(groups == group_levels[i])
        scores_in_group <- unique(total[idx_in_group])
        sub_df <- valid[valid$score %in% scores_in_group, , drop = FALSE]
        if (nrow(sub_df) == 0) next
        csem_mean <- mean(sub_df$CSEM)
        range_str <- paste0(min(scores_in_group), "-", max(scores_in_group))
        n_persons <- length(idx_in_group)
        mean_score <- mean(total[idx_in_group])
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
        binned_df$lwr.ci <- round(lwr_bin, digits)
        binned_df$upr.ci <- round(upr_bin, digits)
      }
    }

    out <- list(CSEM = raw_df)
    if (!is.null(binned_df)) out$binned.CSEM <- binned_df
    return(out)
  }

  # --- Modo con suavizamiento (smooth = TRUE) ---
  # Obtener estimaciones crudas para todos los puntajes con n>=2
  unique_scores <- sort(unique(total))
  raw_list <- list()
  for (s in unique_scores) {
    idx <- which(total == s)
    n_s <- length(idx)
    if (n_s >= 2) {
      csem_val <- sd(diff[idx])
    } else {
      csem_val <- NA_real_
    }
    raw_list[[length(raw_list)+1]] <- data.frame(
      score = s, n = n_s, CSEM.raw = csem_val, stringsAsFactors = FALSE
    )
  }
  raw_df <- do.call(rbind, raw_list)
  raw_df <- raw_df[!is.na(raw_df$CSEM.raw), , drop = FALSE]

  if (nrow(raw_df) < degree + 1)
    stop("Not enough valid scores to fit polynomial of degree ", degree)

  # Ajuste polinómico sobre CSEM^2
  fit <- lm(CSEM.raw^2 ~ poly(score, degree, raw = TRUE), data = raw_df)

  # Rango de evaluación
  if (full.range) {
    if (is.null(score.range))
      stop("full.range = TRUE requires 'score.range' (e.g., c(0,36)).")
    eval_scores <- seq(score.range[1], score.range[2], by = 1)
  } else {
    eval_scores <- raw_df$score
  }

  pred_var <- predict(fit, newdata = data.frame(score = eval_scores))
  pred_var <- pmax(pred_var, 0)
  csem_smooth <- sqrt(pred_var)

  if (full.range) {
    n_vals <- sapply(eval_scores, function(s) sum(total == s))
  } else {
    n_vals <- raw_df$n
  }

  result_df <- data.frame(
    score = eval_scores,
    n = n_vals,
    CSEM.smooth = round(csem_smooth, digits),
    stringsAsFactors = FALSE
  )

  if (ci) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    lwr <- result_df$score - z * result_df$CSEM.smooth
    upr <- result_df$score + z * result_df$CSEM.smooth
    result_df$lwr.ci <- round(lwr, digits)
    result_df$upr.ci <- round(upr, digits)
  }

  # Agrupación por cuantiles usando valores suavizados
  binned_df <- NULL
  if (!is.null(bin.score)) {
    q <- stats::quantile(total, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
    q <- unique(q)
    groups <- cut(total, breaks = q, include.lowest = TRUE, right = TRUE)
    group_levels <- levels(groups)
    bin_list <- list()
    for (i in seq_along(group_levels)) {
      idx_in_group <- which(groups == group_levels[i])
      scores_in_group <- unique(total[idx_in_group])
      pred_grp <- predict(fit, newdata = data.frame(score = scores_in_group))
      pred_grp <- pmax(pred_grp, 0)
      csem_pred <- sqrt(pred_grp)
      csem_mean <- mean(csem_pred, na.rm = TRUE)
      range_str <- paste0(min(scores_in_group), "-", max(scores_in_group))
      n_persons <- length(idx_in_group)
      mean_score <- mean(total[idx_in_group])
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
      binned_df$lwr.ci <- round(lwr_bin, digits)
      binned_df$upr.ci <- round(upr_bin, digits)
    }
  }

  out <- list(CSEM = result_df)
  if (!is.null(binned_df)) out$binned.CSEM <- binned_df
  return(out)
}
