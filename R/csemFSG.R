#' ANOVA method for Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Estimates conditional standard errors of measurement using the variance-components
#' approach described in Feldt, Steffen, & Gupta (1985). This is the method
#' referred to as "ANOVA" in the JASP Reliability module.
#'
#' For each unique total score (or quantile group), the CSEM is computed as:
#' \deqn{CSEM = \sqrt{ \frac{J}{J-1} \sum_{j=1}^J s_j^2 }}
#' where \eqn{J} is the number of items and \eqn{s_j^2} is the sample variance
#' of item \eqn{j} within the group (using divisor \eqn{n_g - 1}).
#'
#' The function offers polynomial smoothing (`smooth = TRUE`) and expansion to
#' a full range of integer scores (`full.range = TRUE`). It also allows
#' aggregation into quantile groups (`bin.score`).
#'
#' @param data A data frame or matrix with item responses (subjects in rows,
#'   items in columns). Items can be dichotomous or polytomous.
#' @param bin.score integer. Number of quantile groups (e.g., 5 for quintiles).
#'   If `NULL` (default), CSEM is reported for each observed score (with `n >= 2`).
#' @param smooth logical. If `TRUE`, applies polynomial smoothing to the squared
#'   CSEM estimates. Default = `FALSE`.
#' @param degree integer. Polynomial degree used when `smooth = TRUE`. Default = 2.
#' @param full.range logical. If `TRUE` and `smooth = TRUE`, evaluates the smoothed
#'   CSEM for every integer score from `score.range[1]` to `score.range[2]`.
#'   Requires `score.range`. Default = `FALSE`.
#' @param ci logical. If `TRUE`, compute confidence intervals for the true score.
#'   Default = `FALSE`.
#' @param conf.level numeric. Confidence level for intervals (default 0.95).
#' @param digits integer. Rounding for CSEM and confidence limits. Default = 3.
#' @param score.range numeric vector of length 2 (min, max). Required when
#'   `full.range = TRUE`. Defines the theoretical score range (e.g., c(0,36)).
#'   Also used to truncate confidence intervals (if ci=TRUE). If `NULL`, the
#'   observed range is used for truncation.
#' @param na.rm logical. If `TRUE` (default), removes rows with any missing values.
#'
#' @return A list with elements:
#' \item{CSEM}{data.frame with columns `score`, `n`,
#'   and either `CSEM` (raw) or `CSEM.smooth` (smoothed). If `ci = TRUE`, also
#'   `lwr.ci` and `upr.ci` (truncated to possible score range).}
#' \item{binned.CSEM}{(if `bin.score` is integer) data.frame with quantile groups:
#'   `group`, `range`, `n`, `mean_score`, `CSEM.mean`, and intervals if `ci = TRUE`.}
#'
#' @references
#' Feldt, L. S., Steffen, M., & Gupta, N. C. (1985).
#' A comparison of five methods for estimating the standard error of measurement
#' at specific score levels. \emph{Applied Psychological Measurement}, 9(4), 351–361.
#'
#' @export
csemFSG <- function(data,
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
  if (!is.data.frame(data) && !is.matrix(data))
    stop("`data` must be a data frame or matrix.")
  data <- as.data.frame(data)

  if (na.rm) data <- stats::na.omit(data)
  if (anyNA(data)) stop("Missing values present. Set na.rm = TRUE to remove them.")

  n_persons <- nrow(data)
  n_items   <- ncol(data)
  if (n_persons < 2) stop("At least 2 persons required.")
  if (n_items < 2) stop("At least 2 items required.")

  total <- rowSums(data)

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

  # --- Helper: last observation carried forward (for raw CSEM with full.range) ---
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

  # --- Raw CSEM for each unique score (only n >= 2) ---
  unique_scores <- sort(unique(total))
  raw_df <- data.frame(score = unique_scores,
                       n = as.integer(table(total)[as.character(unique_scores)]))
  raw_df$CSEM <- sapply(unique_scores, function(s) {
    idx <- which(total == s)
    if (length(idx) < 2) return(NA_real_)
    s2j <- apply(data[idx, , drop = FALSE], 2, stats::var)
    sqrt( (n_items / (n_items - 1)) * sum(s2j) )
  })
  raw_df <- raw_df[!is.na(raw_df$CSEM), , drop = FALSE]
  raw_df$CSEM <- round(raw_df$CSEM, digits)

  # --- Modo sin suavizamiento ---
  if (!smooth) {
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
      # Rellenar NA hacia abajo para mejorar presentación
      result_df$CSEM <- na_locf(result_df$CSEM)
    } else {
      result_df <- raw_df
    }
  } else {
    # --- Modo con suavizamiento ---
    if (nrow(raw_df) < degree + 1)
      stop("Not enough unique scores to fit polynomial of degree ", degree)
    fit <- lm(CSEM^2 ~ poly(score, degree, raw = TRUE), data = raw_df)

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
  }

  # --- Intervalos de confianza (si ci = TRUE) ---
  if (ci) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    if (smooth) {
      csem_vals <- result_df$CSEM.smooth
    } else {
      csem_vals <- result_df$CSEM
    }
    lwr <- result_df$score - z * csem_vals
    upr <- result_df$score + z * csem_vals
    # Truncar al rango de puntajes
    lwr <- pmax(lwr, score_min_teo)
    upr <- pmin(upr, score_max_teo)
    result_df$lwr.ci <- round(lwr, digits)
    result_df$upr.ci <- round(upr, digits)
  }

  # --- Agrupación por cuantiles (bin.score) ---
  binned_df <- NULL
  if (!is.null(bin.score)) {
    # Obtener CSEM para cada score único observado (raw o smooth)
    if (smooth) {
      # Predicciones del modelo para los scores observados
      scores_obs <- sort(unique(total))
      pred_obs <- predict(fit, newdata = data.frame(score = scores_obs))
      pred_obs <- pmax(pred_obs, 0)
      csem_obs <- sqrt(pred_obs)
      temp_df <- data.frame(score = scores_obs, CSEM = csem_obs)
    } else {
      # Usamos raw_df (que ya tiene CSEM)
      temp_df <- raw_df
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
