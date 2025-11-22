#' Feldt's method for Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Calculates the Conditional Standard Error of Measurement (CSEM) following the
#' Feldt, Steffen, & Gupta (1985) variance-components approach. The CSEM is
#' estimated at each score level or score band using item variances within groups.
#'
#' @param data A data frame or matrix with item responses (subjects in rows, items in columns).
#'   Items can be dichotomous or polytomous (Likert). All items contribute to the total score.
#' @param score.group Character. Method for grouping total scores. One of:
#'   \itemize{
#'     \item "all": CSEM for each unique total score.
#'     \item "deciles": 10 bands by score quantiles.
#'     \item "quartiles": 4 bands by score quantiles.
#'     \item "quintiles": 5 bands by score quantiles.
#'     \item "k": custom number of bands defined by \code{k}.
#'   }
#'   Default is "all".
#' @param k Integer. Number of bands if \code{score.group = "k"} (>= 2).
#' @param aggregate Character. How to represent the central score of each band:
#'   \code{"mean"}, \code{"median"}, or \code{"midpoint"} (midpoint of observed score range in the band).
#'   Default is "midpoint".
#' @param min.n Integer. Minimum band size. Bands with \code{n < min.n} are merged with adjacent bands
#'   until all bands have at least \code{min.n} subjects. Default = 5.
#' @param merge.method Character. Rule to merge undersized bands: \code{"nearest"} (default),
#'   \code{"left"}, or \code{"right"}.
#' @param conf.level Numeric vector or NULL. Confidence level(s) for score intervals
#'   (e.g., 0.95 or c(0.90, 0.95)). If NULL (default), confidence intervals are not computed.
#' @param bound.scores Logical. If TRUE (default), score-interval bounds are truncated to the chosen score range.
#' @param integer.scores Logical. If TRUE (default), interval bounds are rounded to integers.
#' @param score.range Numeric length-2 or NULL. If NULL (default), bounds use the \emph{observed}
#'   total-score range. If provided as \code{c(min,max)}, bounds are truncated to this range
#'   (e.g., the theoretical range).
#' @param score.display Character. How to display the \strong{Score} column:
#'   \code{"auto"} (default: numeric for "all", "lo, hi" for bands),
#'   \code{"center"} (always the band center as a number), or
#'   \code{"range"} (always "lo, hi").
#' @param digits.csem Integer. Rounding digits for \code{CSEM}. Default = 4.
#' @param digits.ci Integer. Rounding digits for CI bounds when \code{integer.scores = FALSE}. Default = 2.
#' @param quantile.type Integer in 1:9. Quantile algorithm used to build bands by quantiles.
#'   Default = 7. See \code{?quantile} for details.
#' @param smooth Logical. If TRUE, applies a polynomial smoothing to the CSEM as a function of the
#'   score (or band center). The raw Feldt estimates are kept in \code{CSEM.raw}, and the smoothed
#'   values replace \code{CSEM}. Confidence intervals are based on the smoothed CSEM. Default = FALSE.
#' @param degree Integer. Degree of the polynomial used when \code{smooth = TRUE}. Default = 2.
#' @param scale Numeric length-2 or NULL. If not NULL, defines a linear transformation of the score to a
#'   new scale with mean and SD given by \code{c(M, SD)}. For \code{score.group = "all"}, the function
#'   adds columns with scaled scores and scaled confidence limits. For banded scores, it returns a list
#'   with the original table and a separate table in the scaled metric. Default = NULL.
#' @param digits.scale Integer. Rounding digits for scaled scores and scaled CI bounds. Default = 2.
#' @param na.rm Logical. If TRUE (default), removes rows with any NA across items.
#'
#' @details
#' For each score group \(g\) with \(N_g\) examinees and \(k\) items, the conditional error variance is:
#' \deqn{\widehat{\sigma}^2_{E(cond)}(g) =
#'   \left(\frac{k}{k-1}\right)\left(\frac{N_g}{N_g-1}\right)\sum_{j=1}^k s^2_{jg}}
#' where \(s^2_{jg}\) is the sample variance of item \(j\) in group \(g\) (using divisor \(N_g-1\)).
#' Then \eqn{CSEM(g) = \sqrt{\widehat{\sigma}^2_{E(cond)}(g)}}.
#'
#' When bands are used (deciles/quartiles/quintiles/k), the representative score for centering the CI is
#' determined by \code{aggregate}: mean, median, or midpoint of the observed score range in the band.
#'
#' Confidence intervals for the \emph{true score} are centered on that representative score:
#' \deqn{[\; \text{Score}_\text{center} - z \cdot CSEM,\;\; \text{Score}_\text{center} + z \cdot CSEM \;]}
#' where \(z\) is the standard normal quantile for the requested \code{conf.level}.
#'
#' When \code{smooth = TRUE}, a polynomial regression \eqn{CSEM \sim \text{poly(Score, degree)}} is fitted
#' to the raw Feldt estimates, and the predicted values (truncated at zero) replace \code{CSEM}. The raw
#' estimates are retained in \code{CSEM.raw}.
#'
#' When \code{scale = c(M, SD)}, the observed total scores are linearly transformed to a new metric with
#' mean \code{M} and standard deviation \code{SD}. For \code{score.group = "all"}, the function returns
#' a single data frame with scaled scores and scaled confidence intervals. For banded scores, it returns a
#' list with the original table and a separate table with band centers, ranges, and confidence intervals
#' in the scaled metric.
#'
#' @return
#' If \code{scale = NULL} or \code{score.group = "all"}, a data frame with columns:
#' \itemize{
#'   \item \strong{Score}: numeric if \code{score.group = "all"} (unless \code{score.display = "range"}),
#'         otherwise either a center value or a "lo, hi" string depending on \code{score.display}.
#'   \item \strong{CSEM}: conditional standard error (numeric). If \code{smooth = TRUE}, the smoothed
#'         estimates; the raw Feldt estimates are in \code{CSEM.raw}.
#'   \item \strong{CSEM.raw}: raw Feldt CSEM estimates (only when \code{smooth = TRUE}).
#'   \item \strong{lwr.ci.xx}, \strong{upr.ci.xx}: lower/upper CI bounds for each requested level (if any).
#'   \item \strong{Score.scaled}, \strong{lwr.ci.xx.scaled}, \strong{upr.ci.xx.scaled}:
#'         scaled scores and CIs when \code{scale} is not NULL and \code{score.group = "all"}.
#'   \item \strong{n}: number of examinees in the score or band.
#' }
#'
#' If \code{scale} is not NULL and \code{score.group} is a banded option ("deciles", "quartiles",
#' "quintiles", or "k"), a list with:
#' \itemize{
#'   \item \strong{table}: the original data frame as described above (in the raw score metric).
#'   \item \strong{scaled}: a data frame with band labels, band centers and ranges in the scaled metric,
#'         and scaled confidence intervals.
#' }
#'
#' @references
#' Feldt, L. S., Steffen, M., & Gupta, N. C. (1985).
#' A comparison of five methods for estimating the standard error of measurement at specific score levels.
#' \emph{Applied Psychological Measurement, 9}(4), 351–361. https://doi.org/10.1177/014662168500900402
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(matrix(sample(1:5, 200 * 8, replace = TRUE), ncol = 8))
#'
#' # CSEM for each unique score
#' head(csemFeldt(X, score.group = "all"))
#'
#' # Bands by quartiles with 90% and 95% CIs
#' res.q <- csemFeldt(X, score.group = "quartiles", conf.level = c(0.90, 0.95))
#'
#' # CSEM by deciles, smoothed, with scaled scores (e.g., M = 5, SD = 2)
#' res.d <- csemFeldt(X, score.group = "deciles", conf.level = 0.95,
#'                    smooth = TRUE, degree = 2, scale = c(5, 2))
#'
#' @export
csemFeldt <- function(data,
                      score.group   = c("all", "deciles", "quartiles", "quintiles", "k"),
                      k             = NULL,
                      aggregate     = c("mean", "median", "midpoint"),
                      min.n         = 5,
                      merge.method  = c("nearest","left","right"),
                      conf.level    = NULL,
                      bound.scores  = TRUE,
                      integer.scores= TRUE,
                      score.range   = NULL,
                      score.display = c("auto","center","range"),
                      digits.csem   = 4,
                      digits.ci     = 2,
                      quantile.type = 7,
                      smooth        = FALSE,
                      degree        = 2,
                      scale         = NULL,
                      digits.scale  = 2,
                      na.rm         = TRUE) {

  # --- Validate args ---
  score.group    <- match.arg(score.group)
  aggregate      <- match.arg(aggregate)
  merge.method   <- match.arg(merge.method)
  score.display  <- match.arg(score.display)

  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data.frame or matrix of item responses.")
  }
  data <- as.data.frame(data)

  if (na.rm) {
    data <- stats::na.omit(data)
  } else if (anyNA(data)) {
    stop("Missing values found. Set `na.rm = TRUE` to drop incomplete rows.")
  }

  n_persons <- nrow(data)
  k_items   <- ncol(data)
  if (n_persons < 2) stop("Not enough persons after NA handling.")
  if (k_items   < 2) stop("At least two items are required.")

  total <- rowSums(data)

  # score bounds (observed or user-specified) for original metric
  if (is.null(score.range)) {
    bound_min <- min(total)
    bound_max <- max(total)
  } else {
    if (!is.numeric(score.range) || length(score.range) != 2L || score.range[1] >= score.range[2]) {
      stop("`score.range` must be numeric c(min, max) with min < max.")
    }
    bound_min <- score.range[1]
    bound_max <- score.range[2]
  }

  # mean and sd of total scores (for scaling)
  if (!is.null(scale)) {
    if (!is.numeric(scale) || length(scale) != 2L) {
      stop("`scale` must be numeric c(M, SD).")
    }
    M_scale <- scale[1]
    SD_scale <- scale[2]
    if (!is.finite(SD_scale) || SD_scale <= 0) stop("`scale[2]` (SD) must be > 0.")
    mean_total <- mean(total)
    sd_total   <- stats::sd(total)
    if (sd_total <= 0) stop("Total-score SD is zero; scaling is not defined.")
    transform_score <- function(x) {
      z <- (x - mean_total) / sd_total
      M_scale + SD_scale * z
    }
  }

  # --- Build groups (indices of persons) ---
  cut_by_quantiles <- function(x, m, qtype) {
    probs <- seq(0, 1, length.out = m + 1)
    qs <- unique(stats::quantile(x, probs = probs, type = qtype, names = FALSE))
    if (length(qs) <= 2L) {
      # fallback: all unique scores as groups
      return(split(seq_along(x), factor(x, levels = sort(unique(x)))))
    }
    f <- cut(x, breaks = qs, include.lowest = TRUE, right = TRUE)
    split(seq_along(x), f, drop = TRUE)
  }

  if (score.group == "all") {
    f <- factor(total, levels = sort(unique(total)))
    group_list <- split(seq_along(total), f)
  } else {
    m <- switch(score.group,
                deciles   = 10L,
                quartiles = 4L,
                quintiles = 5L,
                k = {
                  if (is.null(k) || !is.numeric(k) || k < 2L)
                    stop("Provide a valid `k` (>=2) when score.group = 'k'.")
                  as.integer(k)
                })
    group_list <- cut_by_quantiles(total, m, quantile.type)
  }

  # --- Merge undersized groups ---
  merge_small_groups <- function(group_list, tot_scores, min_n, method = "nearest") {
    sizes <- sapply(group_list, length)
    if (length(group_list) == 0) return(group_list)

    while (any(sizes < min_n) && length(group_list) > 1) {
      i <- which(sizes < min_n)[1]

      pick_neighbor <- function(i, method) {
        if (i == 1L) return(2L)
        if (i == length(group_list)) return(length(group_list) - 1L)
        if (method == "left")  return(i - 1L)
        if (method == "right") return(i + 1L)
        # nearest by median of total score
        med_i <- stats::median(tot_scores[group_list[[i]]])
        med_l <- stats::median(tot_scores[group_list[[i-1]]])
        med_r <- stats::median(tot_scores[group_list[[i+1]]])
        if (abs(med_i - med_l) <= abs(med_i - med_r)) i - 1L else i + 1L
      }

      j <- pick_neighbor(i, method)
      new_group <- c(group_list[[i]], group_list[[j]])
      new_index <- sort(c(i, j))

      keep_idx   <- setdiff(seq_along(group_list), new_index)
      insert_pos <- min(new_index)
      new_list <- list()
      if (insert_pos > 1) new_list <- c(new_list, group_list[keep_idx[keep_idx < insert_pos]])
      new_list <- c(new_list, list(new_group))
      if (insert_pos < length(group_list)) new_list <- c(new_list, group_list[keep_idx[keep_idx > insert_pos]])
      group_list <- new_list
      sizes <- sapply(group_list, length)
    }
    group_list
  }

  if (score.group != "all" && min.n > 1) {
    group_list <- merge_small_groups(group_list, total, min.n, merge.method)
  }

  # --- Compute raw Feldt CSEM per group ---
  build_row_raw <- function(idx) {
    Ng <- length(idx)
    if (Ng < 2) {
      s2j <- rep(0, k_items)
    } else {
      s2j <- apply(data[idx, , drop = FALSE], 2, stats::var)
      s2j[is.na(s2j)] <- 0
    }
    sigma2_cond <- (k_items/(k_items - 1)) * (Ng/(Ng - 1)) * sum(s2j)
    sigma2_cond <- max(sigma2_cond, 0)
    CSEM_raw <- sqrt(sigma2_cond)

    scores_g <- total[idx]
    if (score.group == "all") {
      score_lo <- score_hi <- as.numeric(names(table(scores_g)))[1]
      score_center <- score_lo
    } else {
      rng <- range(scores_g)
      score_lo <- rng[1]
      score_hi <- rng[2]
      score_center <- switch(aggregate,
                             mean     = mean(scores_g),
                             median   = stats::median(scores_g),
                             midpoint = mean(rng))
    }

    data.frame(
      .score_lo = score_lo,
      .score_hi = score_hi,
      .center   = score_center,
      CSEM.raw  = CSEM_raw,
      n         = Ng,
      stringsAsFactors = FALSE
    )
  }

  df <- do.call(rbind, lapply(group_list, build_row_raw))

  # order by center
  df <- df[order(df$.center), , drop = FALSE]
  rownames(df) <- NULL

  # --- Optional smoothing ---
  if (smooth) {
    if (!is.numeric(degree) || degree < 1) stop("`degree` must be a positive integer.")
    # fit polynomial on centers
    fit <- stats::lm(CSEM.raw ~ stats::poly(.center, degree, raw = TRUE), data = df)
    CSEM_sm <- stats::predict(fit, newdata = df)
    CSEM_sm[!is.finite(CSEM_sm)] <- NA_real_
    CSEM_sm <- pmax(0, CSEM_sm)
    df$CSEM <- CSEM_sm
  } else {
    df$CSEM <- df$CSEM.raw
  }

  # --- Confidence intervals in original score metric ---
  if (!is.null(conf.level)) {
    clv <- sort(unique(conf.level))
    for (cl in clv) {
      if (!is.numeric(cl) || cl <= 0 || cl >= 1) next
      z <- stats::qnorm(1 - (1 - cl)/2)
      lwr <- df$.center - z * df$CSEM
      upr <- df$.center + z * df$CSEM
      if (bound.scores) {
        lwr <- pmax(lwr, bound_min)
        upr <- pmin(upr, bound_max)
      }
      if (integer.scores) {
        lwr <- round(lwr)
        upr <- round(upr)
      } else {
        lwr <- round(lwr, digits.ci)
        upr <- round(upr, digits.ci)
      }
      tag <- sprintf("%.0f", 100 * cl)
      df[[paste0("lwr.ci.", tag)]] <- lwr
      df[[paste0("upr.ci.", tag)]] <- upr
    }
  }

  # --- Build Score label according to score.display ---
  make_range_label <- function(lo, hi) paste0(lo, ", ", hi)

  if (score.group == "all") {
    # all groups are single scores
    if (score.display == "range") {
      df$Score <- make_range_label(df$.score_lo, df$.score_hi)
    } else {
      # "auto" and "center" -> numeric center
      df$Score <- df$.center
    }
  } else {
    # banded scores
    if (score.display == "center") {
      df$Score <- df$.center
    } else {
      # "auto" and "range" -> "lo, hi"
      df$Score <- make_range_label(df$.score_lo, df$.score_hi)
    }
  }

  # --- Round CSEM and optionally drop CSEM.raw if no smoothing ---
  df$CSEM <- round(df$CSEM, digits.csem)
  if (!smooth) {
    df$CSEM.raw <- NULL
  } else {
    df$CSEM.raw <- round(df$CSEM.raw, digits.csem)
  }

  # drop internal columns
  df$.score_lo <- NULL
  df$.score_hi <- NULL
  df$.center   <- NULL

  # reorder columns: Score, CSEM(.raw), CIs, n
  ci_cols <- grep("^lwr\\.ci\\.|^upr\\.ci\\.", names(df), value = TRUE)
  base_cols <- c("Score", "CSEM")
  if (smooth) base_cols <- c("Score", "CSEM", "CSEM.raw")
  df <- df[, c(base_cols, setdiff(names(df), c(base_cols, ci_cols, "n")), ci_cols, "n"), drop = FALSE]

  # --- Scaling to new metric (if requested) ---
  if (is.null(scale)) {
    # no scaling; return data.frame
    return(df)
  }

  # we have scale = c(M, SD) and transform_score() already defined
  # we will use the center of each group for scaling
  # first, reconstruct numeric centers from Score and group structure
  # For "all", Score is center; for bands, we reconstruct using original totals

  # rebuild centers before we stripped them: we can recompute from row labels? No.
  # But we still have 'total' and 'group_list', and df rows correspond to group_list order.
  # So recompute centers from group_list and aggregate.

  centers <- sapply(group_list, function(idx) {
    scores_g <- total[idx]
    if (score.group == "all") {
      as.numeric(names(table(scores_g)))[1]
    } else {
      rng <- range(scores_g)
      switch(aggregate,
             mean     = mean(scores_g),
             median   = stats::median(scores_g),
             midpoint = mean(rng))
    }
  })
  centers <- centers[order(centers)]  # same order as df

  # scaled center and (if banded) scaled range
  center_scaled <- transform_score(centers)
  center_scaled <- round(center_scaled, digits.scale)

  # scaled CIs
  scaled_ci <- NULL
  if (!is.null(conf.level)) {
    clv <- sort(unique(conf.level))
    scaled_ci <- list()
    for (cl in clv) {
      tag <- sprintf("%.0f", 100 * cl)
      lwr_name <- paste0("lwr.ci.", tag)
      upr_name <- paste0("upr.ci.", tag)
      if (lwr_name %in% names(df) && upr_name %in% names(df)) {
        lwr_sc <- round(transform_score(df[[lwr_name]]), digits.scale)
        upr_sc <- round(transform_score(df[[upr_name]]), digits.scale)
        scaled_ci[[paste0("lwr.ci.", tag, ".scaled")]] <- lwr_sc
        scaled_ci[[paste0("upr.ci.", tag, ".scaled")]] <- upr_sc
      }
    }
    if (length(scaled_ci)) {
      scaled_ci <- as.data.frame(scaled_ci, check.names = FALSE)
    }
  }

  if (score.group == "all") {
    # add scaled columns to same df
    df$Score.scaled <- center_scaled
    if (!is.null(scaled_ci)) {
      df <- cbind(df, scaled_ci)
    }
    return(df)
  }

  # Banded case: build separate scaled table
  # band range in original metric
  band_lo <- sapply(group_list, function(idx) min(total[idx]))
  band_hi <- sapply(group_list, function(idx) max(total[idx]))
  ord_b  <- order(centers)
  band_lo <- band_lo[ord_b]
  band_hi <- band_hi[ord_b]

  band_lo_sc <- round(transform_score(band_lo), digits.scale)
  band_hi_sc <- round(transform_score(band_hi), digits.scale)

  Band <- paste0(band_lo, ", ", band_hi)
  Score.range.scaled <- paste0(band_lo_sc, ", ", band_hi_sc)

  scaled_df <- data.frame(
    Band                = Band,
    Score.center.scaled = center_scaled,
    Score.range.scaled  = Score.range.scaled,
    n                   = df$n,
    stringsAsFactors = FALSE
  )
  if (!is.null(scaled_ci) && nrow(scaled_ci) == nrow(scaled_df)) {
    scaled_df <- cbind(scaled_df, scaled_ci)
  }

  # return both
  list(
    table  = df,
    scaled = scaled_df
  )
}

# Optional wrapper for backward compatibility
#' @rdname csemFeldt
#' @export
feldt.CSEM <- function(...) {
  .Deprecated("csemFeldt")
  csemFeldt(...)
}
