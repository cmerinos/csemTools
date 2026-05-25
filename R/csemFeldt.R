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
#' @param ci Logical. If \code{TRUE}, compute confidence intervals for the true score.
#'   Default \code{FALSE}.
#' @param conf.level Numeric vector or NULL. Confidence level(s) for score intervals
#'   (e.g., 0.95 or c(0.90, 0.95)). If NULL (default), confidence intervals are not computed.
#' @param bound.scores Logical. If TRUE (default), score-interval bounds are truncated to the chosen score range.
#' @param score.range Numeric length-2 or NULL. If NULL (default), bounds use the \emph{observed}
#'   total-score range. If provided as \code{c(min,max)}, bounds are truncated to this range
#'   (e.g., the theoretical range).
#' @param score.display Character. How to display the \strong{Score} column:
#'   \code{"auto"} (default: numeric for "all", "lo, hi" for bands),
#'   \code{"center"} (always the band center as a number), or
#'   \code{"range"} (always "lo, hi").
#' @param digits.csem Integer. Rounding digits for \code{CSEM} and confidence limits. Default = 4.
#' @param smooth Logical. If TRUE, applies a polynomial smoothing to the CSEM as a function of the
#'   score (or band center). The raw Feldt estimates are kept in \code{CSEM.raw}, and the smoothed
#'   values replace \code{CSEM}. Confidence intervals are based on the smoothed CSEM. Default = FALSE.
#' @param degree Integer. Degree of the polynomial used when \code{smooth = TRUE}. Default = 2.
#' @param full.range Logical. If TRUE, CSEM values are reported for every integer score from the
#'   minimum to the maximum possible. Requires \code{smooth = TRUE} for stable estimates;
#'   if smooth=FALSE, linear interpolation is used (warning issued). Default FALSE.
#' @param summary Logical. If TRUE, returns additional components with summary information
#'   (general statistics and, if smoothing, polynomial coefficients and fit measures).
#'   Default FALSE.
#' @param na.rm Logical. If TRUE (default), removes rows with any NA across items.
#'
#' @details
#' For each score group \eqn{g} with \eqn{N_g} examinees and \eqn{k} items, the conditional error variance is:
#' \deqn{\widehat{\sigma}^2_{E(cond)}(g) =
#'   \left(\frac{k}{k-1}\right)\left(\frac{N_g}{N_g-1}\right)\sum_{j=1}^k s^2_{jg}}
#' where \eqn{s^2_{jg}} is the sample variance of item \eqn{j} in group \eqn{g} (using divisor \eqn{N_g-1}).
#' Then \eqn{CSEM(g) = \sqrt{\widehat{\sigma}^2_{E(cond)}(g)}}.
#'
#' When bands are used (deciles/quartiles/quintiles/k), the representative score for centering the CI is
#' determined by \code{aggregate}: mean, median, or midpoint of the observed score range in the band.
#'
#' Confidence intervals for the \emph{true score} are centered on that representative score:
#' \deqn{[\; \text{Score}_\text{center} - z \cdot CSEM,\;\; \text{Score}_\text{center} + z \cdot CSEM \;]}
#' where \eqn{z} is the standard normal quantile for the requested \code{conf.level}.
#'
#' When \code{smooth = TRUE}, a polynomial regression \eqn{CSEM \sim \text{poly}(Score, degree)} is fitted
#' to the raw Feldt estimates, and the predicted values (truncated at zero) replace \code{CSEM}. The raw
#' estimates are retained in \code{CSEM.raw}.
#'
#' If \code{full.range = TRUE}, the function evaluates the CSEM (smoothed or interpolated) at each integer
#' score within the observed (or user-supplied) range. This is particularly useful for creating a complete
#' CSEM table. It is strongly recommended to use \code{smooth = TRUE} together with \code{full.range = TRUE}
#' to avoid erratic estimates.
#'
#' @return
#' A list with components:
#' \item{table}{Data frame with columns:
#'   \itemize{
#'     \item \strong{Score}: numeric if \code{score.group = "all"} (unless \code{score.display = "range"}),
#'           otherwise either a center value or a "lo, hi" string depending on \code{score.display}.
#'     \item \strong{CSEM}: conditional standard error (numeric). If \code{smooth = TRUE}, the smoothed
#'           estimates; the raw Feldt estimates are in \code{CSEM.raw}.
#'     \item \strong{CSEM.raw}: raw Feldt CSEM estimates (only when \code{smooth = TRUE}).
#'     \item \strong{lwr.xx}, \strong{upr.xx}: lower/upper CI bounds for each requested level (if any).
#'     \item \strong{n}: number of examinees in the score or band (for full.range = TRUE, n is NA).
#'   }}
#' \item{summary}{(if \code{summary = TRUE}) a list with:
#'   \itemize{
#'     \item \code{general}: a data frame with method, number of groups, persons, items, score range, etc.
#'     \item \code{polynomial}: (if \code{smooth = TRUE}) a list with coefficients, R-squared, AIC, BIC,
#'           residual standard error, and other fit statistics.
#'   }}
#'
#' @references
#' Feldt, L. S., Steffen, M., & Gupta, N. C. (1985).
#' A comparison of five methods for estimating the standard error of measurement at specific score levels.
#' \emph{Applied Psychological Measurement}, 9(4), 351–361. \doi{10.1177/014662168500900402}
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(matrix(sample(1:5, 200 * 8, replace = TRUE), ncol = 8))
#'
#' # CSEM for each unique score, with 95% CI
#' res <- csemFeldt(X, score.group = "all", ci = TRUE, conf.level = 0.95)
#' head(res$table)
#'
#' # Bands by quartiles, smoothed, full range, with summary
#' res2 <- csemFeldt(X, score.group = "quartiles", smooth = TRUE, degree = 2,
#'                   full.range = TRUE, summary = TRUE)
#' res2$table
#' res2$summary$general
#' res2$summary$polynomial
#'
#' @export
csemFeldt <- function(data,
                      score.group   = c("all", "deciles", "quartiles", "quintiles", "k"),
                      k             = NULL,
                      aggregate     = c("mean", "median", "midpoint"),
                      min.n         = 5,
                      merge.method  = c("nearest", "left", "right"),
                      ci            = FALSE,
                      conf.level    = NULL,
                      bound.scores  = TRUE,
                      score.range   = NULL,
                      score.display = c("auto", "center", "range"),
                      digits.csem   = 4,
                      smooth        = FALSE,
                      degree        = 2,
                      full.range    = FALSE,
                      summary       = FALSE,
                      na.rm         = TRUE) {

  # --- Argument matching ---
  score.group   <- match.arg(score.group)
  aggregate     <- match.arg(aggregate)
  merge.method  <- match.arg(merge.method)
  score.display <- match.arg(score.display)

  # --- Data preparation ---
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
  n_items   <- ncol(data)
  if (n_persons < 2) stop("Not enough persons after NA handling.")
  if (n_items   < 2) stop("At least two items are required.")

  total <- rowSums(data)

  # Score bounds (observed or user-specified)
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

  # --- Build groups (indices of persons) ---
  cut_by_quantiles <- function(x, m) {
    probs <- seq(0, 1, length.out = m + 1)
    qs <- unique(stats::quantile(x, probs = probs, type = 7, names = FALSE))
    if (length(qs) <= 2L) {
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
    group_list <- cut_by_quantiles(total, m)
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
      s2j <- rep(0, n_items)
    } else {
      s2j <- apply(data[idx, , drop = FALSE], 2, stats::var)
      s2j[is.na(s2j)] <- 0
    }
    sigma2_cond <- (n_items/(n_items - 1)) * (Ng/(Ng - 1)) * sum(s2j)
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
  df <- df[order(df$.center), , drop = FALSE]
  rownames(df) <- NULL

  # --- Optional smoothing and polynomial info ---
  poly_info <- NULL
  if (smooth) {
    if (!is.numeric(degree) || degree < 1) stop("`degree` must be a positive integer.")
    fit <- stats::lm(CSEM.raw ~ stats::poly(.center, degree, raw = TRUE), data = df)
    CSEM_sm <- stats::predict(fit, newdata = df)
    CSEM_sm[!is.finite(CSEM_sm)] <- NA_real_
    CSEM_sm <- pmax(0, CSEM_sm)
    df$CSEM <- CSEM_sm

    coef_tab <- summary(fit)$coefficients
    rownames(coef_tab) <- c("(Intercept)", paste0("beta", 1:degree))
    poly_info <- list(
      degree = degree,
      coefficients = coef_tab,
      r.squared = summary(fit)$r.squared,
      adj.r.squared = summary(fit)$adj.r.squared,
      residual.se = summary(fit)$sigma,
      AIC = AIC(fit),
      BIC = BIC(fit),
      logLik = as.numeric(logLik(fit)),
      deviance = deviance(fit),
      df.residual = df.residual(fit)
    )
  } else {
    df$CSEM <- df$CSEM.raw
  }

  # --- Confidence intervals (original metric) ---
  if (ci && !is.null(conf.level)) {
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
      lwr <- round(lwr, digits.csem)
      upr <- round(upr, digits.csem)
      tag <- sprintf("%.0f", 100 * cl)
      df[[paste0("lwr.", tag)]] <- lwr
      df[[paste0("upr.", tag)]] <- upr
    }
  }

  # --- Build Score label according to score.display ---
  make_range_label <- function(lo, hi) paste0(lo, ", ", hi)

  if (score.group == "all") {
    if (score.display == "range") {
      df$Score <- make_range_label(df$.score_lo, df$.score_hi)
    } else {
      df$Score <- df$.center
    }
  } else {
    if (score.display == "center") {
      df$Score <- df$.center
    } else {
      df$Score <- make_range_label(df$.score_lo, df$.score_hi)
    }
  }

  # --- Round CSEM and optionally drop CSEM.raw ---
  df$CSEM <- round(df$CSEM, digits.csem)
  if (!smooth) {
    df$CSEM.raw <- NULL
  } else {
    df$CSEM.raw <- round(df$CSEM.raw, digits.csem)
  }

  # --- Full range evaluation (if requested) ---
  if (full.range) {
    if (!smooth) {
      warning("full.range = TRUE is recommended only with smooth = TRUE. ",
              "Using linear interpolation between observed centers.")
    }
    all_scores <- seq(from = ceiling(bound_min), to = floor(bound_max), by = 1)
    if (smooth) {
      newdata <- data.frame(.center = all_scores)
      pred_CSEM <- stats::predict(fit, newdata = newdata)
      pred_CSEM <- pmax(0, pred_CSEM)
      full_df <- data.frame(
        Score = all_scores,
        CSEM = round(pred_CSEM, digits.csem),
        n = NA
      )
      if (smooth) full_df$CSEM.raw <- NA
      if (ci && !is.null(conf.level)) {
        for (cl in clv) {
          z <- stats::qnorm(1 - (1 - cl)/2)
          lwr <- all_scores - z * pred_CSEM
          upr <- all_scores + z * pred_CSEM
          if (bound.scores) {
            lwr <- pmax(lwr, bound_min)
            upr <- pmin(upr, bound_max)
          }
          lwr <- round(lwr, digits.csem)
          upr <- round(upr, digits.csem)
          tag <- sprintf("%.0f", 100 * cl)
          full_df[[paste0("lwr.", tag)]] <- lwr
          full_df[[paste0("upr.", tag)]] <- upr
        }
      }
      df <- full_df
    } else {
      interp_csem <- stats::approx(x = df$.center, y = df$CSEM, xout = all_scores, rule = 2)
      full_df <- data.frame(
        Score = all_scores,
        CSEM = round(interp_csem$y, digits.csem),
        n = NA
      )
      if (ci && !is.null(conf.level)) {
        for (cl in clv) {
          z <- stats::qnorm(1 - (1 - cl)/2)
          lwr <- all_scores - z * interp_csem$y
          upr <- all_scores + z * interp_csem$y
          if (bound.scores) {
            lwr <- pmax(lwr, bound_min)
            upr <- pmin(upr, bound_max)
          }
          lwr <- round(lwr, digits.csem)
          upr <- round(upr, digits.csem)
          tag <- sprintf("%.0f", 100 * cl)
          full_df[[paste0("lwr.", tag)]] <- lwr
          full_df[[paste0("upr.", tag)]] <- upr
        }
      }
      df <- full_df
    }
  } else {
    # Remove internal columns from original df
    df$.score_lo <- NULL
    df$.score_hi <- NULL
    df$.center   <- NULL
    # Reorder columns: Score, CSEM, CSEM.raw (if exists), CIs, n
    ci_cols <- grep("^lwr\\.|^upr\\.", names(df), value = TRUE)
    base_cols <- c("Score", "CSEM")
    if (smooth) base_cols <- c("Score", "CSEM", "CSEM.raw")
    df <- df[, c(base_cols, setdiff(names(df), c(base_cols, ci_cols, "n")), ci_cols, "n"), drop = FALSE]
  }

  # --- Summary information (if requested) ---
  summary_out <- NULL
  if (summary) {
    general_info <- data.frame(
      parameter = c("score.group", "n_groups", "n_persons", "n_items",
                    "min_score", "max_score", "smooth", "full.range"),
      value = c(score.group, if (full.range) NA else length(group_list), n_persons, n_items,
                bound_min, bound_max, smooth, full.range),
      stringsAsFactors = FALSE
    )
    summary_out <- list(general = general_info)
    if (smooth && !is.null(poly_info)) summary_out$polynomial <- poly_info
  }

  # --- Output ---
  out <- list(table = df)
  if (summary) out$summary <- summary_out
  return(out)
}
