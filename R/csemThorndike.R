#' Thorndike's method for Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Computes the Conditional Standard Error of Measurement (CSEM) using Thorndike's
#' half-test difference method. The CSEM for a given total score is the standard
#' deviation of the difference between two parallel half-tests.
#'
#' @param half1 data.frame or matrix with items from the first half of the test.
#' @param half2 data.frame or matrix with items from the second half of the test.
#'   Must have the same number of rows as \code{half1}.
#' @param bin.score Either \code{NULL} (default) to work with individual scores,
#'   or an integer \code{k} to divide the total scores into \code{k} quantile groups.
#' @param min.n Integer. Minimum number of observations per score group when
#'   \code{bin.score = NULL} and \code{smooth = FALSE}. Groups with fewer than
#'   \code{min.n} observations are merged with the nearest neighbor. Ignored
#'   when \code{smooth = TRUE} or \code{bin.score} is an integer. Default = 20.
#' @param smooth Logical. If \code{TRUE}, applies polynomial smoothing to the
#'   squared CSEM estimates. Default = \code{FALSE}.
#' @param degree Integer. Degree of the polynomial used when \code{smooth = TRUE}.
#'   Default = 2.
#' @param full.range Logical. If \code{TRUE} and \code{smooth = TRUE}, evaluates
#'   the smoothed CSEM for every integer score from \code{score.range[1]} to
#'   \code{score.range[2]}. Requires \code{score.range} to be provided.
#'   Default = \code{FALSE}.
#' @param ci Logical. If \code{TRUE}, compute confidence intervals for the true score.
#'   Default = \code{FALSE}.
#' @param conf.level Numeric. Confidence level for intervals (default 0.95).
#' @param digits Integer. Number of decimal places for CSEM and confidence limits.
#'   Default = 3.
#' @param score.range Optional numeric vector of length 2 (min, max) defining the
#'   theoretical score range. Required when \code{full.range = TRUE}. If \code{NULL}
#'   but \code{full.range = TRUE}, the observed range is used (with a warning).
#'
#' @return A list with elements:
#' \item{CSEM}{data.frame with columns: \code{score} (integer score),
#'   \code{n} (number of subjects, 0 for unobserved, NA for merged scores),
#'   \code{CSEM.raw} (raw estimate, present if \code{smooth = FALSE}),
#'   \code{CSEM.smooth} (smoothed estimate, present if \code{smooth = TRUE}),
#'   and \code{lwr.ci}, \code{upr.ci} if \code{ci = TRUE}.}
#' \item{binned.CSEM}{(only if \code{bin.score} is an integer) data.frame with
#'   quantile group statistics: \code{group}, \code{range}, \code{n},
#'   \code{mean_score} (mean total score in the group),
#'   \code{CSEM.mean} (mean CSEM in the group, raw or smoothed),
#'   and \code{lwr.ci}, \code{upr.ci} if \code{ci = TRUE}.}
#' \item{polynomial}{(only if \code{smooth = TRUE}) list with polynomial coefficients,
#'   fit statistics, and degree.}
#'
#' @export
csemThorndike <- function(half1, half2,
                          bin.score = NULL,
                          min.n = 20,
                          smooth = FALSE,
                          degree = 2,
                          full.range = FALSE,
                          ci = FALSE,
                          conf.level = 0.95,
                          digits = 3,
                          score.range = NULL) {

  # --- Argument validation ---
  if (missing(half1) || missing(half2))
    stop("Both 'half1' and 'half2' must be provided.")
  half1 <- as.data.frame(half1)
  half2 <- as.data.frame(half2)
  if (nrow(half1) != nrow(half2))
    stop("'half1' and 'half2' must have the same number of rows.")
  n_items_total <- ncol(half1) + ncol(half2)

  if (!is.null(bin.score) && (!is.numeric(bin.score) || length(bin.score) != 1 || bin.score < 2))
    stop("If 'bin.score' is provided, it must be an integer >= 2.")
  if (!is.numeric(min.n) || min.n < 1)
    stop("'min.n' must be a positive integer.")
  if (!is.logical(smooth) || length(smooth) != 1)
    stop("'smooth' must be a single logical value.")
  if (smooth && (!is.numeric(degree) || degree < 1))
    stop("'degree' must be a positive integer when smooth = TRUE.")
  if (!is.logical(full.range))
    stop("'full.range' must be logical.")
  if (full.range && !smooth)
    warning("full.range = TRUE is recommended only with smooth = TRUE.")
  if (full.range && is.null(score.range)) {
    warning("score.range not provided. Using observed score range for full.range.")
    score.range <- c(min(total, na.rm = TRUE), max(total, na.rm = TRUE))
  }
  if (!is.null(score.range) && (!is.numeric(score.range) || length(score.range) != 2))
    stop("'score.range' must be a numeric vector of length 2 (min, max).")
  if (!is.logical(ci))
    stop("'ci' must be logical.")
  if (ci && (conf.level <= 0 || conf.level >= 1))
    stop("'conf.level' must be between 0 and 1.")
  if (!is.numeric(digits) || digits < 0)
    stop("'digits' must be a non-negative integer.")

  # --- Compute total scores and differences ---
  total1 <- rowSums(half1, na.rm = TRUE)
  total2 <- rowSums(half2, na.rm = TRUE)
  total <- total1 + total2
  diff <- total1 - total2

  # --- Helper: merge small groups (only when smooth=FALSE and bin.score=NULL) ---
  # This version respects scores with n >= min.n (they stay individual)
  merge_small_groups <- function(scores, diffs, min_n) {
    unique_scores <- sort(unique(scores))
    tbl <- data.frame(score = unique_scores,
                      n = as.integer(table(scores)[as.character(unique_scores)]),
                      stringsAsFactors = FALSE)
    # Raw CSEM for each unique score (only if n>=2)
    tbl$CSEM.raw <- sapply(tbl$score, function(s) {
      idx <- which(scores == s)
      if (length(idx) >= 2) sd(diffs[idx]) else NA_real_
    })
    # Keep scores with n >= min_n AND valid CSEM
    keep <- which(tbl$n >= min_n & !is.na(tbl$CSEM.raw))
    result <- tbl[keep, , drop = FALSE]
    # Scores with 0 < n < min_n (or n>=2 but still < min_n) need merging
    to_merge <- which(tbl$n > 0 & !(1:nrow(tbl) %in% keep))
    if (length(to_merge) > 0) {
      merge_df <- tbl[to_merge, , drop = FALSE]
      # Order by score (they are already ordered)
      # Build contiguous blocks
      blocks <- list()
      current_block <- data.frame()
      for (i in 1:nrow(merge_df)) {
        if (nrow(current_block) == 0) {
          current_block <- merge_df[i, , drop = FALSE]
        } else {
          if (merge_df$score[i] == current_block$score[nrow(current_block)] + 1) {
            current_block <- rbind(current_block, merge_df[i, , drop = FALSE])
          } else {
            blocks <- c(blocks, list(current_block))
            current_block <- merge_df[i, , drop = FALSE]
          }
        }
      }
      if (nrow(current_block) > 0) blocks <- c(blocks, list(current_block))
      # For each block, compute combined CSEM and expand to each score
      for (blk in blocks) {
        scores_blk <- blk$score
        idx_all <- unlist(lapply(scores_blk, function(s) which(scores == s)))
        csem_block <- sd(diffs[idx_all])
        for (s in scores_blk) {
          result <- rbind(result, data.frame(score = s,
                                             n = NA_integer_,
                                             CSEM.raw = csem_block,
                                             stringsAsFactors = FALSE))
        }
      }
    }
    # Add scores with n=0 (not observed) to complete the observed range
    all_observed <- min(scores):max(scores)
    missing <- setdiff(all_observed, result$score)
    if (length(missing) > 0) {
      missing_df <- data.frame(score = missing, n = 0L, CSEM.raw = NA_real_, stringsAsFactors = FALSE)
      result <- rbind(result, missing_df)
    }
    result <- result[order(result$score), ]
    rownames(result) <- NULL
    return(result)
  }

  # --- Step 1: raw CSEM estimates ---
  if (smooth) {
    # For smoothing, compute raw CSEM at each unique score with at least 2 persons.
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
      raw_list[[length(raw_list) + 1]] <- data.frame(
        score = s,
        n = n_s,
        CSEM.raw = csem_raw,
        stringsAsFactors = FALSE
      )
    }
    raw_df <- do.call(rbind, raw_list)
    raw_df <- raw_df[!is.na(raw_df$CSEM.raw), ]
  } else {
    if (is.null(bin.score)) {
      raw_df <- merge_small_groups(total, diff, min.n)
    } else {
      # bin.score integer: quantile groups based on persons.
      unique_scores <- sort(unique(total))
      temp_list <- list()
      for (s in unique_scores) {
        idx <- which(total == s)
        n_s <- length(idx)
        if (n_s >= 2) {
          csem_raw <- sd(diff[idx])
        } else {
          csem_raw <- NA_real_
        }
        temp_list[[length(temp_list) + 1]] <- data.frame(
          score = s,
          n = n_s,
          CSEM.raw = csem_raw,
          stringsAsFactors = FALSE
        )
      }
      temp_df <- do.call(rbind, temp_list)
      temp_df <- temp_df[!is.na(temp_df$CSEM.raw), ]
      # Quantile groups on persons
      q <- stats::quantile(total, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
      q <- unique(q)
      groups <- cut(total, breaks = q, include.lowest = TRUE, right = TRUE)
      group_levels <- levels(groups)
      binned_list <- list()
      for (i in seq_along(group_levels)) {
        idx_in_group <- which(groups == group_levels[i])
        scores_in_group <- unique(total[idx_in_group])
        sub_df <- temp_df[temp_df$score %in% scores_in_group, ]
        if (nrow(sub_df) == 0) next
        csem_mean <- mean(sub_df$CSEM.raw)
        range_str <- paste0(min(scores_in_group), "-", max(scores_in_group))
        n_persons <- length(idx_in_group)
        mean_score <- mean(total[idx_in_group])
        binned_list[[i]] <- data.frame(
          group = i,
          range = range_str,
          n = n_persons,
          mean_score = mean_score,
          CSEM.mean = csem_mean,
          stringsAsFactors = FALSE
        )
      }
      binned_df <- do.call(rbind, binned_list)
      raw_df <- temp_df
    }
  }

  # --- Step 2: Smoothing (if requested) ---
  poly_out <- NULL
  csem_final <- NULL
  if (smooth) {
    fit_df <- raw_df[is.finite(raw_df$CSEM.raw) & !is.na(raw_df$CSEM.raw), ]
    if (nrow(fit_df) < degree + 1)
      stop("Not enough data points to fit a polynomial of degree ", degree, ". Reduce degree or set smooth = FALSE.")
    y <- fit_df$CSEM.raw^2
    x <- fit_df$score
    poly_terms <- paste("I(x^", 1:degree, ")", sep = "", collapse = " + ")
    form <- as.formula(paste("y ~", poly_terms))
    fit <- lm(form, data = data.frame(x = x, y = y))
    pred_var <- predict(fit, newdata = data.frame(x = fit_df$score))
    pred_var <- pmax(pred_var, 0)
    csem_smooth <- sqrt(pred_var)
    smooth_df <- fit_df
    smooth_df$CSEM.smooth <- round(csem_smooth, digits)
    if (full.range) {
      full_min <- score.range[1]
      full_max <- score.range[2]
      full_scores <- seq(full_min, full_max, by = 1)
      pred_full_var <- predict(fit, newdata = data.frame(x = full_scores))
      pred_full_var <- pmax(pred_full_var, 0)
      csem_full <- sqrt(pred_full_var)
      # Build full data frame with n from actual frequencies (0 for unobserved)
      full_n <- sapply(full_scores, function(s) sum(total == s))
      smooth_df <- data.frame(
        score = full_scores,
        n = full_n,
        CSEM.smooth = round(csem_full, digits),
        stringsAsFactors = FALSE
      )
    }
    # Polynomial info
    coef_sum <- summary(fit)$coefficients
    coef_df <- data.frame(
      term = rownames(coef_sum),
      estimate = coef_sum[, "Estimate"],
      std.error = coef_sum[, "Std. Error"],
      t.value = coef_sum[, "t value"],
      p.value = coef_sum[, "Pr(>|t|)"],
      row.names = NULL
    )
    fit_stats <- list(
      r.squared = summary(fit)$r.squared,
      adj.r.squared = summary(fit)$adj.r.squared,
      AIC = AIC(fit),
      BIC = BIC(fit),
      logLik = as.numeric(logLik(fit)),
      residual.se = summary(fit)$sigma,
      df.residual = df.residual(fit)
    )
    poly_out <- list(coefficients = coef_df, fit.statistics = fit_stats, degree = degree)
    csem_final <- smooth_df

    # If bin.score is integer and smooth=TRUE, build binned.CSEM using predicted values
    if (!is.null(bin.score)) {
      q <- stats::quantile(total, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
      q <- unique(q)
      groups <- cut(total, breaks = q, include.lowest = TRUE, right = TRUE)
      group_levels <- levels(groups)
      binned_list2 <- list()
      for (i in seq_along(group_levels)) {
        idx_in_group <- which(groups == group_levels[i])
        scores_in_group <- unique(total[idx_in_group])
        # Predict CSEM.smooth for each score in the group
        pred_var_grp <- predict(fit, newdata = data.frame(x = scores_in_group))
        pred_var_grp <- pmax(pred_var_grp, 0)
        csem_pred <- sqrt(pred_var_grp)
        csem_mean_smooth <- mean(csem_pred, na.rm = TRUE)
        range_str <- paste0(min(scores_in_group), "-", max(scores_in_group))
        n_persons <- length(idx_in_group)
        mean_score <- mean(total[idx_in_group])
        binned_list2[[i]] <- data.frame(
          group = i,
          range = range_str,
          n = n_persons,
          mean_score = mean_score,
          CSEM.mean = csem_mean_smooth,
          stringsAsFactors = FALSE
        )
      }
      binned_df <- do.call(rbind, binned_list2)
    }
  } else {
    # No smoothing
    csem_final <- raw_df
    names(csem_final)[names(csem_final) == "CSEM.raw"] <- "CSEM"
    csem_final$n <- as.integer(csem_final$n)
  }

  # --- Step 3: Confidence intervals (if requested) ---
  if (ci) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    if (smooth) {
      csem_vals <- csem_final$CSEM.smooth
    } else {
      csem_vals <- csem_final$CSEM
    }
    lwr <- csem_final$score - z * csem_vals
    upr <- csem_final$score + z * csem_vals
    csem_final$lwr.ci <- round(lwr, digits)
    csem_final$upr.ci <- round(upr, digits)

    # CIs for binned.CSEM if exists
    if (exists("binned_df") && !is.null(binned_df) && nrow(binned_df) > 0) {
      lwr_bin <- binned_df$mean_score - z * binned_df$CSEM.mean
      upr_bin <- binned_df$mean_score + z * binned_df$CSEM.mean
      binned_df$lwr.ci <- round(lwr_bin, digits)
      binned_df$upr.ci <- round(upr_bin, digits)
    }
  }

  # --- Output ---
  out <- list(CSEM = csem_final)
  if (exists("binned_df") && !is.null(binned_df) && nrow(binned_df) > 0) {
    out$binned.CSEM <- binned_df
  }
  if (smooth && !is.null(poly_out)) {
    out$polynomial <- poly_out
  }
  class(out) <- "csemThorndike"
  return(out)
}
