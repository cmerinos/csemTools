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
#'   the smoothed CSEM for every integer score from the minimum observed to the
#'   maximum observed (or to the theoretical range if \code{score.range} provided).
#'   Default = \code{FALSE}.
#' @param ci Logical. If \code{TRUE}, compute confidence intervals for the true score.
#'   Default = \code{FALSE}.
#' @param conf.level Numeric. Confidence level for intervals (default 0.95).
#' @param digits Integer. Number of decimal places for CSEM and confidence limits.
#'   Default = 3.
#' @param score.range Optional numeric vector of length 2 (min, max) to define the
#'   range of possible scores. Only used when \code{full.range = TRUE} and
#'   \code{smooth = TRUE}. If \code{NULL}, the observed range is used.
#'
#' @return A list with elements:
#' \item{CSEM}{data.frame with columns: \code{score} (integer score),
#'   \code{n} (number of subjects, \code{NA} for predicted scores),
#'   \code{CSEM.raw} (raw estimate, present if \code{smooth = FALSE}),
#'   \code{CSEM.smooth} (smoothed estimate, present if \code{smooth = TRUE}),
#'   and \code{lwr.ci}, \code{upr.ci} if \code{ci = TRUE}.}
#' \item{binned.CSEM}{(only if \code{bin.score} is an integer) data.frame with
#'   quantile group statistics: \code{group}, \code{range}, \code{n},
#'   \code{CSEM.mean}, and \code{lwr.ci}, \code{upr.ci} if \code{ci = TRUE}.}
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
  if (!is.logical(ci))
    stop("'ci' must be logical.")
  if (ci && (conf.level <= 0 || conf.level >= 1))
    stop("'conf.level' must be between 0 and 1.")
  if (!is.numeric(digits) || digits < 0)
    stop("'digits' must be a non-negative integer.")
  if (!is.null(score.range) && (!is.numeric(score.range) || length(score.range) != 2))
    stop("'score.range' must be a numeric vector of length 2 (min, max).")

  # --- Compute total scores and differences ---
  total1 <- rowSums(half1, na.rm = TRUE)
  total2 <- rowSums(half2, na.rm = TRUE)
  total <- total1 + total2
  diff <- total1 - total2

  # --- Helper: merge small score groups (only when smooth=FALSE and bin.score=NULL) ---
  # Returns a data frame with columns: score, n, CSEM.raw
  merge_small_groups <- function(scores, diffs, min_n) {
    # Build a data frame of unique scores with their n and raw CSEM
    unique_scores <- sort(unique(scores))
    raw_list <- list()
    for (s in unique_scores) {
      idx <- which(scores == s)
      n_s <- length(idx)
      if (n_s >= 2) {
        csem_val <- sd(diffs[idx])
      } else {
        csem_val <- NA_real_
      }
      raw_list[[length(raw_list) + 1]] <- data.frame(
        score = s,
        n = n_s,
        CSEM.raw = csem_val,
        stringsAsFactors = FALSE
      )
    }
    raw_df <- do.call(rbind, raw_list)
    raw_df <- raw_df[order(raw_df$score), ]

    # If all scores have n >= min_n, return as is (but keep only those with n>=2)
    if (all(raw_df$n >= min_n | is.na(raw_df$CSEM.raw))) {
      raw_df <- raw_df[!is.na(raw_df$CSEM.raw), ]
      return(raw_df)
    }

    # We'll create intervals (blocks) of consecutive scores.
    # Start with each score as its own block.
    blocks <- list()
    for (i in 1:nrow(raw_df)) {
      blocks[[i]] <- list(
        score_min = raw_df$score[i],
        score_max = raw_df$score[i],
        indices = which(scores == raw_df$score[i]),
        n = raw_df$n[i],
        csem = raw_df$CSEM.raw[i]
      )
    }

    # While any block has n < min_n and we have more than one block, merge with nearest neighbor.
    while (length(blocks) > 1 && any(sapply(blocks, function(b) b$n) < min_n)) {
      # Find first block with n < min_n
      idx_small <- which(sapply(blocks, function(b) b$n) < min_n)[1]
      # Determine neighbor to merge with (by score proximity)
      if (idx_small == 1) {
        idx_neighbor <- 2
      } else if (idx_small == length(blocks)) {
        idx_neighbor <- length(blocks) - 1
      } else {
        left_dist <- abs(blocks[[idx_small]]$score_min - blocks[[idx_small - 1]]$score_max)
        right_dist <- abs(blocks[[idx_small + 1]]$score_min - blocks[[idx_small]]$score_max)
        idx_neighbor <- if (left_dist <= right_dist) idx_small - 1 else idx_small + 1
      }
      # Merge blocks
      new_indices <- c(blocks[[idx_small]]$indices, blocks[[idx_neighbor]]$indices)
      new_min <- min(blocks[[idx_small]]$score_min, blocks[[idx_neighbor]]$score_min)
      new_max <- max(blocks[[idx_small]]$score_max, blocks[[idx_neighbor]]$score_max)
      # Recompute CSEM using all individuals in the merged block
      new_csem <- sd(diffs[new_indices])
      new_block <- list(
        score_min = new_min,
        score_max = new_max,
        indices = new_indices,
        n = length(new_indices),
        csem = new_csem
      )
      # Replace the two blocks with the new one at the position of the smaller index
      keep <- setdiff(seq_along(blocks), c(idx_small, idx_neighbor))
      new_pos <- min(idx_small, idx_neighbor)
      new_blocks <- list()
      if (new_pos > 1) new_blocks <- c(new_blocks, blocks[keep[keep < new_pos]])
      new_blocks <- c(new_blocks, list(new_block))
      if (new_pos <= length(blocks)) new_blocks <- c(new_blocks, blocks[keep[keep >= new_pos]])
      blocks <- new_blocks
    }

    # Expand each block to all integer scores in its range, with n = NA for intermediate scores
    result_list <- list()
    for (blk in blocks) {
      score_seq <- blk$score_min:blk$score_max
      for (s in score_seq) {
        result_list[[length(result_list) + 1]] <- data.frame(
          score = s,
          n = if (blk$score_min == blk$score_max) blk$n else NA_integer_,
          CSEM.raw = blk$csem,
          stringsAsFactors = FALSE
        )
      }
    }
    result <- do.call(rbind, result_list)
    result <- result[order(result$score), ]
    rownames(result) <- NULL
    return(result)
  }

  # --- Step 1: Obtain raw CSEM estimates ---
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
    raw_df <- raw_df[!is.na(raw_df$CSEM.raw), ]  # keep only scores with valid CSEM
  } else {
    if (is.null(bin.score)) {
      # Individual scores with possible merging by min.n
      raw_df <- merge_small_groups(total, diff, min.n)
      # raw_df already has columns score, n, CSEM.raw
    } else {
      # bin.score is integer: quantile groups based on persons.
      # First compute raw CSEM per unique score (with at least 2 persons)
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

      # Create quantile groups on persons (not on unique scores)
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
        binned_list[[i]] <- data.frame(
          group = i,
          range = range_str,
          n = n_persons,
          CSEM.mean = csem_mean,
          stringsAsFactors = FALSE
        )
      }
      binned_df <- do.call(rbind, binned_list)
      raw_df <- temp_df  # individual scores with raw CSEM
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
    # Build polynomial formula
    poly_terms <- paste("I(x^", 1:degree, ")", sep = "", collapse = " + ")
    form <- as.formula(paste("y ~", poly_terms))
    fit <- lm(form, data = data.frame(x = x, y = y))
    # Predict for scores in fit_df
    pred_var <- predict(fit, newdata = data.frame(x = fit_df$score))
    pred_var <- pmax(pred_var, 0)
    csem_smooth <- sqrt(pred_var)
    smooth_df <- fit_df
    smooth_df$CSEM.smooth <- round(csem_smooth, digits)
    if (full.range) {
      if (!is.null(score.range)) {
        full_min <- score.range[1]
        full_max <- score.range[2]
      } else {
        full_min <- min(total, na.rm = TRUE)
        full_max <- max(total, na.rm = TRUE)
      }
      full_scores <- seq(full_min, full_max, by = 1)
      pred_full_var <- predict(fit, newdata = data.frame(x = full_scores))
      pred_full_var <- pmax(pred_full_var, 0)
      csem_full <- sqrt(pred_full_var)
      smooth_df <- data.frame(
        score = full_scores,
        n = NA_integer_,
        CSEM.smooth = round(csem_full, digits),
        stringsAsFactors = FALSE
      )
    }
    # Extract polynomial info
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
    poly_out <- list(
      coefficients = coef_df,
      fit.statistics = fit_stats,
      degree = degree
    )
    csem_final <- smooth_df
  } else {
    # No smoothing: csem_final is raw_df, but rename CSEM.raw to CSEM
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

    # Also add CIs to binned.CSEM if it exists and bin.score is integer
    if (!is.null(bin.score) && exists("binned_df")) {
      z <- stats::qnorm(1 - (1 - conf.level) / 2)
      lwr_bin <- binned_df$CSEM.mean - z * binned_df$CSEM.mean  # CI on the mean? Actually for the true score, we center on the mean score of the group.
      upr_bin <- binned_df$CSEM.mean + z * binned_df$CSEM.mean
      # But more logical: center on the mean total score of the group? The article is not explicit.
      # We'll center on the mean raw score of the group (we don't have it now). Simpler: use the group's CSEM.mean for both center and error.
      # I'll compute the mean total score per group and use that as center.
      # Let's compute mean total score per bin:
      bin_means <- tapply(total, groups, mean)
      bin_means <- as.numeric(bin_means[order(unique(groups))])
      lwr_bin <- bin_means - z * binned_df$CSEM.mean
      upr_bin <- bin_means + z * binned_df$CSEM.mean
      binned_df$lwr.ci <- round(lwr_bin, digits)
      binned_df$upr.ci <- round(upr_bin, digits)
    }
  }

  # --- Step 4: Prepare output list ---
  out <- list(CSEM = csem_final)
  if (!is.null(bin.score) && exists("binned_df")) {
    out$binned.CSEM <- binned_df
  }
  if (smooth && !is.null(poly_out)) {
    out$polynomial <- poly_out
  }

  class(out) <- "csemThorndike"
  return(out)
}
