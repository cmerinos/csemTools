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
#'   the smoothed CSEM for every integer score from 0 to the total number of items.
#'   Default = \code{FALSE}.
#' @param ci Logical. If \code{TRUE}, compute confidence intervals for the true score.
#'   Default = \code{FALSE}.
#' @param conf.level Numeric. Confidence level for intervals (default 0.95).
#' @param digits Integer. Number of decimal places for CSEM and confidence limits.
#'   Default = 3.
#'
#' @return A list with elements:
#' \item{CSEM}{data.frame with columns: \code{score} (integer score),
#'   \code{n} (number of subjects, \code{NA} for predicted scores),
#'   \code{CSEM.raw} (raw estimate, present if \code{smooth = FALSE}),
#'   \code{CSEM.smooth} (smoothed estimate, present if \code{smooth = TRUE}),
#'   and \code{lwr.ci}, \code{upr.ci} if \code{ci = TRUE}.}
#' \item{binned.CSEM}{(only if \code{bin.score} is an integer) data.frame with
#'   quantile group statistics: \code{group}, \code{range}, \code{n}, \code{CSEM.mean}.}
#' \item{polynomial}{(only if \code{smooth = TRUE}) list with polynomial coefficients,
#'   fit statistics, and degree.}
#'
#' @examples
#' \dontrun{
#' # Basic usage with individual scores (no smoothing)
#' res <- csemThorndike(half1, half2, bin.score = NULL, smooth = FALSE)
#'
#' # With quantile grouping and smoothing
#' res2 <- csemThorndike(half1, half2, bin.score = 4, smooth = TRUE, degree = 2)
#' }
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
                          digits = 3) {

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

  # --- Compute total scores and differences ---
  total1 <- rowSums(half1, na.rm = TRUE)
  total2 <- rowSums(half2, na.rm = TRUE)
  total  <- total1 + total2
  diff   <- total1 - total2

  # --------------------------------------------------------------------------
  # Helper: merge small score groups (only used when smooth = FALSE and bin.score = NULL)
  # Returns a data frame with columns: score, n, CSEM.raw
  # For merged ranges, each integer score in the range appears as a separate row
  # with n = NA and CSEM.raw equal to the group's CSEM.
  # --------------------------------------------------------------------------
  merge_small_groups <- function(scores, diffs, min_n) {
    # First, compute raw CSEM for each unique score (with at least 2 persons)
    unique_scores <- sort(unique(scores))
    raw_list <- list()
    for (s in unique_scores) {
      idx <- which(scores == s)
      n_s <- length(idx)
      if (n_s >= 2) {
        csem <- sd(diffs[idx])
      } else {
        csem <- NA_real_
      }
      raw_list[[length(raw_list) + 1]] <- data.frame(
        score = s,
        n = n_s,
        CSEM.raw = csem,
        stringsAsFactors = FALSE
      )
    }
    raw_df <- do.call(rbind, raw_list)
    raw_df <- raw_df[order(raw_df$score), ]

    # Separate scores with valid CSEM (n>=2) and those with insufficient n
    valid <- !is.na(raw_df$CSEM.raw)
    if (all(valid)) {
      # No merging needed
      return(raw_df)
    }

    # We will create blocks of consecutive scores where:
    # - any score with insufficient n is merged with its nearest valid neighbor(s)
    # To simplify: walk from low to high, accumulate a block until total n >= min_n,
    # then compute CSEM for the block using all individuals in that block,
    # and expand to each integer score in the block.
    # This ensures each block has at least min_n persons.
    # If a block contains only scores with n<2, we still compute CSEM if total n>=2.

    # Build a table of all individuals with their scores and diffs
    df_all <- data.frame(score = scores, diff = diffs)
    df_all <- df_all[order(df_all$score), ]

    blocks <- list()
    current_block <- data.frame()
    current_min <- NULL
    current_max <- NULL

    i <- 1
    while (i <= nrow(df_all)) {
      if (nrow(current_block) == 0) {
        # start a new block at current row
        current_block <- df_all[i, , drop = FALSE]
        current_min <- current_max <- df_all$score[i]
        i <- i + 1
      } else {
        # check if adding the next row would cause the block to have too few persons?
        # Actually we want to accumulate until the block size >= min_n OR we reach a score
        # that is not consecutive? The requirement is only about sample size, not about
        # score contiguity. However, merging non-consecutive scores would be weird.
        # We'll keep merging consecutive scores until block size >= min_n.
        # But if the next score is more than 1 away, we might want to keep separate.
        # Simpler: use the approach of merging only the scores that are too small with their
        # nearest neighbor (by score). This is more complex.
        #
        # After reconsideration, the original intention was to merge only those score levels
        # that have less than min_n individuals, merging them with the closest score level
        # (by score value). I will implement that simpler logic.

        # I'll rewrite the merging as follows:
        # 1. Identify all score levels with n < min_n and n >= 2? Actually if n<2, CSEM is NA.
        #    We need to merge them with a neighboring level that has at least min_n persons or at least 2.
        #    We'll keep a list of groups (intervals) and iteratively merge the smallest group
        #    with its nearest neighbor.
      }
    }

    # Alternative robust algorithm:
    # Start with each unique score as its own group.
    # While any group has total n < min_n and more than one group exists:
    #   Find the group with smallest n (and n < min_n)
    #   Merge it with the adjacent group (by score) that is closest in score value.
    #   Update the group's score range and recompute CSEM using all individuals in the merged group.
    # After merging, expand each group to all integer scores within its range.

    # Build initial groups
    groups <- list()
    for (s in unique_scores) {
      idx <- which(scores == s)
      groups[[length(groups) + 1]] <- list(
        scores = s,           # original unique score
        idx = idx,
        n = length(idx),
        min_score = s,
        max_score = s
      )
    }
    # Assign names for easier reference (the original score)
    names(groups) <- as.character(unique_scores)

    # Function to compute CSEM for a group (using its individuals)
    compute_group_csem <- function(grp) {
      if (length(grp$idx) < 2) return(NA_real_)
      sd(diffs[grp$idx])
    }

    # Merge while any group has n < min_n and there is more than one group
    while (any(sapply(groups, function(g) g$n) < min_n) && length(groups) > 1) {
      # Find first group with n < min_n (by smallest score)
      idx_small <- which(sapply(groups, function(g) g$n) < min_n)[1]
      grp_small <- groups[[idx_small]]
      # Determine neighbor to merge with (nearest by score)
      if (idx_small == 1) {
        idx_neighbor <- 2
      } else if (idx_small == length(groups)) {
        idx_neighbor <- length(groups) - 1
      } else {
        left_score <- groups[[idx_small - 1]]$scores  # representative score
        right_score <- groups[[idx_small + 1]]$scores
        current_score <- grp_small$scores
        if (abs(current_score - left_score) <= abs(current_score - right_score)) {
          idx_neighbor <- idx_small - 1
        } else {
          idx_neighbor <- idx_small + 1
        }
      }
      grp_neighbor <- groups[[idx_neighbor]]
      # Merge the two
      new_idx <- c(grp_small$idx, grp_neighbor$idx)
      new_min <- min(grp_small$min_score, grp_neighbor$min_score)
      new_max <- max(grp_small$max_score, grp_neighbor$max_score)
      # The representative score can be the midpoint or mean; we'll keep min and max.
      new_group <- list(
        scores = NULL,   # not a single score anymore
        idx = new_idx,
        n = length(new_idx),
        min_score = new_min,
        max_score = new_max
      )
      # Replace the two groups with the merged one at the position of the smaller index
      if (idx_small > idx_neighbor) {
        keep <- setdiff(seq_along(groups), c(idx_small, idx_neighbor))
        new_pos <- idx_neighbor
      } else {
        keep <- setdiff(seq_along(groups), c(idx_small, idx_neighbor))
        new_pos <- idx_small
      }
      new_groups <- list()
      if (new_pos > 1) new_groups <- c(new_groups, groups[keep[keep < new_pos]])
      new_groups <- c(new_groups, list(new_group))
      if (new_pos <= length(groups)) new_groups <- c(new_groups, groups[keep[keep >= new_pos]])
      groups <- new_groups
      # Rename? Not needed.
    }

    # Now compute CSEM for each final group and expand to integer scores
    result_list <- list()
    for (grp in groups) {
      csem_val <- compute_group_csem(grp)
      if (is.na(csem_val)) next
      score_seq <- grp$min_score:grp$max_score
      for (s in score_seq) {
        result_list[[length(result_list) + 1]] <- data.frame(
          score = s,
          n = if (grp$min_score == grp$max_score) grp$n else NA_integer_,
          CSEM.raw = csem_val,
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
    # For smoothing, we need raw estimates at each unique score with at least 2 persons.
    unique_scores <- sort(unique(total))
    raw_list <- list()
    for (s in unique_scores) {
      idx <- which(total == s)
      if (length(idx) >= 2) {
        csem_raw <- sd(diff[idx])
      } else {
        csem_raw <- NA_real_
      }
      raw_list[[length(raw_list) + 1]] <- data.frame(
        score = s,
        n = length(idx),
        CSEM.raw = csem_raw,
        stringsAsFactors = FALSE
      )
    }
    raw_df <- do.call(rbind, raw_list)
    raw_df <- raw_df[!is.na(raw_df$CSEM.raw), ]  # remove scores with <2 persons
  } else {
    if (is.null(bin.score)) {
      raw_df <- merge_small_groups(total, diff, min.n)
    } else {
      # bin.score is integer: quantile groups based on persons, not unique scores
      # First compute raw CSEM per unique score
      unique_scores <- sort(unique(total))
      temp_list <- list()
      for (s in unique_scores) {
        idx <- which(total == s)
        if (length(idx) >= 2) {
          csem_raw <- sd(diff[idx])
        } else {
          csem_raw <- NA_real_
        }
        temp_list[[length(temp_list) + 1]] <- data.frame(
          score = s,
          n = length(idx),
          CSEM.raw = csem_raw,
          stringsAsFactors = FALSE
        )
      }
      temp_df <- do.call(rbind, temp_list)
      temp_df <- temp_df[!is.na(temp_df$CSEM.raw), ]

      # Create quantile groups on persons
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
    # Use raw_df (score and CSEM.raw)
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
      full_scores <- 0:n_items_total
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
    lwr <- pmax(lwr, 0)
    upr <- pmin(upr, n_items_total)
    csem_final$lwr.ci <- round(lwr, digits)
    csem_final$upr.ci <- round(upr, digits)
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
