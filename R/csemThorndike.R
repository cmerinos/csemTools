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

  # --- Helper: merge small groups (only when smooth = FALSE and bin.score = NULL) ---
  merge_small_groups <- function(scores, diff_scores, min_n) {
    # scores: vector of total scores (integer)
    # diff_scores: vector of differences
    # Returns a list with components:
    #   groups: list of indices for each final group
    #   group_scores: representative score (the integer value if group is a single score,
    #                 otherwise a sequence from min to max)
    #   group_n: number of persons in group
    #   group_csem_raw: CSEM (sd(diff) within group)
    df <- data.frame(score = scores, diff = diff_scores)
    # Start with groups by unique score
    score_vals <- sort(unique(scores))
    groups <- lapply(score_vals, function(s) which(scores == s))
    names(groups) <- as.character(score_vals)

    # Merge while any group has size < min_n and more than one group exists
    while (any(sapply(groups, length) < min_n) && length(groups) > 1) {
      # Find first group with size < min_n
      idx_small <- which(sapply(groups, length) < min_n)[1]
      # Determine neighbor to merge with (nearest in score)
      group_scores <- as.numeric(names(groups))
      if (idx_small == 1) {
        idx_neighbor <- 2
      } else if (idx_small == length(groups)) {
        idx_neighbor <- length(groups) - 1
      } else {
        left_score <- group_scores[idx_small - 1]
        right_score <- group_scores[idx_small + 1]
        current_score <- group_scores[idx_small]
        if (abs(current_score - left_score) <= abs(current_score - right_score))
          idx_neighbor <- idx_small - 1
        else
          idx_neighbor <- idx_small + 1
      }
      # Merge the two groups
      new_indices <- c(groups[[idx_small]], groups[[idx_neighbor]])
      new_name <- paste0(min(group_scores[c(idx_small, idx_neighbor)]), "-",
                         max(group_scores[c(idx_small, idx_neighbor)]))
      # Remove the two old groups and insert the merged one at the lower position
      new_order <- sort(c(idx_small, idx_neighbor))
      groups <- groups[-new_order]
      # Insert new group at the position of the first index
      if (new_order[1] == 1) {
        groups <- c(list(new_indices), groups)
      } else if (new_order[1] == length(groups) + 1) {
        groups <- c(groups, list(new_indices))
      } else {
        groups <- c(groups[1:(new_order[1]-1)], list(new_indices), groups[new_order[1]:length(groups)])
      }
      names(groups)[new_order[1]] <- new_name
    }

    # Build output data frame for each final group
    result_list <- list()
    for (i in seq_along(groups)) {
      idx <- groups[[i]]
      grp_name <- names(groups)[i]
      n_grp <- length(idx)
      csem_raw <- sd(diff_scores[idx])
      # If group name contains a dash, it's a merged range
      if (grepl("-", grp_name)) {
        bounds <- as.numeric(strsplit(grp_name, "-")[[1]])
        score_seq <- bounds[1]:bounds[2]
        for (s in score_seq) {
          result_list[[length(result_list)+1]] <- data.frame(
            score = s,
            n = NA_integer_,
            CSEM.raw = csem_raw,
            stringsAsFactors = FALSE
          )
        }
      } else {
        s <- as.numeric(grp_name)
        result_list[[length(result_list)+1]] <- data.frame(
          score = s,
          n = n_grp,
          CSEM.raw = csem_raw,
          stringsAsFactors = FALSE
        )
      }
    }
    result <- do.call(rbind, result_list)
    result <- result[order(result$score), ]
    rownames(result) <- NULL
    return(result)
  }

  # --- Step 1: Obtain raw CSEM estimates (per unique score, possibly merged) ---
  if (smooth) {
    # For smoothing, we need raw estimates at each unique score (no merging)
    score_vals <- sort(unique(total))
    raw_df <- data.frame(score = score_vals,
                         n = as.integer(table(total)[as.character(score_vals)]),
                         CSEM.raw = sapply(score_vals, function(s) {
                           idx <- which(total == s)
                           if (length(idx) < 2) return(NA_real_)
                           sd(diff[idx])
                         }),
                         stringsAsFactors = FALSE)
    # Remove any rows with NA CSEM (less than 2 persons)
    raw_df <- raw_df[!is.na(raw_df$CSEM.raw), ]
  } else {
    if (is.null(bin.score)) {
      # Individual scores with possible merging by min.n
      raw_df <- merge_small_groups(total, diff, min.n)
    } else {
      # bin.score is an integer: quantile groups (no merging by min.n)
      # We first compute raw CSEM per unique score, then average within quantiles
      score_vals <- sort(unique(total))
      temp_df <- data.frame(score = score_vals,
                            n = as.integer(table(total)[as.character(score_vals)]),
                            CSEM.raw = sapply(score_vals, function(s) {
                              idx <- which(total == s)
                              if (length(idx) < 2) return(NA_real_)
                              sd(diff[idx])
                            }),
                            stringsAsFactors = FALSE)
      temp_df <- temp_df[!is.na(temp_df$CSEM.raw), ]
      # Create quantile groups based on individual scores (not persons)
      # We need to assign each UNIQUE SCORE to a quantile group based on its frequency?
      # Actually, bin.score = k means we want k groups of approximately equal number of PERSONS.
      # Therefore, we should assign each PERSON to a group, then average CSEM.raw per group.
      # Simpler: create quantile groups on the vector of total scores (with repetitions)
      q <- stats::quantile(total, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
      # Ensure unique cut points
      q <- unique(q)
      groups <- cut(total, breaks = q, include.lowest = TRUE, right = TRUE)
      group_levels <- levels(groups)
      # For each group, collect the unique scores within that group and average their CSEM.raw
      # (weighted by the number of persons at that score, or unweighted? The article is not explicit.
      # We'll use unweighted mean of CSEM.raw over the distinct scores in the group.)
      binned_list <- list()
      for (i in seq_along(group_levels)) {
        idx_in_group <- which(groups == group_levels[i])
        scores_in_group <- unique(total[idx_in_group])
        sub_df <- temp_df[temp_df$score %in% scores_in_group, ]
        if (nrow(sub_df) == 0) next
        csem_mean <- mean(sub_df$CSEM.raw)
        range_str <- paste0(min(scores_in_group), "-", max(scores_in_group))
        n_persons <- length(idx_in_group)
        binned_list[[i]] <- data.frame(group = i,
                                       range = range_str,
                                       n = n_persons,
                                       CSEM.mean = csem_mean,
                                       stringsAsFactors = FALSE)
        # Also, for each score in this group, we will later assign the group mean as CSEM.raw?
        # According to our earlier agreement, when bin.score is an integer, we still return
        # $CSEM with individual scores (raw) and additionally $binned.CSEM.
        # So we keep raw_df as temp_df (individual scores).
      }
      binned_df <- do.call(rbind, binned_list)
      raw_df <- temp_df  # individual scores, with raw CSEM
    }
  }

  # --- Step 2: Smoothing (if requested) ---
  poly_out <- NULL
  csem_final <- NULL
  if (smooth) {
    # Use raw_df (which contains score and CSEM.raw) for smoothing.
    # Remove any NA or infinite values
    fit_df <- raw_df[is.finite(raw_df$CSEM.raw) & !is.na(raw_df$CSEM.raw), ]
    if (nrow(fit_df) < degree + 1)
      stop("Not enough data points to fit a polynomial of degree ", degree, ". Reduce degree or set smooth = FALSE.")
    # Fit polynomial on CSEM.raw^2
    y <- fit_df$CSEM.raw^2
    x <- fit_df$score
    # Use raw polynomial (not orthogonal) for interpretability
    poly_form <- as.formula(paste("y ~", paste("I(x^", 1:degree, ")", collapse = " + ")))
    fit <- lm(poly_form)
    # Predicted values for the same scores
    pred_var <- predict(fit, newdata = data.frame(x = fit_df$score))
    pred_var <- pmax(pred_var, 0)
    csem_smooth <- sqrt(pred_var)
    # Build smoothed data frame (merge with raw_df)
    smooth_df <- fit_df
    smooth_df$CSEM.smooth <- round(csem_smooth, digits)
    if (full.range) {
      # Predict for all integer scores from 0 to n_items_total
      full_scores <- 0:n_items_total
      pred_full_var <- predict(fit, newdata = data.frame(x = full_scores))
      pred_full_var <- pmax(pred_full_var, 0)
      csem_full <- sqrt(pred_full_var)
      full_df <- data.frame(score = full_scores,
                            n = NA_integer_,
                            CSEM.smooth = round(csem_full, digits),
                            stringsAsFactors = FALSE)
      # For scores that were observed, we could keep the original n? But it's cleaner to have n=NA for all when full.range=TRUE.
      # To avoid duplication, we'll replace smooth_df with full_df.
      smooth_df <- full_df
    }
    # Extract polynomial information
    coef_sum <- summary(fit)$coefficients
    coef_df <- data.frame(term = rownames(coef_sum),
                          estimate = coef_sum[, "Estimate"],
                          std.error = coef_sum[, "Std. Error"],
                          t.value = coef_sum[, "t value"],
                          p.value = coef_sum[, "Pr(>|t|)"],
                          row.names = NULL)
    fit_stats <- list(
      r.squared = summary(fit)$r.squared,
      adj.r.squared = summary(fit)$adj.r.squared,
      AIC = AIC(fit),
      BIC = BIC(fit),
      logLik = as.numeric(logLik(fit)),
      residual.se = summary(fit)$sigma,
      df.residual = df.residual(fit)
    )
    poly_out <- list(coefficients = coef_df,
                     fit.statistics = fit_stats,
                     degree = degree)
    # Final CSEM output: use smooth_df (with CSEM.smooth)
    csem_final <- smooth_df
    # If raw_df is still needed? We'll keep the raw values in a separate element? Not required by spec.
    # But we can keep raw_df as separate if user wants? For simplicity, we output only smoothed.
  } else {
    # No smoothing: csem_final is raw_df (with CSEM.raw)
    csem_final <- raw_df
    names(csem_final)[names(csem_final) == "CSEM.raw"] <- "CSEM"
    # Ensure n is integer
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
    # Truncate to possible score range (0 to n_items_total)
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
