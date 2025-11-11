#' Feldt's Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Calculates the Conditional Standard Error of Measurement (CSEM) following the
#' Feldt, Steffen, & Gupta (1985) variance-components approach. The CSEM is
#' estimated at each score level or score band using item variances within groups.
#'
#' @param data A data frame or matrix with item responses (subjects in rows, items in columns).
#'             Items can be dichotomous or polytomous (Likert). All items contribute to the total score.
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
#' @param conf.level Numeric vector or NULL. Confidence level(s) for score intervals (e.g., 0.95 or c(0.90, 0.95)).
#'   If NULL (default), confidence intervals are not computed.
#' @param bound.scores Logical. If TRUE (default), score-interval bounds are truncated to the chosen score range.
#' @param integer.scores Logical. If TRUE (default), interval bounds are rounded to integers.
#' @param score.range Numeric length-2 or NULL. If NULL (default), bounds use the \emph{observed} total-score range.
#'   If provided as \code{c(min,max)}, bounds are truncated to this range (e.g., the theoretical range).
#' @param score.display Character. How to display the \strong{Score} column:
#'   \code{"auto"} (default: numeric for "all", "lo, hi" for bands),
#'   \code{"center"} (always the band center as a number), or
#'   \code{"range"} (always "lo, hi").
#' @param digits.csem Integer. Rounding digits for \code{CSEM}. Default = 4.
#' @param digits.ci Integer. Rounding digits for CI bounds when \code{integer.scores = FALSE}. Default = 2.
#' @param quantile.type Integer in 1:9. Quantile algorithm used to build bands by quantiles.
#'   Default = 7. See \code{?quantile} for details.
#' @param na.rm Logical. If TRUE (default), removes rows with any NA across items.
#'
#' @details
#' For each score group \(g\) with \(N_g\) examinees and \(k\) items, the conditional error variance is:
#' \deqn{\widehat{\sigma}^2_{E(cond)}(g) =
#'   \left(\frac{k}{k-1}\right)\left(\frac{N_g}{N_g-1}\right)\sum_{j=1}^k s^2_{jg}}
#' where \(s^2_{jg}\) is the sample variance of item \(j\) in group \(g\) (using divisor \(N_g-1\)).
#' Then \eqn{CSEM(g) = \sqrt{\widehat{\sigma}^2_{E(cond)}(g)}}.
#'
#' When bands are used (deciles/quartiles/quintiles/k), the representative score for centering the CI is determined by
#' \code{aggregate}: mean, median, or midpoint of the observed score range in the band.
#'
#' Confidence intervals for the \emph{true score} are centered on that representative score:
#' \deqn{[\; \text{Score}_\text{center} - z \cdot CSEM,\;\; \text{Score}_\text{center} + z \cdot CSEM \;]}
#' where \(z\) is the standard normal quantile for the requested \code{conf.level}.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \strong{Score}: numeric if \code{score.group = "all"} (unless \code{score.display = "range"}),
#'         otherwise either a center value or a "lo, hi" string depending on \code{score.display}.
#'   \item \strong{CSEM}: conditional standard error (numeric).
#'   \item \strong{lwr.ci.xx}, \strong{upr.ci.xx}: lower/upper CI bounds for each requested level (if any).
#'   \item \strong{n}: number of examinees in the score or band.
#' }
#'
#' @references
#' Feldt, L. S., Steffen, M., & Gupta, N. C. (1985).
#' A comparison of five methods for estimating the standard error of measurement at specific score levels.
#' \emph{Applied Psychological Measurement, 9}(4), 351–361. https://doi.org/10.1177/014662168500900402
#'
#' @examples
#' set.seed(123)
#'
#' #Articifial data
#' exampleData <- data.frame(matrix(sample(1:5, 200 * 8, replace = TRUE), ncol = 8))
#'
#' #Run: every scale single point
#' csemFeldt(exampleData, score.group = "all")
#'
#' #score un quartiles, two confidence levles
#' csemFeldt(exampleData, score.group = "quartiles", conf.level = c(0.90, 0.95))
#'
#' @export
csemFeldt <- function(data,
                       score.group = c("all", "deciles", "quartiles", "quintiles", "k"),
                       k = NULL,
                       aggregate = c("mean", "median", "midpoint"),
                       min.n = 5,
                       merge.method = c("nearest","left","right"),
                       conf.level = NULL,
                       bound.scores = TRUE,
                       integer.scores = TRUE,
                       score.range = NULL,
                       score.display = c("auto","center","range"),
                       digits.csem = 4,
                       digits.ci = 2,
                       quantile.type = 7,
                       na.rm = TRUE) {

  # --- Validate args ---
  score.group  <- match.arg(score.group)
  aggregate    <- match.arg(aggregate)
  merge.method <- match.arg(merge.method)
  score.display<- match.arg(score.display)

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

  # score bounds (observed or user-specified)
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

  # --- Build groups ---
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
                  if (is.null(k) || !is.numeric(k) || k < 2L) stop("Provide a valid `k` (>=2) when score.group = 'k'.")
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

      keep_idx <- setdiff(seq_along(group_list), new_index)
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

  # --- Compute CSEM per group ---
  build_row <- function(idx) {
    Ng <- length(idx)
    if (Ng < 2) {
      s2j <- rep(0, k_items)
    } else {
      s2j <- apply(data[idx, , drop = FALSE], 2, stats::var)
      s2j[is.na(s2j)] <- 0
    }
    sigma2_cond <- (k_items/(k_items - 1)) * (Ng/(Ng - 1)) * sum(s2j)
    sigma2_cond <- max(sigma2_cond, 0)
    CSEM <- sqrt(sigma2_cond)

    scores_g <- total[idx]
    if (score.group == "all") {
      score_lo <- score_hi <- as.numeric(names(table(scores_g)))[1]
      score_center <- score_lo
      score_label_auto <- score_center
      score_label_range <- paste0(score_lo, ", ", score_hi)
    } else {
      rng <- range(scores_g)
      score_lo <- rng[1]; score_hi <- rng[2]
      score_center <- switch(aggregate,
                             mean    = mean(scores_g),
                             median  = stats::median(scores_g),
                             midpoint= mean(rng))
      score_label_auto  <- paste0(score_lo, ", ", score_hi)  # for bands
      score_label_range <- score_label_auto
    }

    # select display
    Score <- switch(score.display,
                    auto   = if (score.group == "all") score_label_auto else score_label_auto,
                    center = score_center,
                    range  = score_label_range)

    # build CI columns
    ci_cols <- list()
    if (!is.null(conf.level)) {
      clv <- sort(unique(conf.level))
      for (cl in clv) {
        if (!is.numeric(cl) || cl <= 0 || cl >= 1) next
        z <- stats::qnorm(1 - (1 - cl)/2)
        lwr <- score_center - z * CSEM
        upr <- score_center + z * CSEM
        if (bound.scores) {
          lwr <- max(lwr, bound_min)
          upr <- min(upr, bound_max)
        }
        if (integer.scores) {
          lwr <- round(lwr)
          upr <- round(upr)
        } else {
          lwr <- round(lwr, digits.ci)
          upr <- round(upr, digits.ci)
        }
        tag <- sprintf("%.0f", 100 * cl)
        ci_cols[[paste0("lwr.ci.", tag)]] <- lwr
        ci_cols[[paste0("upr.ci.", tag)]] <- upr
      }
    }

    # assemble row
    c(list(Score = Score, CSEM = CSEM, n = Ng), ci_cols,
      list(.center_numeric = score_center))
  }

  rows <- lapply(group_list, build_row)
  df <- do.call(rbind, lapply(rows, function(x) as.data.frame(x, check.names = FALSE)))

  # types & order
  df$CSEM <- as.numeric(df$CSEM)
  df$n    <- as.integer(df$n)
  df$.center_numeric <- as.numeric(df$.center_numeric)

  # order by the center
  ord <- order(df$.center_numeric)
  df <- df[ord, , drop = FALSE]
  rownames(df) <- NULL

  # finalize Score type
  if (score.display %in% c("center")) {
    df$Score <- as.numeric(df$Score)
  } else if (score.display == "auto" && score.group == "all") {
    # auto + all -> numeric
    df$Score <- as.numeric(df$Score)
  } else {
    df$Score <- as.character(df$Score)
  }

  # rounding for CSEM
  df$CSEM <- round(df$CSEM, digits.csem)

  # drop helper
  df$.center_numeric <- NULL

  df
}
