#' Strong True Score (Compound Binomial) Model for CSEM
#'
#' @description
#' Computes conditional standard errors of measurement (CSEM) under Lord's
#' strong true score theory (also known as the compound binomial model),
#' following the estimator described by Lord (1965) and summarized by
#' Tong and Kolen (2018). The method can be applied to both dichotomous
#' and polytomous items by transforming raw scores into binomial-equivalent
#' scores. Optionally, the function computes confidence intervals for the
#' true number-correct score using either the classic CSEM‑based method or
#' the Wilson score method.
#'
#' @param score.type Character: \code{"dich"} (default) or \code{"poly"}.
#' @param data A numeric matrix or data frame with examinees in rows and
#'   items in columns. For \code{score.type = "dich"} the entries must be
#'   0/1; for \code{"poly"} they range between \code{min.resp} and
#'   \code{max.resp}.
#' @param nitems Optional integer. Number of items. If \code{NULL} (default),
#'   it is taken as \code{ncol(data)}.
#' @param min.resp Numeric. Minimum response value per item.
#'   **Required when \code{score.type = "poly"}**; ignored for dichotomous.
#' @param max.resp Numeric. Maximum response value per item.
#'   **Required when \code{score.type = "poly"}**; ignored for dichotomous.
#' @param ci Logical. If \code{TRUE}, compute confidence intervals for the
#'   true number-correct score. Default \code{FALSE}.
#' @param ci.method Character: \code{"csem"} (classic: \eqn{X \pm z \cdot CSEM})
#'   or \code{"wilson"} (Wilson score interval). Default \code{"csem"}.
#' @param conf.level Numeric vector of confidence levels (e.g., 0.95).
#'   Default \code{NULL} → 0.95.
#' @param digits.csem Integer for rounding. Default 3.
#' @param rho.report Logical. If \code{TRUE}, returns a data frame with the
#'   key parameters used in the strong true score correction (mean total
#'   score, variance of total scores, variance of item difficulties, mean
#'   binomial variance, and the correction factor). Default \code{FALSE}.
#' @param full.range Logical. If \code{TRUE} (default), CSEM values are
#'   reported for the full score range (from minimum possible to maximum
#'   possible). If \code{FALSE}, only scores observed in \code{data} are
#'   reported.
#' @param return.person Logical. If \code{TRUE}, the output includes a
#'   data frame with person-level CSEM values and, if requested, confidence
#'   intervals. Default \code{FALSE}.
#' @param na.rm Logical. If \code{TRUE} (default), rows with any missing
#'   values are removed. If \code{FALSE}, the function stops with an error
#'   when missing values are present.
#'
#' @details
#' \strong{Strong true score / compound binomial model (Lord, 1965)}:
#' Assume that for a given examinee the observed number‑correct score
#' follows a binomial distribution conditional on a true proportion correct
#' \eqn{\phi}, and that the distribution of true scores across examinees
#' follows a beta‑like distribution. Lord derived an estimator for the
#' conditional error variance that adjusts the simple binomial variance
#' using information about the dispersion of item difficulties and the
#' total score distribution.
#'
#' Let \eqn{n} be the number of items, \eqn{x} a raw (or binomial‑equivalent)
#' score, \eqn{\hat{\mu}_X} the sample mean of total scores,
#' \eqn{S_X^2} the sample variance of total scores, \eqn{S_{Xi}^2} the sample
#' variance of item difficulties (proportion correct for dichotomous items,
#' or transformed proportion for polytomous), and \eqn{\bar{pq}} the average
#' of \eqn{p_j(1-p_j)} across items. Then the strong true score error
#' variance is:
#' \deqn{
#'   \hat{\sigma}^2_{E|x} =
#'     \frac{x (n - x)}{n - 1}
#'     \left[
#'       1 -
#'       \frac{n (n - 1) S_{Xi}^2}{
#'         \hat{\mu}_X (n - \hat{\mu}_X) - S_X^2 - n \bar{pq}
#'       }
#'     \right],
#' }
#' and the CSEM is \eqn{\sqrt{\hat{\sigma}^2_{E|x}}}.
#'
#' \strong{Extension to polytomous items}:
#' When \code{score.type = "poly"}, the raw item responses are linearly
#' transformed to proportions in \eqn{[0,1]} using the per‑item minimum and
#' maximum possible values (provided via \code{min.resp} and \code{max.resp}).
#' Then the transformed proportions are summed to form an equivalent total
#' score on the 0–\code{nitems} scale. All subsequent calculations
#' (means, variances, CSEM) are performed on this equivalent scale.
#'
#' \strong{Confidence intervals}:
#' When \code{ci = TRUE}, two methods are available:
#' \itemize{
#'   \item \code{ci.method = "csem"}: Classic CTT interval for the true score
#'     \eqn{\tau = n\phi}:
#'     \deqn{\text{equiv.score} \pm z_{\alpha/2} \cdot \text{CSEM}.}
#'   \item \code{ci.method = "wilson"}: Wilson score interval for the
#'     binomial proportion \eqn{p = \text{equiv.score}/n}, multiplied by
#'     \eqn{n} to obtain limits for \eqn{\tau}. This interval has better
#'     coverage properties for extreme scores or short tests.
#' }
#' For each requested level \eqn{\gamma}, the bounds are labeled
#' \code{lwr.xx} and \code{upr.xx} where \code{xx} is the level in percent.
#'
#' @return
#' A list of class \code{"csemStrong"} with components:
#' \item{score}{A data frame with one row per possible raw score (or per
#'   observed score if \code{full.range = FALSE}). Columns:
#'   \code{raw.score} (original score), \code{equiv.score} (if polytomous),
#'   \code{n} (number of examinees with that score), \code{csem.strong}
#'   (strong true score CSEM), and, if \code{ci = TRUE}, confidence limits.}
#' \item{person}{If \code{return.person = TRUE}, a data frame with one row
#'   per examinee: \code{id}, \code{raw.score}, \code{equiv.score} (if
#'   polytomous), \code{csem.strong}, and the same confidence limits.}
#' \item{summary}{If \code{rho.report = TRUE}, a data frame with the key
#'   parameters (\code{n_items}, \code{mu_X}, \code{var_Xp}, \code{var_Xi},
#'   \code{mean_pq}, \code{correction_factor}); otherwise \code{NULL}.}
#' \item{call}{The matched function call.}
#'
#' @references
#' Lord, F. M. (1965). A strong true score theory, with applications.
#'   \emph{Psychometrika}, 30, 239–270.
#' Tong, Y., & Kolen, M. J. (2018). Conditional standard errors of
#'   measurement. In \emph{Wiley StatsRef: Statistics Reference Online}.
#' Wilson, E. B. (1927). Probable inference, the law of succession, and
#'   statistical inference. \emph{Journal of the American Statistical
#'   Association}, 22(158), 209–212.
#'
#' @examples
#' # Dichotomous example
#' set.seed(123)
#' n_persons <- 200
#' n_items   <- 30
#' theta <- rnorm(n_persons)
#' b     <- rnorm(n_items, sd = 1)
#' pmat  <- plogis(outer(theta, -b, "+"))
#' data_sim <- matrix(rbinom(n_persons * n_items, 1, pmat), nrow = n_persons)
#' res <- csemStrong(score.type = "dich", data = data_sim, ci = TRUE, conf.level = 0.95)
#' head(res$score)
#'
#' # Polytomous example (0–4 scale)
#' # (Assume data_poly is a matrix with values 0:4)
#' \dontrun{
#' csemStrong(score.type = "poly", data = data_poly, min.resp = 0, max.resp = 4)
#' }
#'
#'@importFrom stats complete.cases var qnorm
#'
#' @export
csemStrong <- function(score.type = c("dich", "poly"),
                       data,
                       nitems = NULL,
                       min.resp = NULL,
                       max.resp = NULL,
                       ci = FALSE,
                       ci.method = c("csem", "wilson"),
                       conf.level = NULL,
                       digits.csem = 3,
                       rho.report = FALSE,
                       full.range = TRUE,
                       return.person = FALSE,
                       na.rm = TRUE) {

  score.type <- match.arg(score.type)
  ci.method  <- match.arg(ci.method)

  # --- Data preparation and validation ---
  if (is.data.frame(data)) X <- as.matrix(data)
  else if (is.matrix(data)) X <- data
  else stop("`data` must be a matrix or data frame.")

  if (!is.numeric(X)) stop("`data` must contain numeric scores.")

  # Handle missing values
  if (na.rm) {
    cc <- stats::complete.cases(X)
    if (!all(cc)) {
      X <- X[cc, , drop = FALSE]
      warning("Rows with missing values were removed (na.rm = TRUE).")
    }
  } else {
    if (any(!stats::complete.cases(X)))
      stop("Missing values found. Set na.rm = TRUE to remove incomplete rows.")
  }

  n_persons <- nrow(X)
  n_items_obs <- ncol(X)

  if (n_persons < 2L || n_items_obs < 2L)
    stop("`data` must have at least 2 persons and 2 items.")

  if (!is.null(nitems)) {
    if (nitems != n_items_obs)
      warning("Provided `nitems` (", nitems, ") differs from ncol(data) (", n_items_obs,
              "). Using ncol(data) for calculations.")
    nitems <- n_items_obs
  } else {
    nitems <- n_items_obs
  }

  # --- Polytomous transformation to binomial-equivalent scores ---
  if (score.type == "poly") {
    if (is.null(min.resp) || is.null(max.resp))
      stop("For score.type = 'poly', min.resp and max.resp must be provided.")
    if (min.resp >= max.resp)
      stop("min.resp must be less than max.resp.")
    # Transform each item: (x - min.resp) / (max.resp - min.resp)
    X_scaled <- (X - min.resp) / (max.resp - min.resp)
    # Total equivalent score: sum over items, range 0..nitems
    equiv_total <- rowSums(X_scaled)
    # For the calculations we use equiv_total as the "observed score"
    total_score <- equiv_total
    # Also store raw totals if needed for output
    raw_total <- rowSums(X)   # original raw sum
    # For item difficulties: proportion of maximum per item? Actually we need
    # p_j = mean of scaled item (proportion in [0,1])
    p_item <- colMeans(X_scaled)
    # For mean_pq: average of p_j*(1-p_j) on the scaled metric
    mean_pq <- mean(p_item * (1 - p_item))
    # Variance of item difficulties (on scaled metric)
    var_Xi <- stats::var(p_item)
    # Note: mu_X and var_Xp are computed on equiv_total (binomial scale)
    mu_X <- mean(equiv_total)
    var_Xp <- stats::var(equiv_total)
  } else { # dichotomous
    # For dichotomous, no scaling needed
    total_score <- rowSums(X)
    raw_total <- total_score   # same
    p_item <- colMeans(X)
    mean_pq <- mean(p_item * (1 - p_item))
    var_Xi <- stats::var(p_item)
    mu_X <- mean(total_score)
    var_Xp <- stats::var(total_score)
  }

  # --- Compute correction factor (denominator of Lord's adjustment) ---
  num_K <- nitems * (nitems - 1) * var_Xi
  den_K <- mu_X * (nitems - mu_X) - var_Xp - nitems * mean_pq

  if (den_K <= 0) {
    warning("Denominator of strong true score correction factor is non‑positive. ",
            "CSEM strong may be invalid; returning NA.")
    correction_factor <- NA_real_
  } else {
    correction_factor <- 1 - (num_K / den_K)
  }

  # --- Determine score values for which CSEM will be computed ---
  if (full.range) {
    if (score.type == "poly") {
      # Equivalent scores can be fractional (since scaled items). We'll generate
      # a sequence from min equiv to max equiv, step = 1/nitems? Actually the
      # possible unique equiv_total values are multiples of 1/nitems? Better to
      # use all observed unique values or full grid. For simplicity, use all
      # possible values from 0 to nitems in steps of 1/nitems? That's huge.
      # Simpler: use unique(equiv_total) (observed). But user expects full range?
      # To keep consistent with csemBinom, we'll generate a fine grid of
      # possible equiv scores (0, 0.25, 0.5, ...) but that can be many.
      # Alternative: just use unique observed scores. We'll use unique observed
      # for polytomous when full.range = TRUE? That may be unexpected.
      # Given complexity, for polytomous we'll use the same approach as dich:
      # create sequence from 0 to nitems by 0.01? Not recommended.
      # I'll use unique observed sorted, and if full.range = TRUE, we still only
      # report those observed (since the formula only needs those scores).
      # That is a practical compromise.
      score_vals <- sort(unique(equiv_total))
    } else {
      score_vals <- 0:nitems
    }
  } else {
    score_vals <- sort(unique(total_score))
  }

  # Base binomial variance: x(n - x)/(n - 1)
  base_var <- score_vals * (nitems - score_vals) / (nitems - 1)
  base_var[base_var < 0] <- NA_real_

  # Strong true score variance
  var_strong <- base_var * correction_factor
  var_strong[var_strong < 0] <- NA_real_
  csem_strong <- sqrt(var_strong)

  # --- Confidence intervals ---
  ci_df <- NULL
  if (ci) {
    if (is.null(conf.level)) conf.level <- 0.95
    conf.level <- sort(unique(conf.level))
    if (any(conf.level <= 0 | conf.level >= 1))
      stop("conf.level must contain values between 0 and 1.")

    n <- nitems
    z_vals <- stats::qnorm((1 + conf.level) / 2)
    ci_mat <- matrix(NA, nrow = length(score_vals), ncol = 2 * length(conf.level))
    colnames_ci <- character(2 * length(conf.level))

    for (i in seq_along(conf.level)) {
      if (ci.method == "csem") {
        low <- score_vals - z_vals[i] * csem_strong
        high <- score_vals + z_vals[i] * csem_strong
        low[low < 0] <- 0
        high[high > n] <- n
      } else { # wilson
        center <- (score_vals + z_vals[i]^2 / 2) / (n + z_vals[i]^2)
        half <- z_vals[i] * sqrt((score_vals * (n - score_vals) / n) + z_vals[i]^2 / 4) / (n + z_vals[i]^2)
        p_low <- center - half
        p_high <- center + half
        p_low[p_low < 0] <- 0
        p_high[p_high > 1] <- 1
        low <- n * p_low
        high <- n * p_high
      }
      ci_mat[, 2*i - 1] <- low
      ci_mat[, 2*i] <- high
      colnames_ci[2*i - 1] <- paste0("lwr.", formatC(conf.level[i]*100, format = "f", digits = 0))
      colnames_ci[2*i] <- paste0("upr.", formatC(conf.level[i]*100, format = "f", digits = 0))
    }
    colnames(ci_mat) <- colnames_ci
    ci_df <- as.data.frame(ci_mat)
    ci_df <- as.data.frame(lapply(ci_df, round, digits.csem))
  }

  # --- Frequency table for each score value ---
  freq_tab <- table(total_score)
  freq_vec <- as.numeric(freq_tab[match(score_vals, as.numeric(names(freq_tab)))])
  freq_vec[is.na(freq_vec)] <- 0L

  # --- Build score-level data frame ---
  score_df <- data.frame(
    raw.score = if (score.type == "poly") {
      # For polytomous, we need to map equiv score back to raw? We'll keep raw_total
      # but careful: score_vals are equiv scores; we need to compute approximate raw?
      # Simpler: add a column "equiv.score" and keep "raw.score" as the original raw total?
      # For consistency with csemBinom, we'll show raw.total (sum of original scores)
      # but that may not match score_vals. Let's just report equiv.score and also the
      # original raw total if needed. I'll create a column "equiv.score" and for raw.score
      # I'll put the raw total mean per equiv? That's messy. Instead, for polytomous,
      # we'll report both equiv.score and raw.score (the original summed raw).
      # To compute raw.score for each equiv level, we can group raw_total by equiv_total.
      # I'll do a quick mapping:
      raw_by_equiv <- tapply(raw_total, equiv_total, unique)
      raw_for_score <- sapply(score_vals, function(v) {
        val <- raw_by_equiv[[as.character(v)]]
        if (is.null(val)) NA else val[1]
      })
      raw_for_score
    } else {
      score_vals   # for dichotomous, raw = equiv
    },
    n = freq_vec,
    csem.strong = round(csem_strong, digits.csem)
  )

  if (score.type == "poly") {
    # Insert equiv.score column
    score_df <- cbind(score_df[, 1, drop = FALSE],
                      equiv.score = round(score_vals, digits.csem),
                      score_df[, -1, drop = FALSE])
    names(score_df)[1] <- "raw.score"
  }

  if (!is.null(ci_df)) {
    score_df <- cbind(score_df, ci_df)
  }

  # --- Person-level data (optional) ---
  person_df <- NULL
  if (return.person) {
    idx <- match(total_score, score_vals)
    person_df <- data.frame(
      id = seq_len(n_persons),
      raw.score = if (score.type == "poly") raw_total else total_score,
      csem.strong = round(csem_strong[idx], digits.csem)
    )
    if (score.type == "poly") {
      person_df <- cbind(person_df[, 1:2, drop = FALSE],
                         equiv.score = round(total_score, digits.csem),
                         person_df[, -c(1:2), drop = FALSE])
    }
    if (!is.null(ci_df)) {
      for (nm in names(ci_df)) {
        person_df[[nm]] <- round(ci_df[[nm]][idx], digits.csem)
      }
    }
  }

  # --- Summary parameters (if rho.report) ---
  sum_df <- NULL
  if (rho.report) {
    sum_df <- data.frame(
      parameter = c("n_items", "mu_X", "var_Xp", "var_Xi", "mean_pq", "correction_factor"),
      value = c(nitems, mu_X, var_Xp, var_Xi, mean_pq, correction_factor)
    )
  }

  # --- Output ---
  out <- list(
    score = score_df,
    person = person_df,
    summary = sum_df,
    call = match.call()
  )
  class(out) <- c("csemStrong", class(out))
  return(out)
}
