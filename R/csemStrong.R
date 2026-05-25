#' Strong True Score (Compound Binomial) Model for CSEM
#'
#' @description
#' Computes conditional standard errors of measurement (CSEM) under Lord's
#' strong true score theory (compound binomial model) for dichotomously
#' scored items. The estimator follows Lord (1965) as summarized by
#' Tong and Kolen (2018). Optionally, the function computes confidence
#' intervals for the true number‑correct score.
#'
#' @param data A numeric matrix or data frame with examinees in rows and
#'   items in columns. For \code{score.type = "dich"} entries must be 0/1;
#'   for \code{"poly"} they range between \code{min.resp} and \code{max.resp}.
#' @param score.type Character: \code{"dich"} (default) or \code{"poly"}.
#' @param nitems Optional integer. Number of items. If \code{NULL} (default),
#'   taken as \code{ncol(data)}.
#' @param min.resp Numeric. Minimum response value per item.
#'   **Required when \code{score.type = "poly"}**.
#' @param max.resp Numeric. Maximum response value per item.
#'   **Required when \code{score.type = "poly"}**.
#' @param ci Logical. If \code{TRUE}, compute confidence intervals. Default \code{FALSE}.
#' @param ci.method Character: \code{"csem"} (classic) or \code{"wilson"}.
#'   Default \code{"csem"}.
#' @param conf.level Numeric vector of confidence levels (e.g., 0.95).
#'   Default \code{NULL} → 0.95.
#' @param digits.csem Integer for rounding. Default 3.
#' @param summary Logical. If \code{TRUE}, returns a data frame with key
#'   parameters (mu_X, var_Xp, var_Xi, mean_pq, correction_factor).
#'   Default \code{FALSE}.
#' @param full.range Logical. If \code{TRUE} (default), CSEM values are
#'   reported for the full score range (all possible scores). For dichotomous,
#'   that is 0:nitems; for polytomous, all observed equivalent scores plus
#'   the extremes 0 and nitems. If \code{FALSE}, only scores observed in
#'   \code{data} are reported.
#' @param return.person Logical. If \code{TRUE}, includes a person‑level
#'   data frame. Default \code{FALSE}.
#' @param na.rm Logical. If \code{TRUE} (default), rows with missing values
#'   are removed. If \code{FALSE}, stops on missing values.
#'
#' @details
#' \strong{Strong true score model (Lord, 1965)}:
#' Let \eqn{n} be the number of items, \eqn{x} a raw (or binomial‑equivalent)
#' score, \eqn{\hat{\mu}_X} the sample mean of total scores,
#' \eqn{S_X^2} the sample variance of total scores, \eqn{S_{Xi}^2} the sample
#' variance of item difficulties (proportion correct for dichotomous,
#' scaled proportion for polytomous), and \eqn{\bar{pq}} the average of
#' \eqn{p_j(1-p_j)}. The strong true score error variance is:
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
#' and the CSEM is its square root.
#'
#' \strong{Polytomous items}:
#' Raw scores are linearly transformed to proportions in \eqn{[0,1]} using
#' \code{min.resp} and \code{max.resp}, then summed to obtain an equivalent
#' total score on the 0–\code{nitems} scale. All subsequent calculations
#' are performed on this equivalent scale.
#'
#' @return A list with components:
#' \item{score}{Data frame with columns: \code{raw.score} (original total),
#'   \code{equiv.score} (if polytomous), \code{n} (frequency), \code{csem.strong},
#'   and confidence limits if requested.}
#' \item{person}{(if \code{return.person = TRUE}) Person‑level data frame.}
#' \item{summary}{(if \code{summary = TRUE}) Data frame with parameters.}
#'
#' @references
#' Lord, F. M. (1965). A strong true score theory, with applications.
#'   \emph{Psychometrika}, 30, 239–270.
#' Tong, Y., & Kolen, M. J. (2018). Conditional standard errors of measurement.
#'   In \emph{Wiley StatsRef: Statistics Reference Online}.
#'
#' @examples
#' \dontrun{
#' # Dichotomous example
#' csemStrong(data = data_u, score.type = "dich", ci = TRUE, summary = TRUE)
#'
#' # Polytomous example
#' csemStrong(data = data_poly, score.type = "poly", min.resp = 1, max.resp = 6)
#' }
#'
#'@importFrom stats complete.cases var qnorm
#'
#' @export
csemStrong <- function(data,
                       score.type = c("dich", "poly"),
                       nitems = NULL,
                       min.resp = NULL,
                       max.resp = NULL,
                       ci = FALSE,
                       ci.method = c("csem", "wilson"),
                       conf.level = NULL,
                       digits.csem = 3,
                       summary = FALSE,
                       full.range = TRUE,
                       return.person = FALSE,
                       na.rm = TRUE) {

  score.type <- match.arg(score.type)
  ci.method  <- match.arg(ci.method)

  # --- Data preparation ---
  if (is.data.frame(data)) X <- as.matrix(data)
  else if (is.matrix(data)) X <- data
  else stop("`data` must be a matrix or data frame.")

  if (!is.numeric(X)) stop("`data` must contain numeric scores.")

  # Missing values
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
  if (n_persons < 2 || n_items_obs < 2)
    stop("`data` must have at least 2 persons and 2 items.")

  if (is.null(nitems)) nitems <- n_items_obs
  else if (nitems != n_items_obs)
    warning("Provided nitems (", nitems, ") differs from ncol(data). Using ncol(data).")
  nitems <- n_items_obs

  # --- Transform to equivalent binomial scale for polytomous ---
  if (score.type == "poly") {
    if (is.null(min.resp) || is.null(max.resp))
      stop("For score.type = 'poly', min.resp and max.resp must be provided.")
    if (min.resp >= max.resp) stop("min.resp must be less than max.resp.")
    # Scale each item to [0,1]
    X_scaled <- (X - min.resp) / (max.resp - min.resp)
    equiv_total <- rowSums(X_scaled)          # ranges 0..nitems
    raw_total   <- rowSums(X)                 # original raw sum
    p_item <- colMeans(X_scaled)
    mean_pq <- mean(p_item * (1 - p_item))
    var_Xi <- stats::var(p_item)
    mu_X <- mean(equiv_total)
    var_Xp <- stats::var(equiv_total)
  } else { # dichotomous
    equiv_total <- rowSums(X)
    raw_total   <- equiv_total
    p_item <- colMeans(X)
    mean_pq <- mean(p_item * (1 - p_item))
    var_Xi <- stats::var(p_item)
    mu_X <- mean(equiv_total)
    var_Xp <- stats::var(equiv_total)
  }

  # --- Correction factor (strong true score) ---
  num_K <- nitems * (nitems - 1) * var_Xi
  den_K <- mu_X * (nitems - mu_X) - var_Xp - nitems * mean_pq
  if (den_K <= 0) {
    warning("Denominator non-positive; CSEM strong may be invalid.")
    correction_factor <- NA_real_
  } else {
    correction_factor <- 1 - (num_K / den_K)
  }

  # --- Score values for which CSEM is computed ---
  if (full.range) {
    if (score.type == "dich") {
      score_vals <- 0:nitems
    } else {
      # For polytomous: take unique observed equiv scores, then add 0 and nitems if not present
      score_vals <- sort(unique(equiv_total))
      if (!any(score_vals == 0)) score_vals <- c(0, score_vals)
      if (!any(score_vals == nitems)) score_vals <- c(score_vals, nitems)
      score_vals <- sort(unique(score_vals))
    }
  } else {
    score_vals <- sort(unique(equiv_total))
  }

  # Base binomial variance and strong CSEM
  base_var <- score_vals * (nitems - score_vals) / (nitems - 1)
  base_var[base_var < 0] <- NA_real_
  var_strong <- base_var * correction_factor
  var_strong[var_strong < 0] <- NA_real_
  csem_strong <- sqrt(var_strong)

  # --- Frequency table (ensure no duplications) ---
  freq_tab <- table(equiv_total)
  freq_vec <- as.numeric(freq_tab[match(score_vals, as.numeric(names(freq_tab)))])
  freq_vec[is.na(freq_vec)] <- 0L

  # --- Map raw scores (for polytomous) ---
  if (score.type == "poly") {
    # For each equiv score value, find the corresponding raw_total(s)
    raw_by_equiv <- tapply(raw_total, equiv_total, unique)
    raw_for_score <- sapply(score_vals, function(v) {
      r <- raw_by_equiv[[as.character(v)]]
      if (is.null(r)) NA else r[1]   # take first if multiple
    })
  } else {
    raw_for_score <- score_vals
  }

  # --- Build score data frame (no duplicates) ---
  score_df <- data.frame(
    raw.score = raw_for_score,
    n = freq_vec,
    csem.strong = round(csem_strong, digits.csem)
  )
  if (score.type == "poly") {
    score_df <- cbind(score_df[, 1, drop = FALSE],
                      equiv.score = round(score_vals, digits.csem),
                      score_df[, -1, drop = FALSE])
  }
  # Remove rows with raw.score NA (if any) – they correspond to unobserved extremes added
  score_df <- score_df[!is.na(score_df$raw.score), ]

  # --- Confidence intervals ---
  if (ci) {
    if (is.null(conf.level)) conf.level <- 0.95
    conf.level <- sort(unique(conf.level))
    if (any(conf.level <= 0 | conf.level >= 1))
      stop("conf.level must be between 0 and 1.")
    n <- nitems
    z_vals <- stats::qnorm((1 + conf.level) / 2)
    ci_mat <- matrix(NA, nrow = nrow(score_df), ncol = 2 * length(conf.level))
    colnames_ci <- character(2 * length(conf.level))
    for (i in seq_along(conf.level)) {
      x <- score_vals  # use the equivalent scores (same length as score_df)
      if (ci.method == "csem") {
        low <- x - z_vals[i] * csem_strong
        high <- x + z_vals[i] * csem_strong
        low[low < 0] <- 0
        high[high > n] <- n
      } else { # wilson
        center <- (x + z_vals[i]^2 / 2) / (n + z_vals[i]^2)
        half <- z_vals[i] * sqrt((x * (n - x) / n) + z_vals[i]^2 / 4) / (n + z_vals[i]^2)
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
    score_df <- cbind(score_df, ci_df)
  }

  # --- Person-level data (optional) ---
  person_df <- NULL
  if (return.person) {
    idx <- match(equiv_total, score_vals)
    person_df <- data.frame(
      id = seq_len(n_persons),
      raw.score = raw_total,
      csem.strong = round(csem_strong[idx], digits.csem)
    )
    if (score.type == "poly") {
      person_df <- cbind(person_df[, 1:2, drop = FALSE],
                         equiv.score = round(equiv_total, digits.csem),
                         person_df[, -c(1:2), drop = FALSE])
    }
    if (ci) {
      for (nm in names(ci_df)) {
        person_df[[nm]] <- round(ci_df[[nm]][idx], digits.csem)
      }
    }
  }

  # --- Summary parameters ---
  sum_df <- NULL
  if (summary) {
    sum_df <- data.frame(
      parameter = c("n_items", "mu_X", "var_Xp", "var_Xi", "mean_pq", "correction_factor"),
      value = round(c(nitems, mu_X, var_Xp, var_Xi, mean_pq, correction_factor), 5)
    )
  }

  # --- Output (no $call, no class attribute for printing call) ---
  out <- list(score = score_df, person = person_df, summary = sum_df)
  # Keep class for potential S3 methods but avoid printing call
  class(out) <- "csemStrong"
  return(out)
}
