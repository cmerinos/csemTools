#' Binomial Model for Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Computes the Conditional Standard Error of Measurement (CSEM) under
#' Lord's (1955, 1957) binomial error model. This formulation treats the test
#' as a random sample of dichotomously scored items measuring a single ability.
#' The method can be applied to both dichotomous and polytomous items by
#' transforming raw scores into binomial‑equivalent scores. Optionally, the
#' function computes confidence intervals for the true number‑correct score
#' using either the classic CSEM‑based method or the Wilson score method.
#'
#' @param score.type Character: \code{"dich"} (default) or \code{"poly"}.
#' @param nitems Integer. Number of items in the test.
#' @param min.resp Numeric. Minimum response value per item.
#'   **Required when \code{score.type = "poly"}**; ignored for dichotomous.
#' @param max.resp Numeric. Maximum response value per item.
#'   **Required when \code{score.type = "poly"}**; ignored for dichotomous.
#' @param csem.method Character: \code{"Lord"} (default), \code{"binom"} (alias, same as Lord),
#'   or \code{"LordKeats"} (adjusts Lord's CSEM using KR‑20/KR‑21 or empirical factor).
#' @param data Optional data frame or matrix of dichotomous (0/1) item responses.
#'   Required for \code{csem.method = "LordKeats"} if \code{rhoxx} is not provided.
#'   Only used when \code{score.type = "dich"}.
#' @param rhoxx Optional numeric value (0 < rhoxx < 1). An estimate of test reliability
#'   (interpreted as KR‑20). Used for \code{csem.method = "LordKeats"} when no \code{data}
#'   is supplied. Applies the empirical adjustment KR‑21 ≈ 0.8 × KR‑20.
#' @param ci Logical. If \code{TRUE}, compute confidence intervals. Default \code{FALSE}.
#' @param ci.method Character: \code{"csem"} (classic) or \code{"wilson"}. Default \code{"csem"}.
#' @param conf.level Numeric vector of confidence levels (e.g., 0.95). Default \code{NULL} → 0.95.
#' @param digits.csem Integer for rounding. Default 3.
#'
#' @details
#' \strong{Binomial error model (Lord, 1955, 1957)}:
#' Each examinee has true proportion‑correct \eqn{\phi}; observed number‑correct \eqn{X}
#' is Binomial(\eqn{n},\eqn{\phi}) with \eqn{n = \code{nitems}}. Lord's estimator of the
#' conditional error variance is
#' \deqn{\widehat{\sigma}^2_{E|X} = \frac{X (n - X)}{n - 1},}
#' and the CSEM is its square root. The options \code{"Lord"} and \code{"binom"} yield
#' numerically identical results.
#'
#' \strong{Polytomous items}:
#' Raw scores are linearly transformed to a proportion and then multiplied by \eqn{n}
#' to obtain a binomial‑equivalent score:
#' \deqn{\text{equiv} = n \times \frac{\text{raw} - n\cdot\text{min.resp}}{n\cdot(\text{max.resp} - \text{min.resp})}.}
#'
#' \strong{Lord‑Keats method} (\code{csem.method = "LordKeats"}):
#' Adjusts Lord's CSEM using the ratio \eqn{\sqrt{(1-\text{KR20})/(1-\text{KR21})}}.
#' If a data matrix is provided, KR‑20 and KR‑21 are computed directly.
#' If only \code{rhoxx} (interpreted as KR‑20) is given, KR‑21 is approximated as
#' \code{0.8 * rhoxx} following Wilson, Downing & Ebel (1979). This approximation works
#' for both dichotomous and polytomous (after transformation) cases.
#'
#' \strong{Confidence intervals}:
#' \itemize{
#'   \item \code{ci.method = "csem"}: \code{equiv.score ± z * CSEM}.
#'   \item \code{ci.method = "wilson"}: Wilson interval for proportion \eqn{p = equiv.score/n}.
#' }
#'
#' @return
#' A data frame with one row per possible raw score. Columns:
#' \code{raw.score}, \code{prop.score}, \code{binom.CSEM}, and for \code{score.type = "poly"}
#' also \code{equiv.score}. If \code{ci = TRUE}, additional columns \code{lwr.xx}, \code{upr.xx}.
#'
#' @references
#' Lord, F. M. (1955). Estimating test reliability. Educational and Psychological Measurement, 15, 325-336.
#'
#' Lord, F.M. (1957). Do Tests of the Same Length Have the Same Standard Errors of Measurement?
#' Educational and Psychological Measurement, 17, 510 - 521.
#'
#' Keats, J. A. (1957). Estimation of error variances of test scores. Psychometrika, 22, 29-41.
#'
#' Wilson, R. A., Downing, S. M., & Ebel, R. L. (1979). *An empirical adjustment of the Kuder‑Richardson 21 reliability coefficient.
#' https://eric.ed.gov/?id=ED173387
#'
#' Frisbie, D. A. (1988). Reliability of scores from teacher‑made tests.
#' *Educational Measurement: Issues and Practice, 7*(1), 25–35.
#' https://doi.org/10.1111/j.1745-3992.1988.tb00422.x
#'
#' @examples
#' # Dichotomous, Lord method
#' csemBinom(score.type = "dich", nitems = 40)
#'
#' # Polytomous (0–4 scale) with Lord‑Keats using rhoxx
#' csemBinom(score.type = "poly", nitems = 10, min.resp = 0, max.resp = 4,
#'           csem.method = "LordKeats", rhoxx = 0.85)
#'
#' @export
csemBinom <- function(score.type = c("dich", "poly"),
                      nitems,
                      min.resp = NULL,
                      max.resp = NULL,
                      csem.method = c("Lord", "binom", "LordKeats"),
                      data = NULL,
                      rhoxx = NULL,
                      ci = FALSE,
                      ci.method = c("csem", "wilson"),
                      conf.level = NULL,
                      digits.csem = 3) {

  # --- Argument matching and validation ---
  score.type <- match.arg(score.type)
  csem.method <- match.arg(csem.method)
  ci.method <- match.arg(ci.method)

  if (!is.numeric(nitems) || length(nitems) != 1 || nitems < 2)
    stop("nitems must be a single integer >= 2.")
  nitems <- as.integer(nitems)

  # Handle polytomous vs dichotomous
  if (score.type == "poly") {
    if (is.null(min.resp) || is.null(max.resp))
      stop("For score.type = 'poly', min.resp and max.resp must be provided.")
    if (min.resp >= max.resp)
      stop("min.resp must be less than max.resp.")
    minscore <- min.resp * nitems
    maxscore <- max.resp * nitems
  } else { # dich
    # For dichotomous, ignore user-supplied min.resp/max.resp (set internally)
    if (!is.null(min.resp) || !is.null(max.resp))
      message("Note: min.resp and max.resp are ignored when score.type = 'dich' (set to 0 and 1).")
    min.resp <- 0
    max.resp <- 1
    minscore <- 0
    maxscore <- nitems
  }

  # Possible raw scores (full range)
  raw.score <- seq(minscore, maxscore, by = 1)

  # Proportion and equivalent score
  prop.score <- (raw.score - minscore) / (maxscore - minscore)
  if (score.type == "dich") {
    equiv.score <- raw.score   # same as raw
  } else {
    equiv.score <- prop.score * nitems
  }

  # --- Base CSEM (Lord estimator) ---
  var_lord <- (equiv.score * (nitems - equiv.score)) / (nitems - 1)
  var_lord[var_lord < 0] <- NA_real_
  csem_lord <- sqrt(var_lord)

  # --- Adjust if LordKeats ---
  csem_final <- csem_lord
  if (csem.method == "LordKeats") {
    if (score.type != "dich" && is.null(rhoxx) && is.null(data)) {
      # For polytomous without data, rhoxx must be provided
      stop("For LordKeats with polytomous items, you must supply rhoxx (reliability estimate).")
    }
    if (!is.null(data)) {
      # Use empirical KR-20 and KR-21 from data (dichotomous only)
      if (score.type != "dich")
        stop("Data-based LordKeats only supported for dichotomous items (score.type = 'dich').")
      # Compute KR-20 and KR-21 (assuming helper functions exist in package)
      if (!exists("kr20", mode = "function") || !exists("kr21", mode = "function"))
        stop("Functions kr20() and kr21() are required for data-based LordKeats.")
      kr20_val <- kr20(data)
      kr21_val <- kr21(data)
      if (kr20_val <= 0 || kr20_val >= 1 || kr21_val <= 0 || kr21_val >= 1)
        stop("KR-20 and KR-21 must be strictly between 0 and 1.")
      scale_factor <- sqrt((1 - kr20_val) / (1 - kr21_val))
    } else if (!is.null(rhoxx)) {
      # Use manual reliability estimate (interpreted as KR-20) and apply empirical correction
      if (rhoxx <= 0 || rhoxx >= 1)
        stop("rhoxx must be between 0 and 1 (exclusive).")
      kr20_est <- rhoxx
      kr21_est <- rhoxx * 0.8   # empirical adjustment from Wilson et al. (1979)
      scale_factor <- sqrt((1 - kr20_est) / (1 - kr21_est))
      message("LordKeats: Using empirical adjustment KR21 ≈ 0.8 * KR20 (Wilson, Downing & Ebel, 1979).")
    } else {
      stop("For csem.method = 'LordKeats', you must provide either 'data' or 'rhoxx'.")
    }
    csem_final <- csem_lord * scale_factor
  }

  # --- Confidence intervals ---
  ci_df <- NULL
  if (ci) {
    if (is.null(conf.level)) conf.level <- 0.95
    conf.level <- sort(unique(conf.level))
    if (any(conf.level <= 0 | conf.level >= 1))
      stop("conf.level must contain values between 0 and 1.")

    n <- nitems
    z_vals <- stats::qnorm((1 + conf.level) / 2)
    ci_mat <- matrix(NA, nrow = length(equiv.score), ncol = 2 * length(conf.level))
    colnames_ci <- character(2 * length(conf.level))

    for (i in seq_along(conf.level)) {
      if (ci.method == "csem") {
        low <- equiv.score - z_vals[i] * csem_final
        high <- equiv.score + z_vals[i] * csem_final
        low[low < 0] <- 0
        high[high > n] <- n
      } else { # wilson
        center <- (equiv.score + z_vals[i]^2 / 2) / (n + z_vals[i]^2)
        half <- z_vals[i] * sqrt((equiv.score * (n - equiv.score) / n) + z_vals[i]^2 / 4) / (n + z_vals[i]^2)
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

  # --- Construct result data frame ---
  result <- data.frame(
    raw.score = raw.score,
    prop.score = round(prop.score, digits.csem),
    binom.CSEM = round(csem_final, digits.csem)
  )
  if (score.type == "poly") {
    result <- cbind(result, equiv.score = round(equiv.score, digits.csem))
  }
  if (!is.null(ci_df)) {
    result <- cbind(result, ci_df)
  }

  # Warning for zero CSEM at extremes
  if (any(csem_final == 0, na.rm = TRUE))
    message("⚠️ Some CSEM values are zero, likely at the scale boundaries.")

  return(result)
}
