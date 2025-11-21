#' Binomial Model for Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Computes the Conditional Standard Error of Measurement (CSEM) under
#' Lord's (1955, 1957) binomial error model. This formulation treats the test
#' as a random sample of dichotomously scored items measuring a single ability.
#' The method can be applied to both dichotomous and polytomous items by
#' transforming raw scores into binomial-equivalent scores. Optionally, the
#' function computes Wilson score confidence intervals for the true
#' number-correct score under the binomial model.
#'
#' @param score Numeric vector (retained for compatibility; not used directly
#'   in the current implementation).
#' @param score.type Character string. Type of item scoring:
#'   \code{"dich"} for dichotomous items (default) or \code{"poly"} for
#'   polytomous items. For polytomous items, raw scores are linearly
#'   transformed into binomial-equivalent scores on the 0--\code{nitems} scale.
#' @param nitems Integer. Number of items in the test.
#' @param min.resp Integer. Minimum response value per item.
#' @param max.resp Integer. Maximum response value per item.
#' @param ci Logical. If \code{TRUE}, compute Wilson score confidence intervals
#'   for the true number-correct score. Default is \code{FALSE}.
#' @param conf.level Numeric vector or \code{NULL}. Confidence level(s) for
#'   Wilson intervals (e.g., \code{0.95} or \code{c(0.90, 0.95)}). If
#'   \code{ci = TRUE} and \code{conf.level = NULL}, a default of \code{0.95}
#'   is used. If \code{ci = FALSE}, this argument is ignored.
#' @param digits.csem Integer. Number of decimal places used to round
#'   proportion scores, binomial-equivalent scores, CSEM, and confidence
#'   interval limits. Default is 3.
#'
#' @details
#' Under the binomial error model (Lord, 1955, 1957), each examinee is assumed
#' to possess a true proportion-correct score \eqn{\phi}, and the observed
#' number-correct score \eqn{X} is binomially distributed with parameters
#' \eqn{n} (number of items) and \eqn{\phi}. The conditional error variance
#' for a person with true proportion-correct \eqn{\phi} is:
#' \deqn{
#'   \sigma^2_{E|X} = n \phi (1 - \phi).
#' }
#'
#' Since \eqn{\phi} is unknown, it is replaced by \eqn{\hat\phi = X / n}, yielding
#' Lord's estimator of the conditional error variance:
#' \deqn{
#'   \widehat{\sigma}^2_{E|X} =
#'   \frac{X (n - X)}{n - 1},
#' }
#' and the CSEM is its square root.
#'
#' For polytomous items, raw scores are first rescaled to proportions in
#' \eqn{[0,1]} and subsequently multiplied by \eqn{n} to obtain the
#' binomial-equivalent score \eqn{X}.
#'
#' \strong{Confidence Intervals}
#'
#' When \code{ci = TRUE}, Wilson score confidence intervals (Wilson, 1927)
#' are computed for the true number-correct score:
#' \deqn{
#'   \tau = n \phi.
#' }
#'
#' Wilson's interval is preferred over the Wald interval because it exhibits
#' superior accuracy for binomial proportions, especially for moderate test
#' lengths or extreme scores. For each requested level \eqn{\gamma},
#' the confidence limits are labeled \code{lwr.xx} and \code{upr.xx}, where
#' \code{xx} denotes the level expressed as a percentage (e.g., 95).
#'
#' @return
#' A data frame with:
#' \itemize{
#'   \item \code{raw.score}: Raw score on the original metric.
#'   \item \code{prop.score}: Raw score rescaled to \eqn{[0,1]}.
#'   \item \code{equiv.score}: Binomial-equivalent score on the 0--\code{nitems}
#'     scale.
#'   \item \code{binom.CSEM}: Binomial conditional standard error of measurement.
#'   \item Additional columns for Wilson interval lower and upper bounds,
#'     if \code{ci = TRUE}.
#' }
#'
#' @references
#' Lord, F. M. (1955). Estimating test reliability.
#'   \emph{Educational and Psychological Measurement}, 15, 325–336.
#'
#' Lord, F. M. (1957). Do tests of the same length have the same standard
#'   error of measurement?
#'   \emph{Educational and Psychological Measurement}, 17, 510–521.
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference.
#'   \emph{Journal of the American Statistical Association}, 22, 209–212.
#'
#' Agresti, A. & Coull, B. A. (1998). Approximate is better than “exact”
#'   for interval estimation of binomial proportions.
#'   \emph{The American Statistician}, 52, 119–126.
#'
#' @examples
#' # Dichotomous example with confidence intervals
#' csemBinom(
#'   score      = c(4, 6, 9),
#'   score.type = "dich",
#'   nitems     = 12,
#'   min.resp   = 0,
#'   max.resp   = 1,
#'   ci         = TRUE,
#'   conf.level = 0.95
#' )
#'
#' # Polytomous example (0–4 scale)
#' csemBinom(
#'   score      = c(15, 20, 25),
#'   score.type = "poly",
#'   nitems     = 12,
#'   min.resp   = 0,
#'   max.resp   = 4
#' )
#'
#' @export
csemBinom <- function(
    score,
    score.type = c("poly", "dich"),
    nitems,
    min.resp,
    max.resp,
    ci          = FALSE,
    conf.level  = NULL,
    digits.csem = 3
) {

  # Match score.type
  score.type <- match.arg(score.type)

  # 1. Calculate minimum and maximum possible raw scores
  minscore <- min.resp * nitems
  maxscore <- max.resp * nitems

  # 2. Validate parameters
  if (!is.numeric(nitems) || !is.numeric(min.resp) || !is.numeric(max.resp)) {
    stop("All arguments nitems, min.resp, and max.resp must be numeric.")
  }

  if (minscore >= maxscore) {
    stop("Minimum score cannot be greater than or equal to the maximum score.")
  }

  if (nitems <= 1) {
    stop("The number of items (nitems) must be greater than 1.")
  }

  # 3. Generate all possible raw scores
  raw.score <- seq(from = minscore, to = maxscore, by = 1)

  # 4. Compute proportion and binomial-equivalent scores
  if (score.type == "dich") {
    # For dichotomous items, the raw score is already 0..nitems
    prop.score  <- (raw.score - minscore) / (maxscore - minscore)
    equiv.score <- raw.score
  } else {
    message("ℹ️ Converting polytomous scores into a binomial-equivalent score.")
    prop.score  <- (raw.score - minscore) / (maxscore - minscore)
    equiv.score <- prop.score * nitems
  }

  # 5. Compute binomial CSEM
  binom.var  <- (equiv.score * (nitems - equiv.score)) / (nitems - 1)
  binom.var[binom.var < 0] <- NA_real_
  binom.CSEM <- sqrt(binom.var)

  # 6. Optional confidence intervals (Wilson) for the true number-correct score
  ci_df <- NULL
  if (isTRUE(ci)) {

    # Default conf.level if ci = TRUE and conf.level is NULL
    if (is.null(conf.level)) {
      conf.level <- 0.95
    }

    conf.level <- sort(unique(conf.level))

    if (any(conf.level <= 0 | conf.level >= 1)) {
      stop("`conf.level` must contain values strictly between 0 and 1 (e.g., 0.95).")
    }

    n <- nitems
    ci_mat <- matrix(NA_real_, nrow = length(equiv.score),
                     ncol = 2L * length(conf.level))
    colnames_ci <- character(2L * length(conf.level))

    for (i in seq_along(conf.level)) {
      cl <- conf.level[i]
      z  <- stats::qnorm((1 + cl) / 2)

      # Wilson interval for proportion p = X / n using binomial-equivalent X
      center <- (equiv.score + z^2 / 2) / (n + z^2)
      half_w <- z * sqrt((equiv.score * (n - equiv.score) / n) + z^2 / 4) /
        (n + z^2)

      p_low  <- center - half_w
      p_high <- center + half_w

      # Truncate to [0, 1]
      p_low[p_low < 0]   <- 0
      p_high[p_high > 1] <- 1

      tau_low  <- n * p_low
      tau_high <- n * p_high

      ci_mat[, (2 * i - 1)] <- tau_low
      ci_mat[, (2 * i)]     <- tau_high

      colnames_ci[(2 * i - 1)] <- paste0("lwr.", formatC(cl * 100, format = "f", digits = 0))
      colnames_ci[(2 * i)]     <- paste0("upr.", formatC(cl * 100, format = "f", digits = 0))
    }

    colnames(ci_mat) <- colnames_ci
    ci_df <- as.data.frame(ci_mat)
    ci_df <- as.data.frame(lapply(ci_df, round, digits = digits.csem))
  }

  # 7. Create final data frame
  result <- data.frame(
    raw.score   = raw.score,
    prop.score  = round(prop.score,  digits.csem),
    equiv.score = round(equiv.score, digits.csem),
    binom.CSEM  = round(binom.CSEM, digits.csem)
  )

  if (!is.null(ci_df)) {
    result <- cbind(result, ci_df)
  }

  # 8. Warning if extreme values
  if (any(binom.CSEM == 0, na.rm = TRUE)) {
    message("⚠️ Some CSEM values are zero, likely at the scale boundaries.")
  }

  return(result)
}
