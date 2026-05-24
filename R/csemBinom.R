#' Binomial Model for Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Computes the Conditional Standard Error of Measurement (CSEM) under
#' Lord's (1955, 1957) binomial error model. This formulation treats the test
#' as a random sample of dichotomously scored items measuring a single ability.
#' The method can be applied to both dichotomous and polytomous items by
#' transforming raw scores into binomial-equivalent scores. Optionally, the
#' function computes confidence intervals for the true number-correct score
#' using either the classic CSEM‑based method or the Wilson score method.
#'
#' @param score.type Character string. Type of item scoring:
#'   \code{"dich"} for dichotomous items (default) or \code{"poly"} for
#'   polytomous items.
#' @param nitems Integer. Number of items in the test.
#' @param min.resp Numeric. Minimum response value per item.
#'   **Required only when \code{score.type = "poly"}. Ignored for dichotomous.**
#' @param max.resp Numeric. Maximum response value per item.
#'   **Required only when \code{score.type = "poly"}. Ignored for dichotomous.**
#' @param ci Logical. If \code{TRUE}, compute confidence intervals for the
#'   true number-correct score. Default is \code{FALSE}.
#' @param ci.method Character. Method for confidence intervals:
#'   \code{"csem"} (classic: \eqn{X \pm z \cdot CSEM}) or
#'   \code{"wilson"} (Wilson score interval). Default is \code{"csem"}.
#'   Only used if \code{ci = TRUE}.
#' @param conf.level Numeric vector or \code{NULL}. Confidence level(s)
#'   (e.g., \code{0.95} or \code{c(0.90, 0.95)}). If \code{ci = TRUE} and
#'   \code{conf.level = NULL}, a default of \code{0.95} is used.
#' @param digits.csem Integer. Number of decimal places used to round
#'   proportion scores, binomial‑equivalent scores, CSEM, and confidence
#'   interval limits. Default is 3.
#'
#' @details
#' \strong{Binomial error model (Lord, 1955, 1957)}:
#' Each examinee is assumed to possess a true proportion‑correct score
#' \eqn{\phi}, and the observed number‑correct score \eqn{X} is binomially
#' distributed with parameters \eqn{n} (number of items) and \eqn{\phi}.
#' The conditional error variance for a person with true proportion‑correct
#' \eqn{\phi} is:
#' \deqn{
#'   \sigma^2_{E|X} = n \phi (1 - \phi).
#' }
#'
#' Since \eqn{\phi} is unknown, it is replaced by \eqn{\hat\phi = X / n},
#' yielding Lord's estimator of the conditional error variance:
#' \deqn{
#'   \widehat{\sigma}^2_{E|X} = \frac{X (n - X)}{n - 1},
#' }
#' and the CSEM is its square root.
#'
#' \strong{Extension to polytomous items}:
#' Raw scores are first rescaled to proportions in \eqn{[0,1]} and subsequently
#' multiplied by \eqn{n} to obtain the binomial‑equivalent score \eqn{X}:
#' \deqn{
#'   \text{equiv.score} = \frac{\text{raw.score} - \text{min.resp} \times n}
#'   {\text{max.resp} \times n - \text{min.resp} \times n} \times n.
#' }
#' This transformation assumes equal scoring weights and that the original
#' score scale is approximately linear with respect to the underlying ability.
#'
#' \strong{Confidence intervals}:
#' When \code{ci = TRUE}, two methods are available:
#' \itemize{
#'   \item \code{"csem"}: Classic CTT interval for the true score
#'     \eqn{\tau = n\phi}:
#'     \deqn{\text{equiv.score} \pm z_{\alpha/2} \cdot \text{CSEM},}
#'     where \eqn{z_{\alpha/2}} is the normal quantile for the given confidence
#'     level. Limits are truncated to \eqn{[0, n]}.
#'   \item \code{"wilson"}: Wilson score interval (Wilson, 1927) for the
#'     binomial proportion \eqn{p = \text{equiv.score}/n}:
#'     \deqn{
#'       \frac{p + \frac{z^2}{2n} \pm z \sqrt{\frac{p(1-p)}{n} + \frac{z^2}{4n^2}}}
#'       {1 + \frac{z^2}{n}},
#'     }
#'     then multiplied by \eqn{n} to obtain limits for the true number‑correct
#'     score. Wilson's interval is preferred over the Wald interval because it
#'     exhibits superior accuracy for binomial proportions, especially for
#'     moderate test lengths or extreme scores.
#' }
#' For each requested level \eqn{\gamma}, the confidence limits are labeled
#' \code{lwr.xx} and \code{upr.xx}, where \code{xx} denotes the level expressed
#' as a percentage (e.g., 95).
#'
#' @return
#' A data frame with one row per possible raw score (from
#' \code{min.resp*nitems} to \code{max.resp*nitems}). Columns include:
#' \itemize{
#'   \item \code{raw.score}: Raw score on the original metric.
#'   \item \code{prop.score}: Raw score rescaled to \eqn{[0,1]}.
#'   \item \code{binom.CSEM}: Binomial conditional standard error of measurement.
#'   \item \code{equiv.score}: (Only for \code{score.type = "poly"})
#'     Binomial‑equivalent score on the 0–\code{nitems} scale.
#'   \item Additional columns for confidence limits, if \code{ci = TRUE}.
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
#' # Dichotomous example (0/1 items) – classic CSEM-based 95% CI
#' csemBinom(
#'   score.type = "dich",
#'   nitems     = 12,
#'   ci         = TRUE,
#'   ci.method  = "csem",
#'   conf.level = 0.95
#' )
#'
#' # Dichotomous example – Wilson 90% and 95% CIs
#' csemBinom(
#'   score.type = "dich",
#'   nitems     = 12,
#'   ci         = TRUE,
#'   ci.method  = "wilson",
#'   conf.level = c(0.90, 0.95)
#' )
#'
#' # Polytomous example (0–4 scale)
#' csemBinom(
#'   score.type = "poly",
#'   nitems     = 12,
#'   min.resp   = 0,
#'   max.resp   = 4
#' )
#'
#' @export
csemBinom <- function(
    score.type = c("dich", "poly"),
    nitems,
    min.resp = NULL,
    max.resp = NULL,
    ci = FALSE,
    ci.method = c("csem", "wilson"),
    conf.level = NULL,
    digits.csem = 3
) {
  score.type <- match.arg(score.type)
  ci.method <- match.arg(ci.method)

  if (score.type == "dich") {
    # For dichotomous, ignore min.resp and max.resp
    if (!is.null(min.resp) || !is.null(max.resp)) {
      message("Note: 'min.resp' and 'max.resp' are ignored for score.type = 'dich' (set to 0 and 1).")
    }
    min.resp <- 0
    max.resp <- 1
  } else {
    # Polytomous: require min.resp and max.resp
    if (is.null(min.resp) || is.null(max.resp)) {
      stop("For score.type = 'poly', you must provide min.resp and max.resp.")
    }
    if (min.resp >= max.resp) stop("min.resp must be less than max.resp.")
  }

  # Validate nitems
  if (!is.numeric(nitems) || length(nitems) != 1 || nitems < 2) {
    stop("nitems must be a single integer >= 2.")
  }
  nitems <- as.integer(nitems)

  minscore <- min.resp * nitems
  maxscore <- max.resp * nitems
  raw.score <- seq(minscore, maxscore, by = 1)

  # Proportion and equivalent score
  prop.score <- (raw.score - minscore) / (maxscore - minscore)
  if (score.type == "dich") {
    equiv.score <- raw.score   # but we won't show it
  } else {
    equiv.score <- prop.score * nitems
  }

  # CSEM
  var <- (equiv.score * (nitems - equiv.score)) / (nitems - 1)
  var[var < 0] <- NA_real_
  csem <- sqrt(var)

  # Prepare result data frame
  result <- data.frame(
    raw.score = raw.score,
    prop.score = round(prop.score, digits.csem),
    binom.CSEM = round(csem, digits.csem)
  )
  if (score.type == "poly") {
    result <- cbind(result, equiv.score = round(equiv.score, digits.csem))
  }

  # Confidence intervals
  if (ci) {
    if (is.null(conf.level)) conf.level <- 0.95
    conf.level <- sort(unique(conf.level))
    if (any(conf.level <= 0 | conf.level >= 1))
      stop("conf.level must be between 0 and 1.")

    n <- nitems
    z <- stats::qnorm((1 + conf.level)/2)
    ci_mat <- matrix(NA, nrow = length(raw.score), ncol = 2*length(conf.level))
    colnames_ci <- character(2*length(conf.level))

    for (i in seq_along(conf.level)) {
      if (ci.method == "csem") {
        low <- equiv.score - z[i] * csem
        high <- equiv.score + z[i] * csem
        low[low < 0] <- 0
        high[high > n] <- n
      } else { # wilson
        center <- (equiv.score + z[i]^2/2) / (n + z[i]^2)
        half <- z[i] * sqrt((equiv.score*(n - equiv.score)/n) + z[i]^2/4) / (n + z[i]^2)
        p_low <- center - half
        p_high <- center + half
        p_low[p_low < 0] <- 0
        p_high[p_high > 1] <- 1
        low <- n * p_low
        high <- n * p_high
      }
      ci_mat[, 2*i - 1] <- low
      ci_mat[, 2*i] <- high
      colnames_ci[2*i - 1] <- paste0("lwr.", formatC(conf.level[i]*100, format="f", digits=0))
      colnames_ci[2*i] <- paste0("upr.", formatC(conf.level[i]*100, format="f", digits=0))
    }
    colnames(ci_mat) <- colnames_ci
    result <- cbind(result, round(ci_mat, digits.csem))
  }

  # Warning for zero CSEM
  if (any(csem == 0, na.rm=TRUE)) message("⚠️ Some CSEM values are zero at boundaries.")

  return(result)
}
