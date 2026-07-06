#' Gilmer-Feldt Congeneric Reliability for Each Test Half
#'
#' @description
#' Computes the Gilmer-Feldt congeneric reliability coefficient for each test half
#' separately, with bootstrap confidence intervals. This is a lightweight alternative
#' to McDonald's Omega, suitable for congeneric measures.
#'
#' @param half1 A numeric matrix or data frame with item scores for the first half.
#' @param half2 A numeric matrix or data frame with item scores for the second half.
#' @param B Integer. Number of bootstrap resamples (default = 500).
#' @param conf Numeric. Confidence level (default = 0.95).
#' @param na.rm Logical. If TRUE, rows with missing values are removed.
#'
#' @details
#' The Gilmer-Feldt coefficient (Gilmer & Feldt, 1983) is a classical reliability
#' estimator for congeneric tests, where items measure the same construct but may
#' have different loadings. It is computed from the covariance matrix of the items
#' using a weighting scheme that maximizes internal consistency.
#'
#' Bootstrap confidence intervals are obtained via the percentile method.
#'
#' @return A data frame with columns:
#'   \code{Half} ("Half 1" or "Half 2"),
#'   \code{Coefficient} ("Gilmer-Feldt"),
#'   \code{Estimate} (estimated reliability),
#'   \code{lwr.ci}, \code{upr.ci} (bootstrap confidence interval).
#'
#' @references
#' Gilmer, J. S., & Feldt, L. S. (1983). Reliability estimation for a test
#'   with parts of unknown lengths. \emph{Psychometrika}, 48(1), 99-111.
#'
#' @examples
#' set.seed(123)
#' \donttest{
#' ## Load data
#' library(EFA.dimensions)
#' data("data_RSE")
#'
#' ## Recode negative items
#' data_RSE[c("Q3", "Q5", "Q8", "Q9", "Q10")] <- 5 - data_RSE[c("Q3", "Q5", "Q8", "Q9", "Q10")]
#'
#' ## Check split: difficulty criteria
#' RSE.namesHalf <- checkSplit(data = data_RSE, method = "difficulty")
#'
#'## check Distribution
#' checkCongeneric(half1 = data_RSE[,RSE.namesHalf$half1],
#'                  half2 = data_RSE[,RSE.namesHalf$half2],
#'                  B = 1000, conf = .95)
#' }
#'
#' @importFrom boot boot boot.ci
#' @export
checkCongeneric <- function(half1, half2, B = 500, conf = 0.95, na.rm = TRUE) {

  # --- Helper: Gilmer-Feldt coefficient for a single data set ---
  gf_coef <- function(x) {
    x <- as.matrix(x)
    # Remove any rows with NA (if any left)
    if (anyNA(x)) x <- x[complete.cases(x), , drop = FALSE]
    if (nrow(x) < 3) return(NA_real_)
    m <- stats::cov(x, use = "pairwise.complete.obs")
    k <- ncol(m)
    if (k < 2) return(NA_real_)
    total <- sum(m)
    nondiag <- rowSums(m) - diag(m)
    # Find row with maximum sum of nondiagonal covariances
    max_row <- which.max(nondiag)
    key_row <- m[max_row, ]
    max_nd <- nondiag[max_row]
    D <- ifelse(nondiag == max_nd,
                1,
                (nondiag - key_row) / (max_nd - key_row))
    W <- sum(D^2)
    Q <- sum(D)^2
    (Q / (Q - W)) * (sum(nondiag) / total)
  }

  # --- Bootstrap function for one half ---
  boot_gf <- function(dat) {
    stat_fun <- function(data, idx) gf_coef(data[idx, , drop = FALSE])
    bt <- boot::boot(data = dat, statistic = stat_fun, R = B)
    theta <- bt$t0
    ci <- tryCatch({
      boot::boot.ci(bt, conf = conf, type = "perc")$percent[4:5]
    }, error = function(e) c(NA, NA))
    list(est = theta, lwr = ci[1], upr = ci[2])
  }

  # --- Data preparation and NA handling ---
  if (na.rm) {
    ok1 <- stats::complete.cases(half1)
    ok2 <- stats::complete.cases(half2)
    if (!all(ok1)) half1 <- half1[ok1, , drop = FALSE]
    if (!all(ok2)) half2 <- half2[ok2, , drop = FALSE]
    if (nrow(half1) < 3 || nrow(half2) < 3) stop("Not enough complete cases.")
  } else {
    if (anyNA(half1) || anyNA(half2)) stop("Missing values found. Set na.rm = TRUE.")
  }

  # --- Compute for half1 and half2 ---
  res1 <- boot_gf(half1)
  res2 <- boot_gf(half2)

  # --- Output data frame ---
  result <- data.frame(
    Half = c("Half 1", "Half 2"),
    Coefficient = "Gilmer-Feldt",
    Estimate = round(c(res1$est, res2$est), 4),
    lwr.ci = round(c(res1$lwr, res2$lwr), 4),
    upr.ci = round(c(res1$upr, res2$upr), 4)
  )
  return(result)
}
