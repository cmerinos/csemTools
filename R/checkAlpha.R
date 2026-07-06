#' Compute Cronbach's Alpha for Two Test Halves
#'
#' @description
#' Calculates Cronbach's alpha reliability coefficient for each test half
#' separately, with bootstrap confidence intervals.
#'
#' @param half1 A numeric matrix or data frame containing item scores for the first test half.
#' @param half2 A numeric matrix or data frame containing item scores for the second test half.
#' @param B Integer. Number of bootstrap resamples (default = 500).
#' @param conf Numeric. Confidence level (default = 0.95).
#' @param na.rm Logical. If TRUE (default), rows with any missing values in either half are removed.
#'
#' @details
#' Cronbach's alpha is computed as:
#' \deqn{\alpha = \frac{k \bar{r}}{1 + (k - 1) \bar{r}}}
#' where \eqn{k} is the number of items and \eqn{\bar{r}} is the average inter-item correlation.
#'
#' Confidence intervals are estimated via percentile bootstrap.
#'
#' @return A data frame with columns:
#'   \code{Half} ("Half 1" or "Half 2"),
#'   \code{Coefficient} ("alpha"),
#'   \code{Estimate} (estimated alpha),
#'   \code{lwr.ci}, \code{upr.ci} (bootstrap confidence interval).
#'
#' @examples
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
#' checkAlpha(half1 = data_RSE[,RSE.namesHalf$half1],
#'                  half2 = data_RSE[,RSE.namesHalf$half2],
#'                  B = 1000, conf = .95)
#' }
#'
#' @importFrom boot boot boot.ci
#'
#' @export
checkAlpha <- function(half1, half2, B = 500, conf = 0.95, na.rm = TRUE) {

  # --- Data validation ---
  if (!is.matrix(half1) && !is.data.frame(half1)) stop("'half1' must be a matrix or data frame.")
  if (!is.matrix(half2) && !is.data.frame(half2)) stop("'half2' must be a matrix or data frame.")
  if (nrow(half1) != nrow(half2)) stop("Both halves must have the same number of rows.")

  # --- Handle missing values ---
  if (na.rm) {
    ok <- stats::complete.cases(half1) & stats::complete.cases(half2)
    if (!all(ok)) {
      half1 <- half1[ok, , drop = FALSE]
      half2 <- half2[ok, , drop = FALSE]
      warning("Rows with missing values were removed (na.rm = TRUE).")
    }
  } else {
    if (anyNA(half1) || anyNA(half2))
      stop("Missing values found. Set na.rm = TRUE to remove incomplete rows.")
  }

  if (nrow(half1) < 3) stop("At least 3 complete rows are required.")

  # --- Function to compute alpha for a given data matrix ---
  compute_alpha <- function(x) {
    k <- ncol(x)
    if (k < 2) return(NA)
    # Use covariance matrix to compute average inter-item correlation
    covmat <- cov(x, use = "pairwise.complete.obs")
    if (anyNA(covmat)) return(NA)
    r_avg <- mean(covmat[lower.tri(covmat)] / sqrt(diag(covmat) %*% t(diag(covmat)))[lower.tri(covmat)])
    # r_avg <- mean(cor(x, use = "pairwise.complete.obs")[lower.tri(...)])  # simpler
    # Actually, use the formula: alpha = (k * r_avg) / (1 + (k-1) * r_avg)
    (k * r_avg) / (1 + (k - 1) * r_avg)
  }

  # --- Alpha for half1 ---
  alpha1 <- compute_alpha(half1)

  # Bootstrap for half1
  boot_alpha1 <- boot::boot(half1, statistic = function(d, i) {
    compute_alpha(d[i, , drop = FALSE])
  }, R = B)

  ci1 <- tryCatch({
    boot::boot.ci(boot_alpha1, conf = conf, type = "perc")$percent[4:5]
  }, error = function(e) c(NA, NA))

  # --- Alpha for half2 ---
  alpha2 <- compute_alpha(half2)

  boot_alpha2 <- boot::boot(half2, statistic = function(d, i) {
    compute_alpha(d[i, , drop = FALSE])
  }, R = B)

  ci2 <- tryCatch({
    boot::boot.ci(boot_alpha2, conf = conf, type = "perc")$percent[4:5]
  }, error = function(e) c(NA, NA))

  # --- Output ---
  result <- data.frame(
    Half = c("Half 1", "Half 2"),
    Coefficient = "alpha",
    Estimate = round(c(alpha1, alpha2), 4),
    lwr.ci = round(c(ci1[1], ci2[1]), 4),
    upr.ci = round(c(ci1[2], ci2[2]), 4)
  )
  return(result)
}
