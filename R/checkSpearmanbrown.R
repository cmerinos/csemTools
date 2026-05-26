#' Compute Spearman-Brown Reliability Coefficient
#'
#' @description
#' Calculates the Spearman-Brown reliability coefficient for two test halves
#' and estimates confidence intervals using bootstrap resampling.
#'
#' @param half1 A numeric matrix or data frame containing item scores for the first test half.
#' @param half2 A numeric matrix or data frame containing item scores for the second test half.
#' @param B Integer. Number of bootstrap resamples (default = 500).
#' @param conf Numeric. Confidence level (default = 0.95).
#' @param na.rm Logical. If TRUE (default), rows with any missing values in either half are removed.
#'
#' @details
#' The Spearman-Brown coefficient is computed as:
#' \deqn{r_{SB} = \frac{2 r}{1 + r}}
#' where \eqn{r} is the Pearson correlation between the total scores of the two halves.
#'
#' Confidence intervals are estimated via percentile bootstrap.
#'
#' @return A data frame with columns:
#'   \code{Coefficient} ("Spearman-Brown"),
#'   \code{Estimate} (estimated reliability),
#'   \code{lwr.ci}, \code{upr.ci} (bootstrap confidence interval).
#'
#' @examples
#' set.seed(123)
#' half1 <- matrix(rnorm(100*5), ncol=5)
#' half2 <- matrix(rnorm(100*5), ncol=5)
#' checkSpearmanBrown(half1, half2)
#'
#' @importFrom boot boot boot.ci
#' @export
checkSpearmanBrown <- function(half1, half2, B = 500, conf = 0.95, na.rm = TRUE) {

  # --- Data validation ---
  if (!is.matrix(half1) && !is.data.frame(half1)) stop("'half1' must be a matrix or data frame.")
  if (!is.matrix(half2) && !is.data.frame(half2)) stop("'half2' must be a matrix or data frame.")
  if (nrow(half1) != nrow(half2)) stop("Both halves must have the same number of rows.")

  # --- Handle missing values ---
  if (na.rm) {
    # Remove rows with any NA in either half
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

  # --- Total scores ---
  total1 <- rowSums(half1, na.rm = FALSE)
  total2 <- rowSums(half2, na.rm = FALSE)

  # --- Point estimate of Spearman-Brown ---
  r_halves <- cor(total1, total2, use = "complete.obs")
  sb_est <- (2 * r_halves) / (1 + r_halves)

  # --- Bootstrap for confidence interval ---
  boot_sb <- boot::boot(data = cbind(total1, total2),
                        statistic = function(data, i) {
                          t1 <- data[i, 1]
                          t2 <- data[i, 2]
                          r <- cor(t1, t2, use = "complete.obs")
                          (2 * r) / (1 + r)
                        },
                        R = B)

  ci <- tryCatch({
    boot::boot.ci(boot_sb, conf = conf, type = "perc")$percent[4:5]
  }, error = function(e) c(NA, NA))

  # --- Output ---
  result <- data.frame(
    Coefficient = "Spearman-Brown",
    Estimate = round(sb_est, 4),
    lwr.ci = round(ci[1], 4),
    upr.ci = round(ci[2], 4)
  )
  return(result)
}
