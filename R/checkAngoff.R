#' Compute Angoff-Feldt Reliability Coefficient for two test halves
#'
#' @description
#' Calculates the Angoff-Feldt reliability coefficient for two test halves,
#' which is appropriate when the halves may have unequal lengths.
#' Confidence intervals are estimated via bootstrap.
#'
#' @param half1 A numeric matrix or data frame containing item scores for the first test half.
#' @param half2 A numeric matrix or data frame containing item scores for the second test half.
#' @param B Integer. Number of bootstrap resamples (default = 500).
#' @param conf Numeric. Confidence level (default = 0.95).
#' @param na.rm Logical. If TRUE (default), rows with any missing values in either half are removed.
#'
#' @details
#' The Angoff coefficient (Angoff, 1953; Feldt, 2002) is defined as:
#' \deqn{\rho = \frac{4 \text{cov}(X_1, X_2)}{\text{Var}(X_1) + \text{Var}(X_2) + 2 \text{cov}(X_1, X_2)}}
#' where \eqn{X_1} and \eqn{X_2} are the total scores of the two halves.
#'
#' This coefficient is a generalization of the Spearman-Brown formula that does not
#' require equal length halves. For equal lengths, it reduces to the usual
#' Spearman-Brown reliability.
#'
#' @return A data frame with columns:
#'   \code{Coefficient} ("Angoff-Feldt"),
#'   \code{Estimate} (estimated reliability),
#'   \code{lwr.ci}, \code{upr.ci} (bootstrap percentile confidence interval).
#'
#' @references
#' Angoff, W. H. (1953). Test reliability and effective test length.
#'   \emph{Psychometrika}, 18(1), 1-14.
#' Feldt, L. S. (2002). Reliability estimation when a test is split into two parts of unknown effective length.
#'   \emph{Measurement in Education}, 15(3), 295-308.
#'
#' @examples
#' \dontest{
#' set.seed(123)
#'
#' half1 <- matrix(rnorm(100*5), ncol=5)
#'
#' half2 <- matrix(rnorm(100*5), ncol=5)
#'
#' checkAngoff(half1, half2)
#' }
#'
#' @importFrom boot boot boot.ci
#' @export
checkAngoff <- function(half1, half2, B = 500, conf = 0.95, na.rm = TRUE) {

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

  # --- Total scores ---
  total1 <- rowSums(half1, na.rm = FALSE)
  total2 <- rowSums(half2, na.rm = FALSE)

  # --- Function to compute Angoff coefficient from a bootstrap sample ---
  angoff_fun <- function(data, i) {
    t1 <- data[i, 1]
    t2 <- data[i, 2]
    cov_val <- cov(t1, t2, use = "complete.obs")
    var1 <- var(t1, na.rm = TRUE)
    var2 <- var(t2, na.rm = TRUE)
    4 * cov_val / (var1 + var2 + 2 * cov_val)
  }

  # --- Point estimate ---
  est <- angoff_fun(cbind(total1, total2), 1:nrow(half1))

  # --- Bootstrap confidence interval ---
  boot_ang <- boot::boot(data = cbind(total1, total2),
                         statistic = angoff_fun,
                         R = B)

  ci <- tryCatch({
    boot::boot.ci(boot_ang, conf = conf, type = "perc")$percent[4:5]
  }, error = function(e) c(NA, NA))

  # --- Output ---
  result <- data.frame(
    Coefficient = "Angoff-Feldt",
    Estimate = round(est, 4),
    lwr.ci = round(ci[1], 4),
    upr.ci = round(ci[2], 4)
  )
  return(result)
}
