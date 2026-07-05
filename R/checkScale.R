#' Compare Scale (Variability) Between Two Test Halves
#'
#' @description
#' Evaluates differences in variability between two test halves using:
#' \itemize{
#'   \item Log-transformed variability ratio (lnVR) with confidence interval.
#'   \item Bonett-Seier test for equality of variances (robust to non-normality).
#' }
#'
#' @param half1 A numeric matrix or data frame with scores from the first half.
#' @param half2 A numeric matrix or data frame with scores from the second half.
#' @param conf Numeric. Confidence level for lnVR interval (default = 0.95).
#' @param na.rm Logical. If TRUE, rows with missing values are removed.
#'
#' @return A data frame with columns:
#'   \code{lnVR} (log ratio of SDs), \code{SE.lnVR}, \code{lwr.ci}, \code{upr.ci},
#'   \code{cor.halves} (Pearson correlation between halves),
#'   \code{bs.stat} (Bonett-Seier chi-square statistic),
#'   \code{bs.p} (p-value for equality of variances),
#'   \code{n} (effective sample size after NA removal).
#'
#' @references
#' Bonett, D. G., & Seier, E. (2002). Confidence intervals for variance and standard deviation ratios.
#'   Computational Statistics & Data Analysis, 40(3), 603-608.
#'
#' @examples
#' # Simulated data
#' \dontest{
#' half1 <- matrix(rnorm(100*5), ncol=5)
#'
#' half2 <- matrix(rnorm(100*5), ncol=5)
#'
#' checkScale(half1, half2)
#' }
#'
#' @export
checkScale <- function(half1, half2, conf = 0.95, na.rm = TRUE) {
  # Compute total scores
  total1 <- rowSums(half1, na.rm = na.rm)
  total2 <- rowSums(half2, na.rm = na.rm)

  # Remove missing if requested
  if (na.rm) {
    ok <- stats::complete.cases(total1, total2)
    total1 <- total1[ok]
    total2 <- total2[ok]
  }
  n <- length(total1)
  if (n < 4) stop("Need at least 4 complete pairs to compute variability.")

  # Basic statistics
  sd1 <- sd(total1)
  sd2 <- sd(total2)
  r <- cor(total1, total2)

  # --- lnVR and CI (Hedges' approximation) ---
  lnVR <- log(sd1 / sd2)
  SE.lnVR <- sqrt((1 - r^2) / (2 * (n - 3)) + (1/(2*(n-1))) * ((sd1^2/sd2^2) + (sd2^2/sd1^2) - 2))
  z <- qnorm(1 - (1 - conf)/2)
  ci_low <- lnVR - z * SE.lnVR
  ci_up <- lnVR + z * SE.lnVR

  # --- Bonett-Seier test for equality of variances (independent samples version) ---
  # Compute kurtosis (excess) for each half
  kurt1 <- function(x) {
    n <- length(x)
    m4 <- mean((x - mean(x))^4)
    m2 <- mean((x - mean(x))^2)
    (n-1)/((n-2)*(n-3)) * ((n+1)*(m4/m2^2 - 3) + 6)  # unbiased sample excess kurtosis
  }
  g1 <- kurt1(total1)
  g2 <- kurt1(total2)

  # Bonett-Seier statistic (chi-square with 1 df)
  # Equation: x^2 = N * [ln(sd1^2) - ln(sd2^2)]^2 / [4 * (1 - 1/N) * (gamma1 + gamma2 + 2)]
  N <- n  # total sample size (same for both halves because paired)
  var1 <- sd1^2
  var2 <- sd2^2
  numerator <- N * (log(var1) - log(var2))^2
  denominator <- 4 * (1 - 1/N) * (g1 + g2 + 2)
  bs.stat <- numerator / denominator
  bs.p <- 1 - pchisq(bs.stat, df = 1)

  # Result data frame
  result <- data.frame(
    lnVR = round(lnVR, 4),
    SE.lnVR = round(SE.lnVR, 4),
    lwr.ci = round(ci_low, 4),
    upr.ci = round(ci_up, 4),
    cor.halves = round(r, 4),
    bs.stat = round(bs.stat, 4),
    bs.p = round(bs.p, 5),
    n = n
  )
  return(result)
}
