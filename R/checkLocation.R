#' Robust Location Comparison for Two Test Halves (Yuen's trimmed t-test)
#'
#' @description
#' Compares the central tendency of two test halves using Yuen's robust paired t-test
#' (trimmed means) and computes the robust effect size (AKP) with confidence interval.
#'
#' @param half1 A numeric matrix or data frame for the first test half.
#' @param half2 A numeric matrix or data frame for the second test half.
#' @param trim Proportion of observations to trim from each tail (default = 0.2).
#' @param conf Numeric. Confidence level (default = 0.95).
#' @param na.rm Logical. If TRUE, rows with missing values are removed.
#' @param install Logical. If TRUE and required package not installed, attempt to install it (default = FALSE).
#'
#' @details
#' The test statistic follows a t distribution with trimmed degrees of freedom.
#' The effect size (AKP) is defined as the trimmed mean difference divided by the
#' winsorized standard deviation of the differences (Algina, Keselman, & Penfield, 2005).
#' The confidence interval for the effect size is obtained using the noncentral t distribution
#' via the \code{MBESS} package. If \code{MBESS} is not installed and \code{install = TRUE},
#' the function attempts to install it. Otherwise, the CI is omitted.
#'
#' @return A data frame with columns:
#'   \code{Test}, \code{Statistic}, \code{p.value},
#'   \code{lwr.ci}, \code{upr.ci} (for trimmed mean difference),
#'   \code{ES}, \code{ES.lwr.ci}, \code{ES.upr.ci}.
#'
#' @export
checkLocation <- function(half1, half2, trim = 0.2, conf = 0.95, na.rm = TRUE, install = FALSE) {

  # --- Ensure MBESS is available ---
  if (!requireNamespace("MBESS", quietly = TRUE)) {
    if (install) {
      message("Package 'MBESS' not found. Attempting to install...")
      install.packages("MBESS")
      if (!requireNamespace("MBESS", quietly = TRUE)) {
        stop("Installation of 'MBESS' failed. Please install it manually.")
      }
    } else {
      warning("Package 'MBESS' not installed. Confidence interval for effect size will be NA. Install with install.packages('MBESS').")
    }
  }

  # --- Resto del código igual (sin cambios) ---
  total1 <- rowSums(half1, na.rm = na.rm)
  total2 <- rowSums(half2, na.rm = na.rm)
  if (na.rm) {
    ok <- stats::complete.cases(total1, total2)
    total1 <- total1[ok]
    total2 <- total2[ok]
  }
  n <- length(total1)
  if (n < 10) warning("Sample size small; Yuen's test may be unstable.")

  diff <- total1 - total2
  g <- floor(trim * n)
  diff_sorted <- sort(diff)
  trim_mean <- mean(diff_sorted[(g+1):(n-g)])
  winsorized <- diff_sorted
  winsorized[1:g] <- diff_sorted[g+1]
  winsorized[(n-g+1):n] <- diff_sorted[n-g]
  winsor_var <- var(winsorized)
  se <- sqrt(winsor_var / (n - 2*g))
  df <- n - 2*g - 1
  t_crit <- qt(1 - (1 - conf)/2, df = df)
  ci_low <- trim_mean - t_crit * se
  ci_up <- trim_mean + t_crit * se
  d_robust <- trim_mean / sqrt(winsor_var)

  # --- CI for effect size (using MBESS if available) ---
  es_ci_low <- es_ci_up <- NA
  if (requireNamespace("MBESS", quietly = TRUE)) {
    ncp <- d_robust * sqrt(n - 2*g)
    ci_lam <- tryCatch({
      MBESS::conf.limits.nct(ncp = ncp, df = df, conf.level = conf)
    }, error = function(e) list(Lower.Limit = NA, Upper.Limit = NA))
    es_ci_low <- ci_lam$Lower.Limit / sqrt(n - 2*g)
    es_ci_up <- ci_lam$Upper.Limit / sqrt(n - 2*g)
  }

  t_stat <- trim_mean / se
  p.value <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)

  result <- data.frame(
    Test = "Yuen's robust t-test",
    Statistic = round(t_stat, 3),
    p.value = round(p.value, 5),
    lwr.ci = round(ci_low, 3),
    upr.ci = round(ci_up, 3),
    ES = round(d_robust, 3),
    ES.lwr.ci = round(es_ci_low, 3),
    ES.upr.ci = round(es_ci_up, 3)
  )
  return(result)
}
