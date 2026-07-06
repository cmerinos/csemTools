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
#'
#' @details
#' The test statistic follows a t distribution with trimmed degrees of freedom.
#' The effect size (AKP) is defined as the trimmed mean difference divided by the
#' winsorized standard deviation of the differences (Algina, Keselman, & Penfield, 2005).
#' The confidence interval for the effect size is obtained using the noncentral t distribution
#' via the \code{MBESS} package. If \code{MBESS} is not installed, the CI is omitted
#' and a warning is issued.
#'
#' @references
#' Yuen, K. K. (1974). The two-sample trimmed t for unequal population variances. Biometrika, 61(1), 165-170.
#'
#' Algina, J., Keselman, H. J., & Penfield, R. D. (2005). An alternative to Cohen's standardized mean difference effect size: a robust parameter and confidence interval in the two independent groups case.
#' Psychological methods, 10(3), 317–328. \doi{10.1037/1082-989X.10.3.317}
#'
#' Algina, J., Keselman, H. J., & Penfield, R. D. (2005). Effect Sizes and their Intervals: The Two-Level Repeated Measures Case.
#' Educational and Psychological Measurement, 65(2), 241-258. \doi{10.1177/0013164404268675}
#'
#' @return A data frame with columns:
#'   \code{Test}, \code{Statistic}, \code{p.value},
#'   \code{lwr.ci}, \code{upr.ci} (for trimmed mean difference),
#'   \code{ES}, \code{ES.lwr.ci}, \code{ES.upr.ci}.
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
#' ## check Location
#' checkLocation(half1 = data_RSE[,RSE.namesHalf$half1],
#'               half2 = data_RSE[,RSE.namesHalf$half2],
#'               conf = .95)
#' }
#'
#' @export
checkLocation <- function(half1, half2, trim = 0.2, conf = 0.95, na.rm = TRUE) {

  # --- Check if MBESS is available for effect size CI ---
  has_MBESS <- requireNamespace("MBESS", quietly = TRUE)
  if (!has_MBESS) {
    warning("Package 'MBESS' is not installed. Confidence interval for effect size will be NA. ",
            "Install it with install.packages('MBESS').")
  }

  # --- Compute total scores ---
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
  winsor_var <- stats::var(winsorized)
  se <- sqrt(winsor_var / (n - 2*g))
  df <- n - 2*g - 1
  t_crit <- stats::qt(1 - (1 - conf)/2, df = df)
  ci_low <- trim_mean - t_crit * se
  ci_up <- trim_mean + t_crit * se

  # --- Robust effect size (AKP) ---
  d_robust <- trim_mean / sqrt(winsor_var)

  # --- CI for effect size using MBESS (if available) ---
  es_ci_low <- es_ci_up <- NA
  if (has_MBESS) {
    ncp <- d_robust * sqrt(n - 2*g)
    ci_lam <- tryCatch({
      MBESS::conf.limits.nct(ncp = ncp, df = df, conf.level = conf)
    }, error = function(e) list(Lower.Limit = NA, Upper.Limit = NA))
    es_ci_low <- ci_lam$Lower.Limit / sqrt(n - 2*g)
    es_ci_up <- ci_lam$Upper.Limit / sqrt(n - 2*g)
  }

  # --- Test statistic and p-value ---
  t_stat <- trim_mean / se
  p.value <- 2 * stats::pt(abs(t_stat), df = df, lower.tail = FALSE)

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
