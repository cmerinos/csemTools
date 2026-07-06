#' Compare Distributions of Two Test Halves
#'
#' @description
#' Evaluates distributional similarity between two test halves using:
#' \itemize{
#'   \item Anderson-Darling test (from package \pkg{kSamples}) - optional,
#'   \item Overlapping Index (OVI) based on kernel density estimation,
#'   \item Kendall's W (concordance) derived from Spearman's rho.
#' }
#'
#' @param half1 A numeric matrix or data frame for the first test half.
#' @param half2 A numeric matrix or data frame for the second test half.
#' @param B Integer. Number of bootstrap resamples for confidence intervals of OVI and Kendall's W (default = 500).
#' @param conf Numeric. Confidence level (default = 0.95).
#' @param na.rm Logical. If TRUE, rows with missing values are removed.
#'
#' @details
#' The **Overlapping Index** is computed as the area under the minimum of the two
#' kernel density estimates (using Sheather-Jones bandwidth). The index ranges from 0
#' (no overlap) to 1 (identical distributions).
#'
#' **Kendall's W** for two judges is derived from Spearman's rank correlation:
#' \eqn{W = (1 + \rho_{Spearman}) / 2}. It ranges from 0 (no agreement) to 1 (perfect agreement).
#'
#' The **Anderson-Darling test** (from \pkg{kSamples}) tests the null hypothesis that the two
#' distributions are identical. If the package is not installed, the test is skipped.
#'
#' Effect sizes (ES):
#' \itemize{
#'   \item Anderson-Darling: \eqn{AD / \sqrt{n}} (standardised statistic).
#'   \item OVI: the Overlapping Index itself.
#'   \item Kendall's W: the W coefficient.
#' }
#'
#' @return A data frame with columns:
#'   \code{Test}, \code{Statistic} (test value, NA when not applicable),
#'   \code{p.value} (only for AD), \code{lwr.ci}, \code{upr.ci} (only for OVI and W),
#'   \code{ES} (effect size).
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
#' checkDistribution(half1 = data_RSE[,RSE.namesHalf$half1],
#'                  half2 = data_RSE[,RSE.namesHalf$half2],
#'                  B = 1000, conf = .95)
#' }
#'
#' @importFrom boot boot boot.ci
#'
#' @export
checkDistribution <- function(half1, half2, B = 500, conf = 0.95, na.rm = TRUE) {

  # --- Helper: Overlap area from kernel densities ---
  overlap_area <- function(x, y) {
    # Estimate densities on a common grid
    dens_x <- density(x, na.rm = TRUE)
    dens_y <- density(y, na.rm = TRUE)
    # Common grid covering both ranges
    grid <- seq(min(dens_x$x, dens_y$x), max(dens_x$x, dens_y$x), length.out = 512)
    fx <- approx(dens_x$x, dens_x$y, xout = grid, rule = 2)$y
    fy <- approx(dens_y$x, dens_y$y, xout = grid, rule = 2)$y
    # Area of overlap = integral of min(fx, fy)
    area <- sum(pmin(fx, fy)) * (grid[2] - grid[1])
    return(area)
  }

  # --- Data preparation ---
  total1 <- rowSums(half1, na.rm = na.rm)
  total2 <- rowSums(half2, na.rm = na.rm)
  if (na.rm) {
    ok <- stats::complete.cases(total1, total2)
    total1 <- total1[ok]
    total2 <- total2[ok]
  }
  n <- length(total1)
  if (n < 3) stop("Need at least 3 complete pairs.")

  # --- 1) Anderson-Darling test (kSamples) ---
  ad_stat <- NA; ad_p <- NA; ad_es <- NA
  if (requireNamespace("kSamples", quietly = TRUE)) {
    ad_test <- kSamples::ad.test(total1, total2)
    ad_stat <- round(ad_test$ad[1,1], 3)
    ad_p <- round(ad_test$ad[1,3], 5)
    ad_es <- round(ad_stat / sqrt(n), 3)   # standardised effect size
  } else {
    warning("Package 'kSamples' not installed. Anderson-Darling test omitted. Install it with install.packages('kSamples').")
  }

  # --- 2) Overlapping Index (OVI) with bootstrap CI ---
  ovi_est <- overlap_area(total1, total2)
  # Bootstrap for OVI
  boot_ovi <- boot::boot(data = cbind(total1, total2),
                         statistic = function(d, i) {
                           overlap_area(d[i,1], d[i,2])
                         }, R = B)
  ovi_ci <- tryCatch(boot::boot.ci(boot_ovi, conf = conf, type = "perc")$percent[4:5],
                     error = function(e) c(NA, NA))

  # --- 3) Kendall's W (from Spearman's rho) with bootstrap CI ---
  rho <- cor(total1, total2, method = "spearman", use = "pairwise.complete.obs")
  w_est <- (1 + rho) / 2
  boot_w <- boot::boot(data = cbind(total1, total2),
                       statistic = function(d, i) {
                         rho_i <- cor(d[i,1], d[i,2], method = "spearman", use = "complete.obs")
                         (1 + rho_i) / 2
                       }, R = B)
  w_ci <- tryCatch(boot::boot.ci(boot_w, conf = conf, type = "perc")$percent[4:5],
                   error = function(e) c(NA, NA))

  # --- Build output data frame (with blank instead of NA) ---
  result <- data.frame(
    Test = c("Anderson-Darling", "Overlapping Index (OVI)", "Kendall's W"),
    Statistic = c(ifelse(is.na(ad_stat), "", as.character(ad_stat)), "", ""),
    p.value = c(ifelse(is.na(ad_p), "", as.character(ad_p)), "", ""),
    lwr.ci = c("", ifelse(is.na(ovi_ci[1]), "", as.character(round(ovi_ci[1], 3))),
               ifelse(is.na(w_ci[1]), "", as.character(round(w_ci[1], 3)))),
    upr.ci = c("", ifelse(is.na(ovi_ci[2]), "", as.character(round(ovi_ci[2], 3))),
               ifelse(is.na(w_ci[2]), "", as.character(round(w_ci[2], 3)))),
    ES = c(ifelse(is.na(ad_es), "", as.character(ad_es)),
           ifelse(is.na(ovi_est), "", as.character(round(ovi_est, 3))),
           ifelse(is.na(w_est), "", as.character(round(w_est, 3)))),
    stringsAsFactors = FALSE
  )
  return(result)
}
