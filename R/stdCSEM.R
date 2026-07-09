#' Standardized Conditional Reliability (Raju et al., 2007)
#'
#' @description
#' Computes examinee-level (or score-level) reliability from conditional standard
#' errors of measurement (CSEM) and the total observed score variance.
#'
#' @param csem Numeric vector of conditional standard errors (e.g., from `csemMF` or `csemThorndike`).
#' @param var_obs Numeric. Variance of the observed total scores.
#' @param na.rm Logical. If TRUE, missing values in `csem` are removed before computation.
#'
#' @details
#' The standardized conditional reliability for each score level (or examinee) is:
#' \deqn{\rho = 1 - \frac{CSEM^2}{\sigma_X^2}}
#' where \eqn{\sigma_X^2} is the variance of observed scores in the sample.
#' Values are truncated to the 0 to 1 interval.
#'
#' @return A numeric vector of the same length as `csem` (or shorter if `na.rm = TRUE`),
#'   with reliability estimates rounded to 3 decimal places.
#'
#' @references
#' Raju, N. S., Price, L. R., Oshima, T. C., & Nering, M. L. (2007).
#' Standardized conditional SEM: A comparison of methods.
#' \emph{Educational and Psychological Measurement}, 67(6), 903-916.
#'
#' @examples
#' \donttest{
#' # Sample of values
#' csem_vals <- c(2.0, 2.5, 3.0)
#' var_obs <- 25
#' stdCSEM(csem_vals, var_obs)
#'
#' #From Strong true score model
#' #' library(psychTools)
#' data(ability)
#' data.ability <- ability[complete.cases(ability),]
#'
#' strongCSEM.out <- csemStrong(score.type = "dich",
#' data = data.ability,
#' nitems = 16,
#' ci = TRUE,
#' summary = TRUE)
#'
#' # Looking ouput
#' strongCSEM.out
#'
#' # Standardized CSEM
#' stdCSEM(csem = strongCSEM.out$CSEM$csem.strong,
#' var_obs = strongCSEM.out$summary$value[3])
#' }
#'
#' @export
stdCSEM <- function(csem, var_obs, na.rm = FALSE) {

  # Checks
  if (!is.numeric(csem) && !is.numeric(var_obs)) {
    stop("Both 'csem' and 'var_obs' must be numeric.")
  }
  if (length(var_obs) != 1L || var_obs <= 0) {
    stop("'var_obs' must be a single positive number.")
  }

  # NAs
  if (na.rm) {
    csem <- csem[!is.na(csem)]
  } else if (anyNA(csem)) {
    stop("Missing values in 'csem'. Use na.rm = TRUE to remove them.")
  }

  # Calculate
  rel <- 1 - (csem^2) / var_obs
  rel <- pmax(pmin(rel, 1), 0)   # truncar a [0,1]

  # Message if truncation occurred
  if (any(rel == 0 | rel == 1)) {
    message("Note: Some reliability values were truncated to the [0,1] range.")
  }

  # Round to 3 decimal places and return a numeric vector (not a data.frame)
  round(rel, 3)
}
