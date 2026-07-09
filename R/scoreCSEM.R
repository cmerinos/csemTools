#' Obtain CSEM for individual scores from a reference table
#'
#' @description
#' Given vectors of reference scores and their CSEM values, this function
#' matches each individual score exactly to a reference score. If an individual
#' score is not found in the reference, a warning is issued and NA is returned.
#' Optionally, it computes confidence intervals assuming normality.
#'
#' @param score.indiv Numeric vector of individual scores.
#' @param score.ref Numeric vector of reference scores (e.g., from a CSEM table).
#' @param csem Numeric vector of CSEM values corresponding to `score.ref`.
#' @param ci Logical. If `TRUE`, compute confidence intervals. Default `FALSE`.
#' @param conf.level Numeric. Confidence level (default 0.95).
#'
#' @return A data frame with the same number of rows as `length(score.indiv)`,
#'   containing columns: `score.indiv`, `CSEM`, and if `ci = TRUE`, `lwr` and `upr`.
#'
#' @examples
#' \donttest{
#' # Use data
#' library(psychTools)
#' data(ability)
#' data.ability <- ability[complete.cases(ability),]
#'
#' # get CSEM deom Binomial model
#' res <- csemStrong(data.ability, score.type = "dich", nitems = 16)
#'
#' # Get CSEM for scores 5, 8, 12 from a strong CSEM table
#' scoreCSEM(score.indiv = c(5, 8, 12),
#'           score.ref = res$CSEM$raw.score,
#'           csem = res$CSEM$csem.strong,
#'           ci = TRUE)
#' }
#'
#' @export
scoreCSEM <- function(score.indiv,
                      score.ref,
                      csem,
                      ci = FALSE,
                      conf.level = 0.95) {

  # --- Validations ---
  if (length(score.ref) != length(csem))
    stop("score.ref and csem must have the same length.")
  if (!is.numeric(score.indiv) || !is.numeric(score.ref) || !is.numeric(csem))
    stop("All inputs must be numeric.")
  if (ci && (conf.level <= 0 || conf.level >= 1))
    stop("conf.level must be between 0 and 1.")

  n <- length(score.indiv)
  csem_out <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    idx <- which(score.ref == score.indiv[i])
    if (length(idx) == 0) {
      warning("Individual score ", score.indiv[i],
              " not found in reference scores. Returning NA.")
      next
    }
    csem_out[i] <- csem[idx[1]]  # take first if multiple (should not happen)
  }

  result <- data.frame(score.indiv = score.indiv, CSEM = csem_out)

  if (ci) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    result$lwr <- score.indiv - z * csem_out
    result$upr <- score.indiv + z * csem_out
  }

  return(result)
}
