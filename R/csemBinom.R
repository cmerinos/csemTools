#' Compute Conditional Standard Error of Measurement (CSEM) Using the Binomial Model
#'
#' @description
#' This function estimates the Conditional Standard Error of Measurement (CSEM)
#' based on the binomial error model. It transforms observed scores into a binomial
#' equivalent, making it applicable for both dichotomous and polytomous items.
#'
#' @param score Numeric vector. Observed total test scores.
#' @param score.type Character string. Indicates whether the score is based on polytomous ("poly")
#' or dichotomous ("dich") items.
#' @param nitems Integer. The total number of items in the test.
#' @param min.resp Integer. The minimum possible response value per item.
#' @param max.resp Integer. The maximum possible response value per item.
#'
#' @details
#' This function transforms observed scores into a proportion and rescales them
#' into a binomial-equivalent score. The standard error of measurement (CSEM) is
#' then computed using the binomial model:
#'
#' \deqn{CSEM = \sqrt{\frac{X (k - X)}{k - 1}}}
#'
#' where:
#' - \eqn{X} is the binomial-equivalent transformed score.
#' - \eqn{k} is the total number of items in the test.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{raw.score}: The observed raw score.
#'   \item \code{prop.score}: The proportion of the total possible score (scaled between 0 and 1).
#'   \item \code{equiv.score}: The binomial-equivalent transformed score.
#'   \item \code{binom.CSEM}: The computed Conditional Standard Error of Measurement (CSEM).
#' }
#'
#' @examples
#' # Example usage with a polytomous scale (0 to 4 per item)
#' csemBinom(score = c(15, 20, 25), score.type = "poly", nitems = 12, min.resp = 0, max.resp = 4)
#'
#' # Example usage with a dichotomous scale (0/1 per item)
#' csemBinom(score = c(5, 7, 9), score.type = "dich", nitems = 12, min.resp = 0, max.resp = 1)
#'
#' @export
csemBinom <- function(score, score.type = c("poly", "dich"), nitems, min.resp, max.resp) {

  # Validate score type
  score.type <- match.arg(score.type)

  # 1. Calculate minimum and maximum possible scores
  minscore <- min.resp * nitems
  maxscore <- max.resp * nitems

  # 2. Validate parameters
  if (!is.numeric(score) || !is.numeric(nitems) || !is.numeric(min.resp) || !is.numeric(max.resp)) {
    stop("⚠️ All arguments must be numeric.")
  }

  if (minscore >= maxscore) {
    stop("⚠️ Minimum score cannot be greater than or equal to the maximum score.")
  }

  if (nitems <= 0) {
    stop("⚠️ The number of items (nitems) must be greater than 0.")
  }

  # 3. If score.type = "dich", use raw scores directly
  if (score.type == "dich") {
    raw.score <- seq(from = minscore, to = maxscore, by = 1)
    equiv.score <- raw.score  # No transformation needed
  } else {
    message("ℹ️ Converting polytomous scores into a binomial-equivalent score.")
    raw.score <- seq(from = minscore, to = maxscore, by = 1)
    prop.score <- (raw.score - minscore) / (maxscore - minscore)
    equiv.score <- prop.score * nitems
  }

  # 4. Compute binomial CSEM
  binom.CSEM <- sqrt((equiv.score * (nitems - equiv.score)) / (nitems - 1))

  # 5. Create final dataframe
  result <- data.frame(
    raw.score = raw.score,
    equiv.score = round(equiv.score, 3),
    binom.CSEM = round(binom.CSEM, 3)
  )

  # 6. Warning if extreme values
  if (any(binom.CSEM == 0)) {
    message("⚠️ Some CSEM values are zero, likely at the scale boundaries.")
  }

  return(result)
}
