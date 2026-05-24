#' KR-20 Reliability Coefficient for Dichotomous Items
#'
#' @description
#' Computes the Kuder-Richardson formula 20 (KR-20) reliability coefficient
#' for a test with dichotomously scored (0/1) items. This is a measure of
#' internal consistency.
#'
#' @param data A data frame or matrix of dichotomous (0/1) item responses,
#'   with rows as persons and columns as items.
#'
#' @details
#' The KR-20 formula is:
#' \deqn{KR20 = \frac{k}{k-1}\left(1 - \frac{\sum p_j q_j}{\sigma_X^2}\right)}
#' where \eqn{k} is the number of items, \eqn{p_j} is the proportion correct
#' on item \eqn{j}, \eqn{q_j = 1-p_j}, and \eqn{\sigma_X^2} is the variance
#' of the total test scores (using population variance divisor \eqn{N}).
#'
#' Rows with missing values are removed using \code{na.exclude}.
#'
#' @return A single numeric value: the KR-20 reliability coefficient.
#' @export
#'
#' @examples
#' # data.u is a hypothetical dataset of 3000 persons and 40 dichotomous items
#' # kr20(data.u)
kr20 <- function(data) {
  data <- stats::na.exclude(data)
  data <- as.matrix(data)
  if (any(!data %in% c(0, 1))) {
    warning("Non-binary values detected; ensure 0/1 scoring for KR-20.")
  }
  k <- ncol(data)
  N <- nrow(data)
  p <- colMeans(data)
  q <- 1 - p
  sum_pq <- sum(p * q)
  total_scores <- rowSums(data)
  var_total <- var(total_scores) * (N - 1) / N  # population variance
  kr20_val <- (k / (k - 1)) * (1 - sum_pq / var_total)
  return(kr20_val)
}

#' KR-21 Reliability Coefficient for Dichotomous Items
#'
#' @description
#' Computes the Kuder-Richardson formula 21 (KR-21) reliability coefficient,
#' a simplified version of KR-20 that assumes equal item difficulties.
#' It is easier to compute but generally underestimates KR-20.
#'
#' @param data A data frame or matrix of dichotomous (0/1) item responses,
#'   with rows as persons and columns as items.
#'
#' @details
#' The KR-21 formula is:
#' \deqn{KR21 = \frac{k}{k-1}\left(1 - \frac{\bar{X}(k - \bar{X})}{k \sigma_X^2}\right)}
#' where \eqn{k} is the number of items, \eqn{\bar{X}} is the mean total score,
#' and \eqn{\sigma_X^2} is the variance of total scores (population variance).
#'
#' Rows with missing values are removed using \code{na.exclude}.
#'
#' @return A single numeric value: the KR-21 reliability coefficient.
#' @export
#'
#' @examples
#' # data.u is a hypothetical dataset of 3000 persons and 40 dichotomous items
#' # kr21(data.u)
kr21 <- function(data) {
  data <- stats::na.exclude(data)
  data <- as.matrix(data)
  if (any(!data %in% c(0, 1))) {
    warning("Non-binary values detected; ensure 0/1 scoring for KR-21.")
  }
  k <- ncol(data)
  N <- nrow(data)
  total_scores <- rowSums(data)
  mean_total <- mean(total_scores)
  var_total <- var(total_scores) * (N - 1) / N  # population variance
  kr21_val <- (k / (k - 1)) * (1 - (mean_total * (k - mean_total)) / (k * var_total))
  return(kr21_val)
}
