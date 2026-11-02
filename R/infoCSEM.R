#' Approximate Information Function Based on Standardized CSEM
#'
#' @description
#' Computes an approximate conditional information function from standardized
#' conditional SEM values as defined by Raju et al. (2007). In Raju's
#' terminology, the standardized CSEM (std-CSEM) is the conditional reliability
#' coefficient:
#' \deqn{stdCSEM(x) = \rho_{cond}(x) = 1 - \frac{CSEM(x)^2}{\sigma_X^2}}
#' The information function is then obtained as the conditional signal-to-noise
#' ratio:
#' \deqn{I(x) = \frac{\rho_{cond}(x)}{1 - \rho_{cond}(x)} =
#'       \frac{stdCSEM(x)}{1 - stdCSEM(x)}}
#'
#' This function provides a CTT-based analogue to the IRT information function.
#' While the IRT information function is defined as \eqn{I(\theta) = 1/SEM(\theta)^2},
#' the present formulation uses the standardized conditional reliability as the
#' basis for the information metric, yielding a measure that is:
#' \itemize{
#'   \item \strong{Adimensional}: expressed as a signal-to-noise ratio, making
#'         it comparable across tests with different score scales.
#'   \item \strong{Interpretable}: values greater than 1 indicate more signal
#'         (true-score variance) than noise (error variance) at that score level;
#'         values less than 1 indicate the opposite.
#'   \item \strong{Consistent with classical theory}: it is the exact conditional
#'         analogue of the global ratio \eqn{\rho/(1-\rho)}.
#' }
#'
#' @details
#' \strong{Theoretical Background}
#'
#' Raju et al. (2007) proposed the standardized conditional SEM as a measure of
#' conditional reliability:
#' \deqn{\rho_{cond}(x) = 1 - \frac{CSEM(x)^2}{\sigma_X^2}}
#' where \eqn{\sigma_X^2} is the total observed-score variance. This formulation
#' ensures that the average of conditional reliabilities equals the group-level
#' reliability coefficient.
#'
#' Nicewander (2018a, 2018b) showed that conditional reliability and the
#' information function are fundamentally related. For a number-correct score
#' \eqn{X} at a fixed value of the latent variable \eqn{\theta}:
#' \deqn{\rho(X, X'|\theta) = \frac{I(X, \theta)}{I(X, \theta) + 1}}
#' which implies:
#' \deqn{I(X, \theta) = \frac{\rho(X, X'|\theta)}{1 - \rho(X, X'|\theta)}}
#'
#' The present function applies this transformation directly to the
#' standardized conditional SEM of Raju et al. (2007), producing a CTT-based
#' approximate information function. The resulting quantity is interpretable
#' as the conditional signal-to-noise ratio at each score level.
#'
#' \strong{Rationale and Utility}
#'
#' Traditional conditional standard errors of measurement (CSEMs) can be
#' misleading when used as the sole indicator of conditional precision.
#' Nicewander (2018a) showed that CSEMs tend to be smaller at the extremes of
#' the score distribution, suggesting greater precision, while the corresponding
#' conditional reliabilities and information functions indicate greater precision
#' in the middle of the distribution. This apparent paradox arises because the
#' CSEM reflects only the absolute error variance, whereas the signal-to-noise
#' ratio (information) reflects the error variance relative to the true-score
#' variance at each score level.
#'
#' The information function computed by \code{infoCSEM} addresses this limitation
#' by providing a measure that:
#' \itemize{
#'   \item Reflects the relative precision of measurement across the score scale.
#'   \item Is on a common metric that facilitates comparisons across different
#'         tests or testing conditions.
#'   \item Complements the std-CSEM by offering a different perspective on
#'         conditional precision.
#' }
#'
#' When used alongside \code{\link{stdCSEM}}, this function provides a more
#' complete picture of measurement precision across the score scale, consistent
#' with the recommendation of Nicewander (2018a) that multiple indicators of
#' conditional precision should be considered jointly.
#'
#' @param stdCSEM Numeric vector. Standardized conditional SEM values as returned
#'        by \code{\link{stdCSEM}} (i.e., conditional reliability estimates
#'        following Raju et al., 2007). Values must be in the \code{[0, 1]} interval.
#' @param score Numeric vector of the same length as \code{stdCSEM}. The score
#'        values corresponding to each std-CSEM estimate.
#' @param bin.score Integer or \code{NULL}. Number of quantile groups for
#'        aggregating the information estimates. If \code{NULL} (default),
#'        information is reported for each individual score level. If an integer
#'        is provided (e.g., 5 for quintiles), scores are grouped into quantiles
#'        and the mean information within each group is reported.
#' @param na.rm Logical. If \code{TRUE}, removes missing values from both
#'        \code{stdCSEM} and \code{score} before computation. Default = \code{FALSE}.
#' @param digits Integer. Number of decimal places for rounding the output.
#'        Default = 3.
#' @param max.info Numeric. Upper bound for information values to avoid infinite
#'        values when \code{stdCSEM} is exactly 1. Default = 999.
#'
#' @return A list with the following components:
#' \item{info}{A \code{data.frame} with columns:
#'   \itemize{
#'     \item \code{score}: the score values.
#'     \item \code{n}: number of observations at each score level (set to
#'           \code{NA} by default; users can merge with their own data).
#'     \item \code{stdCSEM}: the input standardized conditional SEM values
#'           (i.e., conditional reliabilities).
#'     \item \code{info}: the information values
#'           (\eqn{stdCSEM / (1 - stdCSEM)}).
#'   }
#' }
#' \item{binned.info}{If \code{bin.score} is provided, a \code{data.frame}
#'   with columns:
#'   \itemize{
#'     \item \code{group}: quantile group number.
#'     \item \code{range}: score range of the group (min-max).
#'     \item \code{n}: number of observations in the group.
#'     \item \code{mean_score}: mean score within the group.
#'     \item \code{mean_info}: mean information within the group.
#'     \item \code{mean_stdCSEM}: mean std-CSEM within the group.
#'   }
#' }
#'
#' @references
#' Nicewander, W. A. (2018a). Conditional reliability coefficients for test
#' scores. \emph{Psychological Methods}, 23(2), 351-362.
#' \doi{10.1037/met0000132}
#'
#' Nicewander, W. A. (2018b). Conditional precision of measurement for test
#' scores: Are conditional standard errors sufficient? \emph{Educational and
#' Psychological Measurement}, 79(1), 5-18.
#' \doi{10.1177/0013164418758373}
#'
#' Raju, N. S., Price, L. R., Oshima, T. C., & Nering, M. L. (2007).
#' Standardized conditional SEM: A case for conditional reliability.
#' \emph{Applied Psychological Measurement}, 31(3), 169-180.
#' \doi{10.1177/0146621606291569}
#'
#' @seealso \code{\link{stdCSEM}} for computing standardized conditional
#'          reliability from CSEMs.
#'
#' @examples
#' \donttest{
#' # Direct use with known conditional reliabilities
#' rho_cond <- c(0.84, 0.75, 0.64)
#' scores <- 1:3
#'
#' infoCSEM(stdCSEM = rho_cond, score = scores)
#' # info = c(5.25, 3.00, 1.78)
#'
#' # With binning into quantile groups
#' infoCSEM(stdCSEM = rho_cond, score = scores, bin.score = 2)
#'
#' # Integration with stdCSEM
#' csem_vals <- c(2.0, 2.5, 3.0)
#' var_obs <- 25
#' rho_cond <- stdCSEM(csem = csem_vals, var_obs = var_obs)  # returns rho_cond
#' infoCSEM(stdCSEM = rho_cond, score = 1:3)
#' }
#'
#' @export
infoCSEM <- function(stdCSEM,
                     score,
                     bin.score = NULL,
                     na.rm = FALSE,
                     digits = 3,
                     max.info = 999) {

  # --- Input validation ---
  if (!is.numeric(stdCSEM)) {
    stop("'stdCSEM' must be a numeric vector.")
  }
  if (!is.numeric(score)) {
    stop("'score' must be a numeric vector.")
  }
  if (length(stdCSEM) != length(score)) {
    stop("'stdCSEM' and 'score' must have the same length.")
  }

  # --- Handle missing values ---
  if (na.rm) {
    complete_idx <- stats::complete.cases(stdCSEM, score)
    stdCSEM <- stdCSEM[complete_idx]
    score <- score[complete_idx]
  } else if (anyNA(stdCSEM) || anyNA(score)) {
    stop("Missing values detected. Use na.rm = TRUE to remove them.")
  }

  if (length(stdCSEM) == 0) {
    stop("No valid observations after removing missing values.")
  }

  # --- Validate stdCSEM range (must be in [0, 1] as it is a reliability) ---
  if (any(stdCSEM < 0, na.rm = TRUE)) {
    warning("Negative stdCSEM values detected. These will be set to 0.")
    stdCSEM <- pmax(stdCSEM, 0)
  }
  if (any(stdCSEM > 1, na.rm = TRUE)) {
    warning("stdCSEM values > 1 detected. These will be set to 1.")
    stdCSEM <- pmin(stdCSEM, 1)
  }

  # --- Compute information ---
  # stdCSEM here is rho_cond (Raju's standardized CSEM = conditional reliability)
  # info = rho_cond / (1 - rho_cond)
  info <- stdCSEM / (1 - stdCSEM)
  info[is.infinite(info) | info > max.info] <- max.info
  info[is.nan(info)] <- NA_real_

  # --- Prepare output data frame ---
  result_df <- data.frame(
    score = score,
    n = NA_integer_,
    stdCSEM = round(stdCSEM, digits),
    info = round(info, digits),
    stringsAsFactors = FALSE
  )

  # --- Binning by quantiles ---
  binned_df <- NULL
  if (!is.null(bin.score)) {
    if (!is.numeric(bin.score) || length(bin.score) != 1 || bin.score < 2) {
      stop("'bin.score' must be an integer >= 2.")
    }
    bin.score <- as.integer(bin.score)

    q <- stats::quantile(score, probs = seq(0, 1, length.out = bin.score + 1),
                         type = 7, na.rm = TRUE)
    q <- unique(q)
    groups <- cut(score, breaks = q, include.lowest = TRUE, right = TRUE)
    group_levels <- levels(groups)

    bin_list <- list()
    for (i in seq_along(group_levels)) {
      idx_in_group <- which(groups == group_levels[i])
      if (length(idx_in_group) == 0) next

      scores_in_group <- score[idx_in_group]
      info_in_group <- info[idx_in_group]
      std_in_group <- stdCSEM[idx_in_group]

      bin_list[[i]] <- data.frame(
        group = i,
        range = paste0(min(scores_in_group), "-", max(scores_in_group)),
        n = length(idx_in_group),
        mean_score = round(mean(scores_in_group, na.rm = TRUE), digits),
        mean_info = round(mean(info_in_group, na.rm = TRUE), digits),
        mean_stdCSEM = round(mean(std_in_group, na.rm = TRUE), digits),
        stringsAsFactors = FALSE
      )
    }
    binned_df <- do.call(rbind, bin_list)
  }

  # --- Output ---
  out <- list(info = result_df)
  if (!is.null(binned_df)) {
    out$binned.info <- binned_df
  }
  return(out)
}
