#' Calculate MF-CSEM (Mollenkopf-Feldt Conditional Standard Error of Measurement)
#'
#' @description
#' This function calculates the Conditional Standard Error of Measurement (CSEM)
#' using the Mollenkopf-Feldt (MF) procedure. The CSEM is computed for each possible score
#' level, based on the selected halves of items.
#'
#' @param half1 A vector specifying the names or indices of the columns corresponding to the first half of the items.
#' @param half2 A vector specifying the names or indices of the columns corresponding to the second half of the items.
#' @param data A data frame containing the item responses for all examinees.
#' @param reliability.coef A numeric value indicating the reliability coefficient to be used in the calculation
#'        (e.g., Cronbach's α, Omega, Angoff, or Spearman-Brown).
#' @param n.items Integer. The total number of items in the test.
#' @param min.score.item Integer. The minimum possible score for a single item.
#' @param max.score.item Integer. The maximum possible score for a single item.
#' @param conf.level Numeric. The confidence level for the confidence intervals (default = 95).
#'
#' @details
#' The MF-CSEM is calculated using a multi-step procedure:
#' 1. Compute squared adjusted differences.
#' 2. Fit a polynomial regression model to predict the squared differences from the total scores.
#' 3. Calculate the MF-CSEM for each possible score.
#' 4. Adjust the CSEM values based on the chosen reliability coefficient.
#' 5. Compute confidence intervals for the CSEM.
#'
#' @return A list containing three data frames:
#' \itemize{
#'   \item \strong{MF.CSEM.parameters}: Summary statistics of the MF-CSEM.
#'   \item \strong{MF.CSEM.score}: MF-CSEM calculated with observed scores.
#'   \item \strong{MF.CSEM.ETS}: MF-CSEM calculated with estimated true scores (ETS).
#' }
#'
#' @examples
#' # Example with a simulated dataset
#' set.seed(123)
#'
#' # Artificial data
#' data <- data.frame(matrix(sample(3:5, 100 * 10, replace = TRUE), ncol = 10))
#'
#' Run
#' MF.CSEM(half1 = 1:5, half2 = 6:10, data = data, reliability.coef = 0.85,
#'         n.items = 10, min.score.item = 1, max.score.item = 5)
#'
#' @export
MF.CSEM <- function(half1, half2, data,
                    reliability.coef,
                    n.items, min.score.item, max.score.item,
                    conf.level = 95) {

  # 1️⃣ Extract the halves
  half1.data <- data[, half1, drop = FALSE]
  half2.data <- data[, half2, drop = FALSE]

  # 2️⃣ Validate dimensions
  if (ncol(half1.data) != ncol(half2.data)) {
    stop("Error: Both halves must have the same number of items.")
  }

  # 3️⃣ Calculate adjusted differences
  message("🔍 Calculating squared adjusted differences...")
  squared_diff <- calc.adjusted.diff(half1.data, half2.data)

  # 4️⃣ Fit polynomial regression model
  message("⚙️ Fitting polynomial regression model...")
  model <- calc.step3.regression(squared_diff, half1.data, half2.data)

  # Validate model
  if (any(is.na(coef(model)))) {
    stop("Error: Model coefficients contain NAs. Check the input data.")
  }

  # 5️⃣ Calculate possible scores
  scores <- seq(from = min.score.item * n.items,
                to = max.score.item * n.items,
                by = 1)

  # 6️⃣ Calculate CSEM for each possible score (passing model explicitly)
  CSEM_values <- sapply(scores, function(x, model) {
    X1 <- x
    X2 <- x^2
    X3 <- x^3
    pred <- predict(model, newdata = data.frame(X1 = X1, X2 = X2, X3 = X3))
    sqrt(max(pred, 0))
  }, model = model)

  # 7️⃣ Calculate Estimated True Scores (ETS)
  total_score <- rowSums(half1.data, na.rm = TRUE) + rowSums(half2.data, na.rm = TRUE)
  mean_total <- mean(total_score)

  ETS_scores <- sapply(scores, function(x) {
    ((x - mean_total) * reliability.coef) + mean_total
  })

  # 8️⃣ Calculate confidence intervals
  CI_band <- CSEM_values * qnorm(1 - (1 - conf.level) / 2)

  # Adjust observed scores CI (lower bound ≥ 0)
  lwr_ci_observed <- pmax(0, scores - CI_band)
  upr_ci_observed <- scores + CI_band

  # Adjust ETS scores CI (lower bound ≥ 0)
  lwr_ci_ETS <- pmax(0, ETS_scores - CI_band)
  upr_ci_ETS <- ETS_scores + CI_band

  # 🔔 Warn if any adjustments occurred
  if (any(lwr_ci_observed < 0) || any(lwr_ci_ETS < 0)) {
    message("⚠️ Some lower confidence intervals were below zero and were adjusted to 0.")
  }

  # 9️⃣ Create dataframes
  MF.CSEM.parameters <- data.frame(
    Score = scores,
    CSEM = round(CSEM_values, 2)
  )

  MF.CSEM.score <- data.frame(
    Score = scores,
    CSEM = round(CSEM_values, 2),
    CI_band = round(CI_band, 2),
    lwr.ci = round(lwr_ci_observed, 3),
    upp.ci = round(upr_ci_observed, 2)
  )

  MF.CSEM.ETS <- data.frame(
    Score = scores,
    ETS = round(ETS_scores, 2),
    CI_band = round(CI_band, 2),
    lwr.ci = round(lwr_ci_ETS, 2),
    upp.ci = round(upr_ci_ETS, 2)
  )

  # 🔟 Return the results
  result <- list(
    MF.CSEM.parameters = MF.CSEM.parameters,
    MF.CSEM.score = MF.CSEM.score,
    MF.CSEM.ETS = MF.CSEM.ETS
  )

  message("✅ Process completed.")
  return(result)
}
