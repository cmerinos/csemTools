#' Selecting the Polynomial Degree for CSEM Smoothing (Linear Regression)
#'
#' @description
#' Fits polynomial models of degree 1 to `max_degree` to the data `(x, y)`
#' and returns a table containing information criteria (AIC, BIC), adjusted R^2,
#' and the significance of the term of the highest degree. This helps in choosing the appropriate degree
#' to use in `csemFeldt(..., smooth = TRUE, degree = ...)`.
#'
#' @param x Numeric vector of scores (cluster centers).
#' @param y Numeric vector of CSEM.raw (unsmoothed estimates).
#' @param min_degree Minimum degree to evaluate (default 1).
#' @param max_degree Maximum degree to evaluate (default 6).
#' @param show.coefs If TRUE, displays the coefficients for each model.
#' @param plot If TRUE, generates a plot of AIC and BIC versus degree.
#'
#' @return A list containing:
#' \item{summary}{A data.frame with degrees of freedom, AIC, BIC, adjusted R^2, p-value of the
#'   term with the highest degrees of freedom (compared to the previous model), and whether the degrees of freedom
#'   are "better" according to AIC or BIC.}
#' \item{models}{A list of `lm` models for each degrees of freedom.}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- 0:20
#' y <- 2 + 0.5*x - 0.02*x^2 + rnorm(21, sd=0.2)
#'
#' res <- polynomialDegree(x, y, max_degree = 4, plot = TRUE)
#'
#' print(res$summary)
#' }
#'
#' @importFrom stats lm cor sd predict AIC BIC aggregate anova approx density pchisq complete.cases weighted.mean setNames
#' @importFrom graphics lines legend
#' @importFrom utils head
#'
#' @export
polynomialDegree <- function(x, y,
                             min_degree = 1,
                             max_degree = 6,
                             show.coefs = FALSE,
                             plot = TRUE) {
  if (length(x) != length(y)) stop("x and y must have the same length.")
  if (min_degree < 1) stop("min_degree must be at least 1.")
  if (max_degree > length(unique(x)) - 1) {
    warning("max_degree reduced to ", length(unique(x))-1, " due to insufficient unique x values.")
    max_degree <- length(unique(x)) - 1
  }
  if (min_degree > max_degree) stop("min_degree > max_degree.")

  degrees <- seq(min_degree, max_degree, by = 1)
  results <- data.frame(
    degree = degrees,
    AIC = NA_real_,
    BIC = NA_real_,
    adjR2 = NA_real_,
    p_improve = NA_real_,  # p-value of F-test comparing to previous degree (NA for min_degree)
    best_AIC = FALSE,
    best_BIC = FALSE,
    stringsAsFactors = FALSE
  )

  models_list <- list()
  prev_model <- NULL

  for (d in degrees) {
    # Uncentered polynomial (raw = TRUE) to facilitate interpretation

    form <- as.formula(paste("y ~ poly(x,", d, ", raw = TRUE)"))
    mod <- lm(form, data = data.frame(x = x, y = y))
    models_list[[as.character(d)]] <- mod
    s <- summary(mod)
    adjR2 <- s$adj.r.squared
    aic <- AIC(mod)
    bic <- BIC(mod)

    # F-test for improvement over the previous model (if applicable)
    p_improve <- NA
    if (d > min_degree && !is.null(prev_model)) {
      ftest <- anova(prev_model, mod)
      p_improve <- ftest$`Pr(>F)`[2]
    }

    results[results$degree == d, c("AIC", "BIC", "adjR2", "p_improve")] <-
      c(aic, bic, adjR2, p_improve)

    prev_model <- mod
  }

  # Identify the best models based on AIC and BIC (the lowest value)
  results$best_AIC <- results$AIC == min(results$AIC, na.rm = TRUE)
  results$best_BIC <- results$BIC == min(results$BIC, na.rm = TRUE)

  if (plot) {
    graphics::par(mar = c(5, 4, 4, 4))
    plot(results$degree, results$AIC, type = "b", col = "red", pch = 19,
         xlab = "Degree of the polynomial", ylab = "AIC / BIC",
         main = "Information Criteria for Grade Selection")
    lines(results$degree, results$BIC, type = "b", col = "blue", pch = 17)
    legend("topright", legend = c("AIC", "BIC"), col = c("red", "blue"),
           lty = 1, pch = c(19, 17), bty = "n")
    graphics::grid()
  }

  if (show.coefs) {
    for (d in degrees) {
      cat("\n=== Grado", d, "===\n")
      print(summary(models_list[[as.character(d)]])$coefficients)
    }
  }

  # Summary Message
  best_aic_deg <- results$degree[results$best_AIC][1]
  best_bic_deg <- results$degree[results$best_BIC][1]
  cat("\n=== Summary of Polynomial Models ===\n")
  print(results[, c("degree", "AIC", "BIC", "adjR2", "p_improve")])
  cat("\nBest degree according to AIC:", best_aic_deg, "(AIC =", round(min(results$AIC), 2), ")")
  cat("\nBest degree according to BIC:", best_bic_deg, "(BIC =", round(min(results$BIC), 2), ")")
  cat("\n")

  invisible(list(summary = results, models = models_list))
}
