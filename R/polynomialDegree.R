#' Selección del grado polinomial para suavizado de CSEM (regresión lineal)
#'
#' @description
#' Ajusta modelos polinomiales de grado 1 a `max_degree` a los datos `(x, y)`
#' y devuelve una tabla con criterios de información (AIC, BIC), R² ajustado,
#' y significación del término de mayor grado. Ayuda a elegir el grado apropiado
#' para usar en `csemFeldt(..., smooth = TRUE, degree = ...)`.
#'
#' @param x Vector numérico de puntuaciones (centros de grupo).
#' @param y Vector numérico de CSEM.raw (estimaciones sin suavizar).
#' @param min_degree Grado mínimo a evaluar (por defecto 1).
#' @param max_degree Grado máximo a evaluar (por defecto 6).
#' @param show.coefs Si TRUE, muestra los coeficientes de cada modelo.
#' @param plot Si TRUE, genera gráfico de AIC y BIC frente al grado.
#'
#' @return Una lista con:
#' \item{summary}{data.frame con grados, AIC, BIC, R² ajustado, p-valor del
#'   término de mayor grado (comparación con modelo anterior), y si el grado
#'   es "mejor" según AIC o BIC.}
#' \item{models}{Lista de modelos `lm` para cada grado.}
#'
#' @examples
#' set.seed(123)
#' x <- 0:20
#' y <- 2 + 0.5*x - 0.02*x^2 + rnorm(21, sd=0.2)
#' res <- polynomialDegree(x, y, max_degree = 4, plot = TRUE)
#' print(res$summary)
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
    # Polinomio sin centrar (raw = TRUE) para facilitar interpretación
    form <- as.formula(paste("y ~ poly(x,", d, ", raw = TRUE)"))
    mod <- lm(form, data = data.frame(x = x, y = y))
    models_list[[as.character(d)]] <- mod
    s <- summary(mod)
    adjR2 <- s$adj.r.squared
    aic <- AIC(mod)
    bic <- BIC(mod)

    # Test F de mejora respecto al modelo anterior (si existe)
    p_improve <- NA
    if (d > min_degree && !is.null(prev_model)) {
      ftest <- anova(prev_model, mod)
      p_improve <- ftest$`Pr(>F)`[2]
    }

    results[results$degree == d, c("AIC", "BIC", "adjR2", "p_improve")] <-
      c(aic, bic, adjR2, p_improve)

    prev_model <- mod
  }

  # Identificar mejores grados según AIC y BIC (el menor valor)
  results$best_AIC <- results$AIC == min(results$AIC, na.rm = TRUE)
  results$best_BIC <- results$BIC == min(results$BIC, na.rm = TRUE)

  if (plot) {
    graphics::par(mar = c(5, 4, 4, 4))
    plot(results$degree, results$AIC, type = "b", col = "red", pch = 19,
         xlab = "Grado del polinomio", ylab = "AIC / BIC",
         main = "Criterios de información para selección de grado")
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

  # Mensaje resumen
  best_aic_deg <- results$degree[results$best_AIC][1]
  best_bic_deg <- results$degree[results$best_BIC][1]
  cat("\n=== Resumen de modelos polinomiales ===\n")
  print(results[, c("degree", "AIC", "BIC", "adjR2", "p_improve")])
  cat("\nMejor grado según AIC:", best_aic_deg, "(AIC =", round(min(results$AIC), 2), ")")
  cat("\nMejor grado según BIC:", best_bic_deg, "(BIC =", round(min(results$BIC), 2), ")")
  cat("\n")

  invisible(list(summary = results, models = models_list))
}
