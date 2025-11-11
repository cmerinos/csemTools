#' Calculate Conditional Standard Error of Measurement (CSEM) Using Thorndike's Method
#'
#' @description This function estimates the Conditional Standard Error of Measurement (CSEM) 
#' using Thorndike's method, which relies on the standard deviation of half-test differences.
#'
#' @param half1 Data frame with items from the first half of the test.
#' @param half2 Data frame with items from the second half of the test.
#' @param n.groups Integer. Number of groups to divide the total scores into (e.g., 2 for halves, 4 for quartiles, 10 for deciles, etc.). 
#' If `bin.width` is specified, this argument will be ignored.
#' @param bin.width Numeric. If specified, divides the total scores into intervals of the specified width.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{Group}: The group or interval of scores.
#'   \item \code{N}: The number of observations in the group.
#'   \item \code{Mean_Total}: The mean of the total scores in the group.
#'   \item \code{CSEM}: The estimated Conditional Standard Error of Measurement.
#' }
#'
#' @examples
#' # Example usage:
#' Thorndike.CSEM(half1 = dhalf1, half2 = dhalf2, n.groups = 4)
#' Thorndike.CSEM(half1 = dhalf1, half2 = dhalf2, bin.width = 3)
#'
#' @export
Thorndike.CSEM <- function(half1, half2, n.groups = 10, bin.width = NULL) {
  
  # 1️⃣ Calcular diferencias y puntajes totales
  diff_scores <- rowSums(half1) - rowSums(half2)
  total_scores <- rowSums(half1) + rowSums(half2)
  df <- data.frame(total = total_scores, diff = diff_scores)
  
  # 2️⃣ Validar conflictos de argumentos
  if (!is.null(bin.width) && !missing(n.groups)) {
    warning("⚠️ Se ignorará 'n.groups' porque se ha especificado 'bin.width'.")
  }
  
  # 3️⃣ Crear grupos según n.groups o bin.width
  if (!is.null(bin.width)) {
    df$group <- cut(df$total, breaks = seq(min(df$total), max(df$total), by = bin.width), include.lowest = TRUE)
    grouping <- paste("Intervalos de", bin.width, "puntos")
  } else {
    # Dividir en n.groups
    probs <- seq(0, 1, length.out = n.groups + 1)
    df$group <- cut(df$total, breaks = quantile(df$total, probs = probs, na.rm = TRUE), include.lowest = TRUE)
    grouping <- paste("Grupos basados en", n.groups, "partes iguales")
  }
  
  # 4️⃣ Calcular CSEM para cada grupo
  result <- aggregate(diff ~ group, data = df, FUN = sd, na.rm = TRUE)
  counts <- table(df$group)
  means <- aggregate(total ~ group, data = df, FUN = mean, na.rm = TRUE)
  
  # 5️⃣ Combinar resultados
  final_result <- data.frame(
    Group = result$group,
    N = counts,
    Mean_Total = round(means$total, 2),
    CSEM = round(result$diff, 3)
  )
  
  # 6️⃣ Agregar mensaje informativo
  message("📊 Thorndike CSEM calculado usando: ", grouping)
  
  return(final_result)
}
