#' Strong True Score Model for Conditional Standard Errors of Measurement
#'
#' @description
#' Computes conditional standard errors of measurement (CSEM) under Lord's
#' strong true score theory, also known as the compound binomial model,
#' for dichotomously scored items. The function implements the estimator
#' described by Lord (1965) and summarized by Tong and Kolen (2018), and
#' can optionally provide Wilson score confidence intervals for the true
#' number-correct score.
#'
#' @param data A numeric matrix or data frame with examinees in rows and
#'   dichotomous items (0/1) in columns.
#' @param na.rm Logical. If \code{TRUE} (default), rows with any missing
#'   values are removed before computing the CSEM. If \code{FALSE}, the
#'   function stops with an error when missing values are present.
#' @param digits.csem Integer. Number of decimal places used to round the
#'   CSEM and, by default, the confidence interval limits. Default is 4.
#' @param full.range Logical. If \code{TRUE} (default), CSEM values are
#'   reported for the full score range from 0 to the total number of items.
#'   If \code{FALSE}, CSEM values are reported only for the raw scores
#'   actually observed in the data.
#' @param return.person Logical. If \code{TRUE}, the output includes a
#'   data frame with person-level CSEM values and, if requested, confidence
#'   intervals. Default is \code{FALSE}.
#' @param conf.level Numeric vector or \code{NULL}. Confidence level(s) for
#'   Wilson score intervals for the true number-correct score (e.g.,
#'   \code{0.95} or \code{c(0.90, 0.95)}). If \code{NULL} (default), no
#'   confidence intervals are computed.
#' @param ci.method Character string indicating the method used to compute
#'   confidence intervals for the true score. Currently only
#'   \code{"wilson"} (Wilson score interval) is implemented.
#'
#' @details
#' The strong true score / compound binomial model assumes that, for a
#' given examinee, observed number-correct scores are distributed as a
#' binomial random variable conditional on an underlying true proportion
#' correct, and that the distribution of true scores across examinees
#' follows a beta-like distribution. Under this framework, Lord (1965)
#' derived an estimator of the conditional error variance that adjusts the
#' binomial error model using information about the distribution of item
#' difficulties and total test scores.
#'
#' Let \eqn{n} be the number of items, \eqn{x} a raw number-correct score,
#' \eqn{\hat{\mu}_X} the sample mean of total scores, \eqn{S_X^2} the sample
#' variance of total scores across examinees, and \eqn{S_{Xi}^2} the sample
#' variance of item difficulties (proportion correct across items). Let
#' \eqn{S_{Xi}} denote the average item-level binomial variance
#' \eqn{p_j (1 - p_j)} across items. The strong true score CSEM is based on
#' Lord's adjustment to the binomial error variance:
#' \deqn{
#'   \hat{\sigma}^2_{E|x} =
#'     \frac{x (n - x)}{n - 1}
#'     \left[
#'       1 -
#'       \frac{n (n - 1) S_{Xi}^2}{
#'         \hat{\mu}_X (n - \hat{\mu}_X) - S_X^2 - n S_{Xi}
#'       }
#'     \right],
#' }
#' and the conditional standard error of measurement is
#' \eqn{\sqrt{\hat{\sigma}^2_{E|x}}} for each raw score \eqn{x}.
#'
#' When \code{conf.level} is not \code{NULL}, the function also computes
#' confidence intervals for the true number-correct score \eqn{\tau} using
#' Wilson's (1927) score interval for a binomial proportion. For a given raw
#' score \eqn{x} and test length \eqn{n}, Wilson's interval for the true
#' proportion correct is obtained and then multiplied by \eqn{n} to yield an
#' interval for \eqn{\tau}. For each requested level \eqn{\gamma}, the
#' resulting bounds are reported as \code{lwr.xx} and \code{upr.xx}, where
#' \code{xx} is the confidence level expressed in percent (e.g., 95).
#'
#' The method is intended for tests composed of dichotomously scored items
#' and assumes that the compound binomial model provides a reasonable
#' approximation to the behavior of parallel or closely matched test forms.
#'
#' @return
#' A list of class \code{"csemStrong"} with the following components:
#' \item{score}{A data frame with one row per raw score, containing:
#'   \itemize{
#'     \item \code{raw.score}: Raw number-correct score.
#'     \item \code{n}: Number of examinees attaining that score in \code{data}.
#'     \item \code{csem.strong}: Strong true score CSEM for that score.
#'     \item Additional columns \code{lwr.xx} and \code{upr.xx} for each
#'       requested confidence level, if \code{conf.level} is not \code{NULL}.
#'   }}
#' \item{person}{If \code{return.person = TRUE}, a data frame with one row
#'   per examinee, containing:
#'   \itemize{
#'     \item \code{id}: Row index of the examinee in \code{data}.
#'     \item \code{raw.score}: Raw number-correct score.
#'     \item \code{csem.strong}: CSEM for that examinee's score.
#'     \item The same confidence interval columns as in \code{score}, if
#'       requested.
#'   }
#'   Otherwise, \code{NULL}.}
#' \item{summary}{A data frame summarizing key statistics used in the
#'   computation (number of persons and items, mean and variance of total
#'   scores, variance of item difficulties, and the correction factor).}
#' \item{call}{The matched function call.}
#'
#' @references
#' Lord, F. M. (1965). A strong true score theory, with applications.
#'   \emph{Psychometrika}, 30, 239–270.
#'
#' Tong, Y., & Kolen, M. J. (2018). Conditional standard errors of
#'   measurement. In \emph{Wiley StatsRef: Statistics Reference Online}.
#'   John Wiley & Sons.
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and
#'   statistical inference. \emph{Journal of the American Statistical
#'   Association}, 22(158), 209–212.
#'
#' Agresti, A., & Coull, B. A. (1998). Approximate is better than "exact"
#'   for interval estimation of binomial proportions. \emph{The American
#'   Statistician}, 52(2), 119–126.
#'
#' @examples
#' # Simulated dichotomous test data: 200 persons, 30 items
#' set.seed(123)
#' n_persons <- 200
#' n_items   <- 30
#' # Simple Rasch-like simulation of item difficulties
#' theta <- rnorm(n_persons)
#' b     <- rnorm(n_items, sd = 1)
#' pmat  <- plogis(outer(theta, -b, "+"))
#' data_sim <- matrix(rbinom(n_persons * n_items, size = 1, prob = c(pmat)),
#'                    nrow = n_persons, ncol = n_items)
#'
#' # Strong true score CSEM without confidence intervals
#' res_strong <- csemStrong(data_sim)
#' head(res_strong$score)
#'
#' # Strong true score CSEM with 95% Wilson intervals for the true score
#' res_strong_ci <- csemStrong(data_sim, conf.level = 0.95)
#' head(res_strong_ci$score)
#'
#' @export
csemStrong <- function(
    data,
    na.rm        = TRUE,
    digits.csem  = 4,
    full.range   = TRUE,
    return.person = FALSE,
    conf.level   = NULL,
    ci.method    = c("wilson")
) {
  
  ci.method <- match.arg(ci.method)
  
  ## ---- 1. Comprobaciones básicas y preparación de datos ----
  if (is.data.frame(data)) {
    X <- as.matrix(data)
  } else if (is.matrix(data)) {
    X <- data
  } else {
    stop("`data` must be a matrix or data.frame with persons in rows and items in columns.")
  }
  
  if (!is.numeric(X)) {
    stop("`data` must contain numeric item scores (e.g., 0/1).")
  }
  
  # Manejo de NA
  if (na.rm) {
    cc <- stats::complete.cases(X)
    if (!all(cc)) {
      X <- X[cc, , drop = FALSE]
      warning("Rows with missing values were removed (na.rm = TRUE).")
    }
  } else {
    if (any(!stats::complete.cases(X))) {
      stop("Missing values found in `data`. Set na.rm = TRUE to remove incomplete rows.")
    }
  }
  
  n_persons <- nrow(X)
  n_items   <- ncol(X)
  
  if (n_persons < 2L || n_items < 2L) {
    stop("`data` must have at least 2 persons and 2 items.")
  }
  
  ## ---- 2. Estadísticos globales necesarios para la ecuación (6) ----
  # Puntajes totales por persona
  total_score <- rowSums(X)
  
  # Media y varianza de los puntajes totales (muestra)
  mu_X  <- mean(total_score)
  var_Xp <- stats::var(total_score)
  
  # Dificultades de ítem (proporción correcta)
  p_item <- colMeans(X)
  
  # Varianza de las dificultades de ítem (S^2_Xi)
  var_Xi <- stats::var(p_item)
  
  # Componente medio de varianza binomial por ítem p_j(1 - p_j)
  # Se usa como S_Xi en el denominador de la ecuación (6).
  S_Xi <- mean(p_item * (1 - p_item))
  
  # Numerador y denominador del factor de corrección
  num_K <- n_items * (n_items - 1) * var_Xi
  den_K <- mu_X * (n_items - mu_X) - var_Xp - n_items * S_Xi
  
  if (den_K <= 0) {
    warning("Denominator in the strong true score correction factor is non-positive. ",
            "CSEM strong cannot be computed; returning NA for CSEM.")
    correction_factor <- NA_real_
  } else {
    correction_factor <- 1 - (num_K / den_K)
  }
  
  ## ---- 3. Cálculo de CSEM strong por puntaje ----
  if (full.range) {
    score_vals <- 0:n_items
  } else {
    score_vals <- sort(unique(total_score))
  }
  
  # CSEM binomial base: x(n-x)/(n-1)
  base_var <- score_vals * (n_items - score_vals) / (n_items - 1)
  
  # Varianza strong true score (ec. 6)
  var_strong <- base_var * correction_factor
  
  # Evitar negativos numéricos por redondeo
  var_strong[var_strong < 0] <- NA_real_
  
  csem_strong <- sqrt(var_strong)
  
  ## ---- 4. Intervalos de confianza del true score (Wilson) ----
  ci_cols <- NULL
  if (!is.null(conf.level)) {
    conf.level <- sort(unique(conf.level))
    if (any(conf.level <= 0 | conf.level >= 1)) {
      stop("`conf.level` must contain values strictly between 0 and 1 (e.g., 0.95).")
    }
    
    # Matriz para almacenar los límites (filas: scores, columnas: 2 * length(conf.level))
    ci_mat <- matrix(NA_real_, nrow = length(score_vals),
                     ncol = 2L * length(conf.level))
    
    colnames_ci <- character(2L * length(conf.level))
    
    for (i in seq_along(conf.level)) {
      cl <- conf.level[i]
      z  <- stats::qnorm((1 + cl) / 2)
      
      # Wilson para la proporción p = x/n_items
      # Fórmula clásica de Wilson; luego se multiplica por n_items para obtener τ.
      center <- (score_vals + z^2 / 2) / (n_items + z^2)
      half_w <- z * sqrt((score_vals * (n_items - score_vals) / n_items) + z^2 / 4) /
        (n_items + z^2)
      
      p_low  <- center - half_w
      p_high <- center + half_w
      
      # Asegurar que estén en [0,1]
      p_low[p_low < 0]   <- 0
      p_high[p_high > 1] <- 1
      
      tau_low  <- n_items * p_low
      tau_high <- n_items * p_high
      
      ci_mat[, (2*i - 1)] <- tau_low
      ci_mat[, (2*i)]     <- tau_high
      
      colnames_ci[(2*i - 1)] <- paste0("lwr.", formatC(cl * 100, format = "f", digits = 0))
      colnames_ci[(2*i)]     <- paste0("upr.", formatC(cl * 100, format = "f", digits = 0))
    }
    
    colnames(ci_mat) <- colnames_ci
    ci_cols <- as.data.frame(ci_mat)
  }
  
  ## ---- 5. Tabla por puntaje ----
  freq_tab <- table(total_score)
  freq_vec <- as.numeric(freq_tab[match(score_vals, as.numeric(names(freq_tab)))])
  freq_vec[is.na(freq_vec)] <- 0L
  
  score_df <- data.frame(
    raw.score   = score_vals,
    n           = freq_vec,
    csem.strong = round(csem_strong, digits.csem)
  )
  
  if (!is.null(ci_cols)) {
    # Redondear IC con el mismo número de dígitos que el CSEM (o podrías usar otro)
    ci_cols <- as.data.frame(lapply(ci_cols, round, digits = digits.csem))
    score_df <- cbind(score_df, ci_cols)
  }
  
  ## ---- 6. Tabla por persona (opcional) ----
  person_df <- NULL
  if (isTRUE(return.person)) {
    # Mapear CSEM e IC por puntaje a cada persona
    idx <- match(total_score, score_vals)
    person_df <- data.frame(
      id         = seq_len(n_persons),
      raw.score  = total_score,
      csem.strong = round(csem_strong[idx], digits.csem)
    )
    
    if (!is.null(ci_cols)) {
      for (nm in names(ci_cols)) {
        person_df[[nm]] <- round(ci_cols[[nm]][idx], digits.csem)
      }
    }
  }
  
  ## ---- 7. Resumen de parámetros usados ----
  summary_df <- data.frame(
    parameter = c("n_persons", "n_items", "mu_X", "var_Xp", "var_Xi", "S_Xi", "correction_factor"),
    value     = c(n_persons, n_items, mu_X, var_Xp, var_Xi, S_Xi, correction_factor),
    row.names = NULL
  )
  
  out <- list(
    score   = score_df,
    person  = person_df,
    summary = summary_df,
    call    = match.call()
  )
  
  class(out) <- c("csemStrong", class(out))
  return(out)
}
