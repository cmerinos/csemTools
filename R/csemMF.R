#' Mollenkopf-Feldt Conditional Standard Error of Measurement (MF-CSEM)
#'
#' @description
#' Computes the Conditional Standard Error of Measurement (CSEM) using the
#' Mollenkopf-Feldt (MF) procedure. This method uses polynomial regression
#' to smooth the squared adjusted differences between two test halves, and
#' optionally centers confidence intervals on observed scores or estimated
#' true scores (ETS). The regression can be forced through the origin
#' (constant = 0) as recommended in the literature.
#'
#' @param half1 A data frame or matrix with items from the first half of the test.
#' @param half2 A data frame or matrix with items from the second half of the test.
#'   Must have the same number of rows as \code{half1}.
#' @param degree Integer. Degree of the polynomial regression (default = 2).
#' @param constant Character. Whether to force the polynomial regression through
#'   the origin: \code{"0"} (default) or \code{"free"} (allow intercept).
#' @param center.on Character. Where to center confidence intervals:
#'   \code{"observed"} (default) centers on the observed score;
#'   \code{"true"} centers on the estimated true score (ETS) using reliability.
#' @param reliability Optional numeric value. Reliability coefficient (e.g., alpha)
#'   used to compute ETS. If \code{NULL} and \code{center.on = "true"}, reliability
#'   is estimated internally using Cronbach's alpha on the full test.
#' @param ci Logical. If \code{TRUE}, compute confidence intervals for the true score.
#'   Default \code{FALSE}.
#' @param conf.level Numeric. Confidence level for intervals (default 0.95).
#' @param digits.csem Integer. Number of decimal places for CSEM and confidence limits.
#'   Default 3.
#' @param full.range Logical. If \code{TRUE}, the smoothed CSEM curve is evaluated
#'   for every integer score from 0 to the total number of items (theoretical range).
#'   Requires \code{smooth = TRUE}? This function always smooths. Default \code{FALSE}.
#' @param score.range Numeric vector of length 2, optional. Overrides the score range
#'   for evaluation (e.g., \code{c(0, 40)}). Only used when \code{full.range = FALSE}.
#' @param table.person Logical. If \code{TRUE}, returns a data frame with CSEM values
#'   for each individual in the original data. Default \code{FALSE}.
#' @param summary Logical. If \code{TRUE}, returns additional information about the
#'   polynomial fit (coefficients, R-squared, AIC, BIC, etc.). Default \code{FALSE}.
#' @param na.rm Logical. If \code{TRUE} (default), removes rows with any missing
#'   values in \code{half1} or \code{half2}. If \code{FALSE}, stops on missing values.
#'
#' @details
#' The Mollenkopf-Feldt procedure (Mollenkopf, 1949; Feldt et al., 1985) consists of:
#' \enumerate{
#'   \item Divide the test into two parallel halves (same number of items, similar difficulty).
#'   \item For each examinee, compute the total score \eqn{X = X_1 + X_2} and the
#'         adjusted squared difference \eqn{Y = ((X_1 - X_2) - (\bar{X}_1 - \bar{X}_2))^2}.
#'   \item Fit a polynomial regression of \eqn{Y} on \eqn{X} of degree \code{degree},
#'         optionally forcing the intercept to 0.
#'   \item The predicted values \eqn{\hat{Y}} are estimates of the conditional error variance.
#'         The CSEM for a score \eqn{s} is \eqn{\sqrt{\max(\hat{Y}(s), 0)}}.
#'   \item If \code{center.on = "true"}, confidence intervals for the true score are
#'         centered on the Estimated True Score (ETS) using the reliability coefficient:
#'         \eqn{ETS = r_{xx} \cdot X + (1 - r_{xx}) \cdot \bar{X}}.
#'   \item Otherwise, intervals are centered on the observed score.
#' }
#'
#' @return
#' A list with components:
#' \item{table_score}{Data frame with columns: \code{Score}, \code{CSEM}, and
#'   confidence limits (if \code{ci = TRUE}).}
#' \item{table_person}{(if \code{table.person = TRUE}) Data frame with:
#'   \code{id}, \code{total_score}, \code{CSEM}, and confidence limits (if \code{ci = TRUE}).}
#' \item{summary}{(if \code{summary = TRUE}) A list with:
#'   \code{general} (parameters of the call) and \code{polynomial} (coefficients,
#'   R-squared, AIC, BIC, residual SE, etc.).}
#' \item{reliability}{(if \code{center.on = "true"} and \code{reliability} was
#'   computed internally) The reliability coefficient used for ETS.}
#'
#' @references
#' Mollenkopf, W. G. (1949). Variation of the standard error of measurement of scores.
#'   \emph{Psychometrika}, 14(3), 189–229.
#' Feldt, L. S., Steffen, M., & Gupta, N. C. (1985). A comparison of five methods
#'   for estimating the standard error of measurement at specific score levels.
#'   \emph{Applied Psychological Measurement}, 9(4), 351–361.
#'
#' @examples
#' \dontrun{
#' # Simulate data: 200 persons, 10 items per half (total 20 items)
#' set.seed(123)
#' half1 <- matrix(rbinom(200*10, 1, 0.6), ncol=10)
#' half2 <- matrix(rbinom(200*10, 1, 0.6), ncol=10)
#'
#' # Basic MF-CSEM with default settings
#' res <- csemMF(half1, half2, degree = 2, constant = "0")
#' head(res$table_score)
#'
#' # With confidence intervals centered on observed scores
#' res2 <- csemMF(half1, half2, ci = TRUE, conf.level = 0.95,
#'                center.on = "observed", full.range = TRUE)
#'
#' # With confidence intervals centered on true scores (reliability estimated)
#' res3 <- csemMF(half1, half2, ci = TRUE, center.on = "true", full.range = TRUE)
#'
#' # Return person-level CSEM and polynomial summary
#' res4 <- csemMF(half1, half2, table.person = TRUE, summary = TRUE)
#' }
#'
#' @export
csemMF <- function(half1, half2,
                   degree = 2,
                   constant = c("0", "free"),
                   center.on = c("observed", "true"),
                   reliability = NULL,
                   ci = FALSE,
                   conf.level = 0.95,
                   digits.csem = 3,
                   full.range = FALSE,
                   score.range = NULL,
                   table.person = FALSE,
                   summary = FALSE,
                   na.rm = TRUE) {

  # --- Argument matching and validation ---
  constant <- match.arg(constant)
  center.on <- match.arg(center.on)

  if (missing(half1) || missing(half2))
    stop("Both 'half1' and 'half2' must be provided.")
  half1 <- as.data.frame(half1)
  half2 <- as.data.frame(half2)

  if (nrow(half1) != nrow(half2))
    stop("'half1' and 'half2' must have the same number of rows.")
  if (ncol(half1) != ncol(half2))
    stop("Both halves must have the same number of items.")

  # Handle missing values
  if (na.rm) {
    cc1 <- stats::complete.cases(half1)
    cc2 <- stats::complete.cases(half2)
    ok <- cc1 & cc2
    if (!all(ok)) {
      half1 <- half1[ok, , drop = FALSE]
      half2 <- half2[ok, , drop = FALSE]
      warning("Rows with missing values were removed (na.rm = TRUE).")
    }
  } else {
    if (anyNA(half1) || anyNA(half2))
      stop("Missing values found. Set na.rm = TRUE to remove incomplete rows.")
  }

  n_persons <- nrow(half1)
  n_items_total <- ncol(half1) + ncol(half2)
  if (n_persons < 3) stop("At least 3 persons are required for regression.")
  if (degree < 1) stop("degree must be at least 1.")
  if (degree > 5) warning("Degree > 5 may lead to overfitting. Use with caution.")

  # --- Compute total scores and adjusted squared differences ---
  X1 <- rowSums(half1)
  X2 <- rowSums(half2)
  total_score <- X1 + X2
  diff_score <- X1 - X2
  mean_diff <- mean(diff_score)
  Y <- (diff_score - mean_diff)^2   # adjusted squared difference

  # --- Polynomial regression ---
  # Create polynomial terms (raw, not orthogonal for interpretability)
  poly_terms <- poly(total_score, degree, raw = TRUE)
  colnames(poly_terms) <- paste0("poly", 1:degree)
  df_model <- data.frame(Y = Y, poly_terms)

  if (constant == "0") {
    # Force through origin: remove intercept
    form <- as.formula(paste("Y ~ ", paste(colnames(poly_terms), collapse = " + "), " - 1"))
    fit <- lm(form, data = df_model)
  } else {
    form <- as.formula(paste("Y ~ ", paste(colnames(poly_terms), collapse = " + ")))
    fit <- lm(form, data = df_model)
  }

  # --- Determine evaluation range for scores ---
  if (full.range) {
    eval_min <- 0
    eval_max <- n_items_total
  } else {
    if (!is.null(score.range)) {
      if (!is.numeric(score.range) || length(score.range) != 2 || score.range[1] >= score.range[2])
        stop("score.range must be numeric c(min, max) with min < max.")
      eval_min <- score.range[1]
      eval_max <- score.range[2]
    } else {
      eval_min <- floor(min(total_score))
      eval_max <- ceiling(max(total_score))
    }
  }
  if (eval_max <= eval_min) stop("Invalid evaluation range.")
  score_seq <- seq(from = eval_min, to = eval_max, by = 1)

  # --- Predict Y_hat and compute CSEM ---
  # Build newdata matrix with same polynomial terms
  newdata <- as.data.frame(poly(score_seq, degree, raw = TRUE))
  colnames(newdata) <- colnames(poly_terms)
  Y_hat <- predict(fit, newdata = newdata)
  Y_hat <- pmax(Y_hat, 0)   # variance cannot be negative
  CSEM_seq <- sqrt(Y_hat)

  # --- Confidence intervals ---
  if (ci) {
    if (conf.level <= 0 || conf.level >= 1) stop("conf.level must be between 0 and 1.")
    z <- stats::qnorm(1 - (1 - conf.level) / 2)

    if (center.on == "observed") {
      center_seq <- score_seq
    } else { # true score centering
      # Compute reliability if not provided
      if (is.null(reliability)) {
        # Use Cronbach's alpha on the full test (both halves combined)
        full_test <- cbind(half1, half2)
        # Alpha formula: (k/(k-1)) * (1 - sum(var_items)/var_total)
        k <- ncol(full_test)
        item_vars <- apply(full_test, 2, var, na.rm = TRUE)
        total_var <- var(rowSums(full_test))
        reliability <- (k/(k-1)) * (1 - sum(item_vars)/total_var)
        # Ensure within bounds
        reliability <- max(0, min(1, reliability))
      }
      mean_total <- mean(total_score)
      center_seq <- reliability * score_seq + (1 - reliability) * mean_total
    }

    lwr <- center_seq - z * CSEM_seq
    upr <- center_seq + z * CSEM_seq
    # Truncate to possible score range (for observed) or to theoretical range? Keep general.
    # For observed centering, lower bound can be below 0; we'll keep as is.
    lwr <- round(lwr, digits.csem)
    upr <- round(upr, digits.csem)
  }

  # --- Build score-level table ---
  score_df <- data.frame(Score = score_seq,
                         CSEM = round(CSEM_seq, digits.csem))
  if (ci) {
    score_df$lwr.ci <- lwr
    score_df$upr.ci <- upr
  }

  # --- Person-level table (if requested) ---
  person_df <- NULL
  if (table.person) {
    # Predict Y_hat for each person's total score
    # Create polynomial matrix for original total scores
    X_poly <- as.data.frame(poly(total_score, degree, raw = TRUE))
    colnames(X_poly) <- colnames(poly_terms)
    Y_hat_person <- predict(fit, newdata = X_poly)
    Y_hat_person <- pmax(Y_hat_person, 0)
    CSEM_person <- sqrt(Y_hat_person)

    person_df <- data.frame(id = seq_len(n_persons),
                            total_score = total_score,
                            CSEM = round(CSEM_person, digits.csem))
    if (ci) {
      if (center.on == "observed") {
        center_person <- total_score
      } else {
        center_person <- reliability * total_score + (1 - reliability) * mean(total_score)
      }
      lwr_p <- center_person - z * CSEM_person
      upr_p <- center_person + z * CSEM_person
      person_df$lwr.ci <- round(lwr_p, digits.csem)
      person_df$upr.ci <- round(upr_p, digits.csem)
    }
  }

  # --- Summary information (if requested) ---
  summary_out <- NULL
  if (summary) {
    # General info
    gen_info <- data.frame(
      parameter = c("degree", "constant", "center.on", "n_persons", "n_items_total",
                    "eval_min", "eval_max", "ci", "conf.level", "full.range"),
      value = c(degree, constant, center.on, n_persons, n_items_total,
                eval_min, eval_max, ci, conf.level, full.range),
      stringsAsFactors = FALSE
    )
    # Polynomial fit statistics
    s <- summary(fit)
    coef_tab <- s$coefficients
    if (constant == "0") rownames(coef_tab) <- paste0("beta", 1:degree)
    poly_info <- list(
      degree = degree,
      coefficients = coef_tab,
      r.squared = s$r.squared,
      adj.r.squared = s$adj.r.squared,
      residual.se = s$sigma,
      AIC = AIC(fit),
      BIC = BIC(fit),
      logLik = as.numeric(logLik(fit)),
      deviance = deviance(fit),
      df.residual = df.residual(fit)
    )
    summary_out <- list(general = gen_info, polynomial = poly_info)
  }

  # --- Output ---
  out <- list(table_score = score_df)
  if (table.person) out$table_person <- person_df
  if (summary) out$summary <- summary_out
  if (center.on == "true" && is.null(reliability)) out$reliability <- reliability

  return(out)
}
