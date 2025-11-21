#' Thorndike's method for Conditional Standard Error of Measurement (CSEM)
#'
#' @description
#' Computes the Conditional Standard Error of Measurement (CSEM) using
#' Thorndike's half-test difference method. Optionally smooths the
#' CSEM^2 values with a polynomial regression (Thorndike-style) and
#' provides confidence intervals at the group and score levels.
#'
#' @param half1 Data frame or matrix with items from the first half of the test.
#' @param half2 Data frame or matrix with items from the second half of the test.
#' @param n.groups Integer. Number of groups (quantiles) for total scores.
#'   Ignored if \code{bin.width} is not \code{NULL}.
#' @param bin.width Numeric. If specified, divides total scores into intervals
#'   of this fixed width. Overrides \code{n.groups}.
#' @param smooth Logical. If \code{TRUE}, fits a polynomial regression on CSEM^2
#'   (group-level) and returns smoothed CSEM for each possible total score.
#' @param degree Integer. Polynomial degree used when \code{smooth = TRUE}.
#' @param ci Logical. If \code{TRUE}, computes confidence intervals for
#'   group-level CSEM and, when \code{smooth = TRUE}, for score-level CSEM.
#' @param conf.level Numeric. Confidence level for CIs (default 0.95).
#' @param digits Integer. Number of decimals for CSEM estimates in the output.
#' @param min.score Numeric or integer, optional. Minimum total score at which
#'   the smoothed CSEM curve will be evaluated. If \code{NULL}, the minimum
#'   observed total score (floored) is used.
#' @param max.score Numeric or integer, optional. Maximum total score at which
#'   the smoothed CSEM curve will be evaluated. If \code{NULL}, the maximum
#'   observed total score (ceiled) is used.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{by_group}: data.frame with CSEM by score group (raw),
#'         including CIs if requested.
#'   \item \code{by.score}: data.frame with smoothed CSEM by total score,
#'         including CIs if requested and \code{smooth = TRUE}. \code{NULL} otherwise.
#'   \item \code{params}: data.frame with polynomial coefficients (if \code{smooth = TRUE}),
#'         including standard errors, t-values and p-values. \code{NULL} otherwise.
#'   \item \code{settings}: list with the main arguments used.
#'   \item \code{model}: fitted \code{lm} object for the polynomial on CSEM^2
#'         (if \code{smooth = TRUE}), otherwise \code{NULL}.
#' }
#'
#' @examples
#' # Thorndike CSEM with 10 quantile groups, no smoothing
#' # csemThorndike(half1 = dhalf1, half2 = dhalf2, n.groups = 10, smooth = FALSE)
#'
#' # Thorndike CSEM with smoothing (degree 2) and 95% CIs
#' # csemThorndike(half1 = dhalf1, half2 = dhalf2,
#' #               n.groups = 10, smooth = TRUE,
#' #               ci = TRUE, conf.level = 0.95)
#'
#' # Thorndike CSEM with smoothing over full possible score range
#' # (e.g., dichotomous items, 10 items in total: 0 to 10)
#' # csemThorndike(half1, half2,
#' #               n.groups = 10, smooth = TRUE,
#' #               min.score = 0, max.score = 10)
#'
#' @export
csemThorndike <- function(
    half1, half2,
    n.groups = 10,
    bin.width = NULL,
    smooth = TRUE,
    degree = 2,
    ci = FALSE,
    conf.level = 0.95,
    digits = 3,
    min.score = NULL,
    max.score = NULL
) {
  # --- 0) Checks básicos ---
  if (missing(half1) || missing(half2)) {
    stop("Both 'half1' and 'half2' must be provided.")
  }
  half1 <- as.data.frame(half1)
  half2 <- as.data.frame(half2)

  if (nrow(half1) != nrow(half2)) {
    stop("'half1' and 'half2' must have the same number of rows (subjects).")
  }

  if (!is.null(bin.width) && !is.numeric(bin.width)) {
    stop("'bin.width' must be numeric if provided.")
  }
  if (!is.numeric(n.groups) || n.groups < 2) {
    stop("'n.groups' must be an integer >= 2.")
  }
  if (!is.numeric(degree) || degree < 1) {
    stop("'degree' must be a positive integer.")
  }
  if (!is.logical(smooth) || length(smooth) != 1) {
    stop("'smooth' must be a single logical value.")
  }
  if (!is.logical(ci) || length(ci) != 1) {
    stop("'ci' must be a single logical value.")
  }
  if (conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a number between 0 and 1.")
  }
  if (!is.null(min.score) && !is.numeric(min.score)) {
    stop("'min.score' must be numeric if provided.")
  }
  if (!is.null(max.score) && !is.numeric(max.score)) {
    stop("'max.score' must be numeric if provided.")
  }

  # --- 1) Totales y diferencias ---
  total1 <- rowSums(half1, na.rm = TRUE)
  total2 <- rowSums(half2, na.rm = TRUE)
  total  <- total1 + total2
  diffD  <- total1 - total2

  df <- data.frame(total = total, diff = diffD)

  # --- 2) Agrupación de puntajes ---
  used_scheme <- NULL

  if (!is.null(bin.width)) {
    if (!missing(n.groups)) {
      warning("'n.groups' is ignored because 'bin.width' was specified.")
    }
    minT <- min(df$total, na.rm = TRUE)
    maxT <- max(df$total, na.rm = TRUE)
    breaks <- seq(minT, maxT + bin.width, by = bin.width)
    if (length(breaks) < 2) {
      stop("Not enough range in total scores for the specified 'bin.width'.")
    }
    df$group <- cut(df$total, breaks = breaks, include.lowest = TRUE, right = FALSE)
    used_scheme <- paste0("fixed width (", bin.width, ")")
  } else {
    probs <- seq(0, 1, length.out = n.groups + 1)
    qtls  <- unique(stats::quantile(df$total, probs = probs, na.rm = TRUE, type = 7))
    if (length(qtls) < 2) {
      stop("Not enough distinct total scores to form the requested quantile groups.")
    }
    df$group <- cut(df$total, breaks = qtls, include.lowest = TRUE, right = TRUE)
    used_scheme <- paste0("quantiles (n.groups = ", n.groups, ")")
  }

  # --- 3) Estadísticos por grupo ---
  levs <- levels(df$group)
  if (is.null(levs)) {
    stop("No valid groups were formed. Check 'n.groups' or 'bin.width'.")
  }

  Nbin <- tapply(df$diff,  df$group, function(x) sum(!is.na(x)))
  mTot <- tapply(df$total, df$group, mean, na.rm = TRUE)
  varD <- tapply(df$diff,  df$group, stats::var, na.rm = TRUE)

  Nbin <- Nbin[levs]
  mTot <- mTot[levs]
  varD <- varD[levs]

  # CSEM crudo por grupo: sqrt(Var(D)/4) = sd(D)/2
  CSEM_raw <- rep(NA_real_, length(levs))
  okN <- !is.na(Nbin) & (Nbin > 1) & !is.na(varD)
  CSEM_raw[okN] <- sqrt(pmax(varD[okN], 0) / 4)

  # --- 4) IC por grupo: para el puntaje verdadero ---
  lwr_group <- upp_group <- rep(NA_real_, length(levs))
  if (ci) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    # IC para el puntaje verdadero del grupo: Mean_Total ± z * CSEM
    lwr_group[okN] <- mTot[okN] - z * CSEM_raw[okN]
    upp_group[okN] <- mTot[okN] + z * CSEM_raw[okN]
  }

  by_group <- data.frame(
    Group      = factor(levs, levels = levs),
    N          = as.integer(Nbin),
    Mean_Total = round(as.numeric(mTot), 2),
    CSEM       = round(as.numeric(CSEM_raw), digits = digits),
    stringsAsFactors = FALSE
  )

  if (ci) {
    by_group$lwr.ci <- round(lwr_group, digits = digits)
    by_group$upp.ci <- round(upp_group, digits = digits)
  }

  # --- 5) Suavizado polinómico sobre CSEM^2 (si smooth = TRUE) ---
  by_score  <- NULL
  params_df <- NULL
  fit_model <- NULL

  if (smooth) {
    x  <- by_group$Mean_Total
    y2 <- CSEM_raw^2

    ok_fit <- is.finite(x) & is.finite(y2) & !is.na(x) & !is.na(y2)
    x_fit  <- x[ok_fit]
    y2_fit <- y2[ok_fit]

    if (length(x_fit) < (degree + 1)) {
      stop(
        "Not enough groups with valid CSEM to fit a polynomial of degree = ",
        degree, ".\nIncrease 'n.groups' or reduce 'degree', or set smooth = FALSE."
      )
    }

    Xmat <- data.frame(x = x_fit)
    if (degree >= 2) {
      for (d in 2:degree) {
        Xmat[[paste0("x", d)]] <- x_fit^d
      }
    }

    form <- stats::as.formula(
      paste("y2_fit ~", paste(colnames(Xmat), collapse = " + "))
    )
    dat_fit <- cbind(y2_fit = y2_fit, Xmat)
    fit_model <- stats::lm(form, data = dat_fit)

    # --- rango de scores para la curva CSEM ---
    if (is.null(min.score)) {
      min_s <- floor(min(total, na.rm = TRUE))
    } else {
      min_s <- min.score
    }
    if (is.null(max.score)) {
      max_s <- ceiling(max(total, na.rm = TRUE))
    } else {
      max_s <- max.score
    }
    if (max_s <= min_s) {
      stop("'max.score' must be greater than 'min.score'.")
    }

    score_seq <- seq(from = min_s, to = max_s, by = 1)

    newX <- data.frame(x = score_seq)
    if (degree >= 2) {
      for (d in 2:degree) {
        newX[[paste0("x", d)]] <- score_seq^d
      }
    }

    # Predicción puntual de CSEM^2 y CSEM
    var_hat  <- pmax(stats::predict(fit_model, newdata = newX), 0)
    CSEM_hat <- sqrt(var_hat)

    if (ci) {
      z <- stats::qnorm(1 - (1 - conf.level) / 2)
      # IC para el puntaje verdadero individual: Score ± z * CSEM_hat
      lwr_score <- score_seq - z * CSEM_hat
      upp_score <- score_seq + z * CSEM_hat

      by_score <- data.frame(
        Score  = score_seq,
        CSEM   = round(CSEM_hat, digits = digits),
        lwr.ci = round(lwr_score, digits = digits),
        upp.ci = round(upp_score, digits = digits)
      )
    } else {
      by_score <- data.frame(
        Score = score_seq,
        CSEM  = round(CSEM_hat, digits = digits)
      )
    }

    sum_fit <- summary(fit_model)
    coefs   <- sum_fit$coefficients
    params_df <- data.frame(
      term     = rownames(coefs),
      estimate = coefs[, "Estimate"],
      se       = coefs[, "Std. Error"],
      t.value  = coefs[, "t value"],
      p.value  = coefs[, "Pr(>|t|)"],
      row.names = NULL
    )
  }

  # --- 6) Settings y mensaje ---
  settings <- list(
    n.groups   = n.groups,
    bin.width  = bin.width,
    smooth     = smooth,
    degree     = degree,
    ci         = ci,
    conf.level = conf.level,
    digits     = digits,
    grouping   = used_scheme,
    min.score  = min.score,
    max.score  = max.score
  )

  message("csemThorndike: CSEM computed. Grouping: ", used_scheme,
          "; Smoothing: ", if (smooth) paste0("polynomial (degree = ", degree, ")") else "none",
          "; CI (true score): ", if (ci) paste0("TRUE (", conf.level * 100, "%)") else "FALSE")

  list(
    by_group = by_group,
    by_score = by_score,
    params   = params_df,
    settings = settings,
    model    = fit_model
  )
}
