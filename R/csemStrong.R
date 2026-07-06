#' Strong True Score (Compound Binomial) Model for CSEM
#'
#' @description
#' Computes conditional standard errors of measurement (CSEM) under Lord's
#' strong true score theory (compound binomial model) for dichotomously
#' scored items. The estimator follows Lord (1965) as summarized by
#' Tong and Kolen (2018). Optionally, the function computes confidence
#' intervals for the true number‑correct score.
#'
#' @param data A numeric matrix or data frame with examinees in rows and
#'   items in columns. For \code{score.type = "dich"} entries must be 0/1;
#'   for \code{"poly"} they range between \code{min.resp} and \code{max.resp}.
#' @param score.type Character: \code{"dich"} (default) or \code{"poly"}.
#' @param nitems Optional integer. Number of items. If \code{NULL} (default),
#'   taken as \code{ncol(data)}.
#' @param min.resp Numeric. Minimum response value per item.
#'   **Required when \code{score.type = "poly"}**.
#' @param max.resp Numeric. Maximum response value per item.
#'   **Required when \code{score.type = "poly"}**.
#' @param bin.score Integer. Number of quantile groups (e.g., 5 for quintiles).
#'   If `NULL` (default), no binning is performed. If provided, the function
#'   returns a data frame \code{binned.CSEM} with average CSEM per quantile group.
#' @param ci Logical. If \code{TRUE}, compute confidence intervals. Default \code{FALSE}.
#' @param ci.method Character: \code{"csem"} (classic) or \code{"wilson"}.
#'   Default \code{"csem"}.
#' @param conf.level Numeric vector of confidence levels (e.g., 0.95).
#'   Default \code{NULL} → 0.95.
#' @param digits Integer for rounding CSEM and confidence limits. Default 3.
#' @param full.range Logical. If \code{TRUE} (default), CSEM values are
#'   reported for the full score range (all possible scores). For dichotomous,
#'   that is 0:nitems; for polytomous, all observed equivalent scores plus
#'   the extremes 0 and nitems. If \code{FALSE}, only scores observed in
#'   \code{data} are reported.
#' @param score.range Optional numeric vector of length 2 (min, max). If provided,
#'   confidence intervals are truncated to this range. If NULL, the observed
#'   score range (or theoretical range for full.range) is used.
#' @param summary Logical. If \code{TRUE}, returns a data frame with key
#'   parameters (mu_X, var_Xp, var_Xi, mean_pq, correction_factor).
#'   Default \code{FALSE}.
#' @param na.rm Logical. If \code{TRUE} (default), rows with missing values
#'   are removed. If \code{FALSE}, stops on missing values.
#'
#' @return A list with components:
#' \item{CSEM}{Data frame with columns: \code{raw.score} (original total),
#'   \code{equiv.score} (if polytomous), \code{n} (frequency), \code{csem.strong},
#'   and confidence limits if requested.}
#' \item{binned.CSEM}{(if \code{bin.score} is provided) Data frame with
#'   quantile groups: \code{group}, \code{range}, \code{n}, \code{mean_score},
#'   \code{CSEM.mean}, and confidence limits if \code{ci = TRUE}.}
#' \item{summary}{(if \code{summary = TRUE}) Data frame with parameters.}
#'
#' @references
#' Lord, F. M. (1965). A strong true score theory, with applications.
#'   \emph{Psychometrika}, 30, 239–270.
#'
#' Tong, Y., & Kolen, M. J. (2018). Conditional standard errors of measurement.
#'   In \emph{Wiley StatsRef: Statistics Reference Online}.
#'
#' @examples
#' # Dichotomous,  strong true score
#' \donttest{
#' library(psychTools)
#' data(ability)
#' data.ability <- ability[complete.cases(ability),]
#'
#' csemStrong(score.type = "dich",
#' data = data.ability,
#' nitems = 16,
#' ci = T,
#' summary = F)
#'
#' # Dichotomous,  Compound Binomial Model (strong true score), more summary and binned score
#' csemStrong(score.type = "dich",
#' data = data.ability,
#' nitems = 16,
#' ci = F,
#' summary = T,
#' bin.score = 5)
#'
#' # Polytomous items,
#' ## Load data
#' library(EFA.dimensions)
#' data("data_RSE")
#'
#' ## Recode negative items
#' data_RSE[c("Q3", "Q5", "Q8", "Q9", "Q10")] <- 5 - data_RSE[c("Q3", "Q5", "Q8", "Q9", "Q10")]
#'
#' csemStrong(score.type = "poly",
#' data = data_RSE,
#' min.resp = 1,
#' max.resp = 4,
#' ci = F,
#' summary = T,
#' full.range = T,
#' bin.score = 4)
#' }
#'
#' @export
csemStrong <- function(data,
                       score.type = c("dich", "poly"),
                       nitems = NULL,
                       min.resp = NULL,
                       max.resp = NULL,
                       bin.score = NULL,
                       ci = FALSE,
                       ci.method = c("csem", "wilson"),
                       conf.level = NULL,
                       digits = 3,
                       full.range = TRUE,
                       score.range = NULL,
                       summary = FALSE,
                       na.rm = TRUE) {

  score.type <- match.arg(score.type)
  ci.method  <- match.arg(ci.method)

  # --- Data preparation ---
  if (is.data.frame(data)) X <- as.matrix(data)
  else if (is.matrix(data)) X <- data
  else stop("`data` must be a matrix or data frame.")

  if (!is.numeric(X)) stop("`data` must contain numeric scores.")

  # Missing values
  if (na.rm) {
    cc <- stats::complete.cases(X)
    if (!all(cc)) {
      X <- X[cc, , drop = FALSE]
      warning("Rows with missing values were removed (na.rm = TRUE).")
    }
  } else {
    if (any(!stats::complete.cases(X)))
      stop("Missing values found. Set na.rm = TRUE to remove incomplete rows.")
  }

  n_persons <- nrow(X)
  n_items_obs <- ncol(X)
  if (n_persons < 2 || n_items_obs < 2)
    stop("`data` must have at least 2 persons and 2 items.")

  if (is.null(nitems)) nitems <- n_items_obs
  else if (nitems != n_items_obs)
    warning("Provided nitems (", nitems, ") differs from ncol(data). Using ncol(data).")
  nitems <- n_items_obs

  # --- Transform to equivalent binomial scale for polytomous ---
  if (score.type == "poly") {
    if (is.null(min.resp) || is.null(max.resp))
      stop("For score.type = 'poly', min.resp and max.resp must be provided.")
    if (min.resp >= max.resp) stop("min.resp must be less than max.resp.")
    X_scaled <- (X - min.resp) / (max.resp - min.resp)
    equiv_total <- rowSums(X_scaled)
    equiv_total <- round(equiv_total, digits = 10)
    raw_total   <- rowSums(X)
    p_item <- colMeans(X_scaled)
    mean_pq <- mean(p_item * (1 - p_item))
    var_Xi <- stats::var(p_item)
    mu_X <- mean(equiv_total)
    var_Xp <- stats::var(equiv_total)
  } else {
    equiv_total <- rowSums(X)
    raw_total   <- equiv_total
    p_item <- colMeans(X)
    mean_pq <- mean(p_item * (1 - p_item))
    var_Xi <- stats::var(p_item)
    mu_X <- mean(equiv_total)
    var_Xp <- stats::var(equiv_total)
  }

  # --- Correction factor (strong true score) ---
  num_K <- nitems * (nitems - 1) * var_Xi
  den_K <- mu_X * (nitems - mu_X) - var_Xp - nitems * mean_pq
  if (den_K <= 0) {
    warning("Denominator non-positive; CSEM strong may be invalid.")
    correction_factor <- NA_real_
  } else {
    correction_factor <- 1 - (num_K / den_K)
  }

  # --- Score values for which CSEM is computed ---
  if (full.range) {
    if (score.type == "dich") {
      score_vals <- 0:nitems
    } else {
      score_vals <- sort(unique(equiv_total))
      if (!any(abs(score_vals - 0) < 1e-9)) score_vals <- c(0, score_vals)
      if (!any(abs(score_vals - nitems) < 1e-9)) score_vals <- c(score_vals, nitems)
      score_vals <- sort(unique(round(score_vals, digits = 10)))
    }
  } else {
    score_vals <- sort(unique(equiv_total))
  }

  # Base binomial variance and strong CSEM
  base_var <- score_vals * (nitems - score_vals) / (nitems - 1)
  base_var[base_var < 0] <- NA_real_
  var_strong <- base_var * correction_factor
  var_strong[var_strong < 0] <- NA_real_
  csem_strong <- sqrt(var_strong)

  # Frequency table
  freq_tab <- table(equiv_total)
  freq_vec <- as.numeric(freq_tab[match(score_vals, as.numeric(names(freq_tab)))])
  freq_vec[is.na(freq_vec)] <- 0L

  # Map raw scores (for polytomous)
  if (score.type == "poly") {
    raw_by_equiv <- tapply(raw_total, equiv_total, unique)
    raw_for_score <- sapply(score_vals, function(v) {
      r <- raw_by_equiv[[as.character(v)]]
      if (is.null(r)) NA else r[1]
    })
  } else {
    raw_for_score <- score_vals
  }

  # Build CSEM data frame
  score_df <- data.frame(
    raw.score = raw_for_score,
    n = freq_vec,
    csem.strong = round(csem_strong, digits)
  )
  if (score.type == "poly") {
    score_df <- cbind(score_df[, 1, drop = FALSE],
                      equiv.score = round(score_vals, digits),
                      score_df[, -1, drop = FALSE])
  }
  # Remove rows where raw.score is NA (extremes added without data)
  score_df <- score_df[!is.na(score_df$raw.score), ]
  if (!full.range) {
    score_df <- score_df[score_df$n > 0, ]
  }

  # --- Confidence intervals (if ci = TRUE) ---
  if (ci) {
    if (is.null(conf.level)) conf.level <- 0.95
    conf.level <- sort(unique(conf.level))
    if (any(conf.level <= 0 | conf.level >= 1))
      stop("conf.level must be between 0 and 1.")
    n <- nitems
    x_vals <- if (score.type == "poly") score_df$equiv.score else score_df$raw.score
    z_vals <- stats::qnorm((1 + conf.level) / 2)
    ci_mat <- matrix(NA, nrow = nrow(score_df), ncol = 2 * length(conf.level))
    colnames_ci <- character(2 * length(conf.level))
    for (i in seq_along(conf.level)) {
      x <- x_vals
      if (ci.method == "csem") {
        low <- x - z_vals[i] * score_df$csem.strong
        high <- x + z_vals[i] * score_df$csem.strong
        # Truncate to theoretical range (or score.range if provided)
        if (!is.null(score.range)) {
          low <- pmax(low, score.range[1])
          high <- pmin(high, score.range[2])
        } else {
          low <- pmax(low, 0)
          high <- pmin(high, n)
        }
      } else { # wilson
        center <- (x + z_vals[i]^2 / 2) / (n + z_vals[i]^2)
        half <- z_vals[i] * sqrt((x * (n - x) / n) + z_vals[i]^2 / 4) / (n + z_vals[i]^2)
        p_low <- center - half
        p_high <- center + half
        p_low[p_low < 0] <- 0
        p_high[p_high > 1] <- 1
        low <- n * p_low
        high <- n * p_high
        if (!is.null(score.range)) {
          low <- pmax(low, score.range[1])
          high <- pmin(high, score.range[2])
        }
      }
      ci_mat[, 2*i - 1] <- low
      ci_mat[, 2*i] <- high
      colnames_ci[2*i - 1] <- paste0("lwr.", formatC(conf.level[i]*100, format = "f", digits = 0))
      colnames_ci[2*i] <- paste0("upr.", formatC(conf.level[i]*100, format = "f", digits = 0))
    }
    colnames(ci_mat) <- colnames_ci
    ci_df <- as.data.frame(ci_mat)
    ci_df <- as.data.frame(lapply(ci_df, round, digits))
    score_df <- cbind(score_df, ci_df)
  }

  # --- Binning (if bin.score is integer) ---
  binned_df <- NULL
  if (!is.null(bin.score)) {
    if (!is.numeric(bin.score) || length(bin.score) != 1 || bin.score < 2)
      stop("bin.score must be an integer >= 2.")
    # Use only observed scores (n > 0) for binning
    obs_df <- score_df[score_df$n > 0, ]
    if (nrow(obs_df) == 0) stop("No observed scores for binning.")
    # Expand to person-level scores (weighted by frequency)
    person_scores <- rep(obs_df$raw.score, times = obs_df$n)
    # Quantile groups on persons
    q <- stats::quantile(person_scores, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
    q <- unique(q)
    groups <- cut(person_scores, breaks = q, include.lowest = TRUE, right = TRUE)
    group_levels <- levels(groups)
    bin_list <- list()
    for (i in seq_along(group_levels)) {
      idx_in_group <- which(groups == group_levels[i])
      scores_in_group <- unique(person_scores[idx_in_group])
      sub_df <- obs_df[obs_df$raw.score %in% scores_in_group, , drop = FALSE]
      if (nrow(sub_df) == 0) next
      csem_mean <- weighted.mean(sub_df$csem.strong, w = sub_df$n)
      range_str <- paste0(min(scores_in_group), "-", max(scores_in_group))
      n_persons <- length(idx_in_group)
      mean_score <- mean(person_scores[idx_in_group])
      bin_list[[i]] <- data.frame(
        group = i,
        range = range_str,
        n = n_persons,
        mean_score = mean_score,
        CSEM.mean = csem_mean,
        stringsAsFactors = FALSE
      )
    }
    binned_df <- do.call(rbind, bin_list)
    # Confidence intervals for binned CSEM
    if (ci && !is.null(binned_df)) {
      cl_main <- conf.level[1]
      z <- stats::qnorm(1 - (1 - cl_main) / 2)
      lwr_bin <- binned_df$mean_score - z * binned_df$CSEM.mean
      upr_bin <- binned_df$mean_score + z * binned_df$CSEM.mean
      if (!is.null(score.range)) {
        lwr_bin <- pmax(lwr_bin, score.range[1])
        upr_bin <- pmin(upr_bin, score.range[2])
      } else {
        lwr_bin <- pmax(lwr_bin, 0)
        upr_bin <- pmin(upr_bin, nitems)
      }
      binned_df$lwr.ci <- round(lwr_bin, digits)
      binned_df$upr.ci <- round(upr_bin, digits)
    }
  }

  # --- Summary parameters ---
  sum_df <- NULL
  if (summary) {
    sum_df <- data.frame(
      parameter = c("n_items", "mu_X", "var_Xp", "var_Xi", "mean_pq", "correction_factor"),
      value = round(c(nitems, mu_X, var_Xp, var_Xi, mean_pq, correction_factor), 5)
    )
  }

  # --- Output ---
  out <- list(CSEM = score_df)
  if (!is.null(binned_df)) out$binned.CSEM <- binned_df
  if (summary) out$summary <- sum_df
  class(out) <- "csemStrong"
  return(out)
}
