#' Bootstrap Conditional Standard Error of Measurement (CSEM) based on Item Sampling
#'
#' Estimates the conditional standard error of measurement for each total score
#' using a non‑parametric bootstrap of items (tasks). This method simulates
#' parallel test forms by resampling items with replacement and computes the
#' standard deviation of the resulting bootstrap scores for each person.
#' The individual standard errors are then averaged within each observed total
#' score to obtain the conditional SEM. The procedure follows the logic of
#' Colton, Gao, & Kolen (1996) but is adapted to continuous total scores
#' (sum of Likert items) rather than discrete performance levels.
#'
#' @param data A data frame or matrix containing the item responses.
#'   Rows represent persons, columns represent items. Items can be
#'   dichotomous or polytomous (e.g., Likert scales).
#' @param B Integer. Number of bootstrap replications (default = 500).
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'   If `cores = 1` (default), the bootstrap runs sequentially with a progress bar.
#'   If `cores > 1`, a parallel cluster is created and the bootstrap is distributed
#'   across the cores. A progress bar is still shown (using `pbapply`).
#' @param bin.score Integer or NULL. If an integer (e.g., 5), persons are grouped
#'   into that many quantile groups based on their total score, and the average
#'   CSEM within each group is reported. If `NULL` (default), CSEM is reported
#'   for each distinct observed total score (with at least 2 persons).
#' @param smooth Logical. If `TRUE`, a polynomial regression of the squared CSEM
#'   on total score is fitted, and the smoothed CSEM (square root of the fitted
#'   values) is reported. This stabilizes estimates at score levels with few
#'   observations. Default = `FALSE`.
#' @param degree Integer. Polynomial degree used when `smooth = TRUE`. Default = 2.
#' @param full.range Logical. If `TRUE` and `smooth = TRUE`, the smoothed CSEM
#'   is evaluated for every integer score in `score.range`. If `smooth = FALSE`
#'   and `full.range = TRUE`, missing values are filled by carrying forward the
#'   last observed CSEM (`na.locf`). Requires `score.range` to be provided.
#'   Default = `FALSE`.
#' @param ci Logical. If `TRUE`, confidence intervals for the true score are
#'   computed as `score ± z * CSEM`, assuming normally distributed errors.
#'   The intervals are truncated to `score.range` (if supplied) or to the
#'   observed range. Default = `FALSE`.
#' @param conf.level Numeric. Confidence level for the intervals (default = 0.95).
#' @param digits Integer. Number of decimal places for rounding the output.
#'   Default = 3.
#' @param score.range Numeric vector of length 2 (min, max). Required when
#'   `full.range = TRUE`. Defines the theoretical score range (e.g., c(10,50)).
#'   Also used to truncate confidence intervals (if `ci = TRUE`). If `NULL`,
#'   the observed minimum and maximum are used for truncation.
#' @param na.rm Logical. If `TRUE` (default), rows (persons) with any missing
#'   item responses are removed before analysis.
#'
#' @return A list with up to two elements:
#' \describe{
#'   \item{CSEM}{A data frame with columns:
#'     \itemize{
#'       \item `score` : total score (integer).
#'       \item `n` : number of persons with that score.
#'       \item `CSEM` (or `CSEM.smooth`) : conditional standard error of measurement.
#'       \item `lwr.ci`, `upr.ci` : confidence limits (if `ci = TRUE`).
#'     }
#'   }
#'   \item{binned.CSEM}{If `bin.score` is provided, a data frame with:
#'     \itemize{
#'       \item `group` : quantile group number.
#'       \item `range` : score range of the group.
#'       \item `n` : number of persons in the group.
#'       \item `mean_score` : average total score in the group.
#'       \item `CSEM.mean` : average CSEM within the group.
#'       \item `lwr.ci`, `upr.ci` : confidence limits (if `ci = TRUE`).
#'     }
#'   }
#' }
#'
#' @details
#' \subsection{Methodological flow (adapted from Colton et al., 1996)}{
#' 1. The observed total score is computed for each person as the sum of item responses.
#' 2. For each bootstrap replication `b = 1 … B`:
#'    \itemize{
#'      \item A bootstrap sample of **items** (columns) is drawn with replacement,
#'            preserving the original number of items.
#'      \item For every person, the total score on the bootstrap item set is computed.
#'    }
#' 3. After all replications, each person has `B` bootstrap total scores.
#'    The standard deviation of these `B` scores is the **individual standard error
#'    of measurement** (\eqn{\hat{\sigma}(\Delta_p)}).
#' 4. Persons are grouped by their **observed total score** (from the original data).
#'    Within each score level, the individual error variances are averaged, and the
#'    square root is taken: \eqn{CSEM(x) = \sqrt{ \frac{1}{n_x} \sum_{p: X_p=x} \hat{\sigma}^2(\Delta_p) }}.
#' 5. Optionally, a polynomial smoothing is applied to the squared CSEM values to
#'    reduce sampling variability at score levels with few persons.
#' 6. If `bin.score` is provided, the persons are divided into quantile groups
#'    (based on the original total scores) and the average CSEM inside each group
#'    is reported.
#' }
#'
#' \subsection{Why bootstrap items (not persons)?}{
#' The standard error of measurement refers to the variability of a person’s score
#' across hypothetically parallel test forms. Resampling persons would simulate
#' sampling from the population, not parallel forms. Resampling **items** creates
#' different versions of the test, which is the correct approach to estimate the
#' conditional error variance. This is the essence of the method used by
#' Colton et al. (1996) for performance level scores, here extended to continuous
#' total scores.
#' }
#'
#' \subsection{Parallel processing and progress bar}{
#' When `cores = 1`, the bootstrap runs sequentially and a progress bar is shown
#' via `pbapply::pblapply`. When `cores > 1`, a cluster is created using
#' `parallel::makeCluster`, the bootstrap is distributed with `pbapply::pblapply`
#' (which automatically displays a progress bar for parallel tasks), and the
#' cluster is stopped afterwards. The number of cores should not exceed the
#' available CPU cores.
#' }
#'
#' \subsection{Handling of missing data and small groups}{
#' Persons with any missing item responses are removed if `na.rm = TRUE`.
#' Score levels with fewer than 2 persons receive a missing CSEM (NA) in the raw
#' output; they are omitted from smoothing and from the `full.range` expansion
#' (except when `full.range = TRUE` and `smooth = FALSE`, where NAs are filled
#' by last observation carried forward).
#' }
#'
#' @references
#' Colton, D. A., Gao, X., & Kolen, M. J. (1997). Assessing the reliability of
#' performance level scores using bootstrapping. ACT Research Report Series,
#' 97-3. Iowa City, 1A: ACT.
#'
#' @examples
#' \donttest{
#' # Simulate data: 500 persons, 10 Likert items (1-5)
#' set.seed(123)
#' resp <- matrix(sample(1:5, 500*10, replace = TRUE), nrow = 500)
#'
#' # Basic usage (sequential, no smoothing)
#' result <- csemBoots(resp, B = 100, cores = 1)
#' head(result$CSEM)
#'
#' # With smoothing, full range, and confidence intervals
#' result2 <- csemBoots(resp, B = 200, smooth = TRUE, full.range = TRUE,
#'                      score.range = c(10, 50), ci = TRUE, cores = 2)
#' print(result2$CSEM)
#'
#' # Quantile groups (quintiles)
#' result3 <- csemBoots(resp, B = 100, bin.score = 5)
#' print(result3$binned.CSEM)
#' }
#'
#' @importFrom pbapply pblapply
#' @importFrom parallel makeCluster stopCluster
#' @export
csemBoots <- function(data,
                      B = 500,
                      cores = 1,
                      bin.score = NULL,
                      smooth = FALSE,
                      degree = 2,
                      full.range = FALSE,
                      ci = FALSE,
                      conf.level = 0.95,
                      digits = 3,
                      score.range = NULL,
                      na.rm = TRUE) {

  # ---------- 1. Validations and preparation ----------
  if (!is.data.frame(data) && !is.matrix(data))
    stop("'data' must be a data frame or matrix.")
  data <- as.data.frame(data)

  if (na.rm) data <- stats::na.omit(data)
  if (anyNA(data)) stop("Missing values present. Set na.rm = TRUE to remove them.")

  n_persons <- nrow(data)
  n_items   <- ncol(data)
  if (n_persons < 2) stop("At least 2 persons required.")
  if (n_items < 2) stop("At least 2 items required.")

  # Observed total scores
  total <- rowSums(data)

  # Theoretical score range for truncation and full.range
  if (!is.null(score.range)) {
    if (!is.numeric(score.range) || length(score.range) != 2)
      stop("score.range must be a numeric vector of length 2.")
    score_min_teo <- score.range[1]
    score_max_teo <- score.range[2]
  } else {
    score_min_teo <- min(total)
    score_max_teo <- max(total)
  }

  # ---------- 2. Bootstrap of items (parallel or sequential) ----------
  # We will create a matrix of bootstrap scores: persons x B
  # But to avoid memory overhead, we can compute CSEM per person on the fly?
  # Better to compute all bootstrap scores first (n_persons x B).
  # Using a list of length B, each element a vector of length n_persons.

  # Function to generate one bootstrap replication (returns a vector of total scores)
  boot_one <- function(b, data, n_items) {
    item_indices <- sample(1:n_items, size = n_items, replace = TRUE)
    boot_data <- data[, item_indices, drop = FALSE]
    boot_totals <- rowSums(boot_data)
    return(boot_totals)
  }

  # Set up parallel backend if cores > 1
  if (cores > 1) {
    if (!requireNamespace("parallel", quietly = TRUE))
      stop("Package 'parallel' is required for cores > 1.")
    cl <- parallel::makeCluster(cores)
    # Export needed objects to cluster (though boot_one sees data from closure)
    parallel::clusterExport(cl, c("data", "n_items"), envir = environment())
    # Use pbapply::pblapply with cluster
    boot_list <- pbapply::pblapply(1:B, function(b) boot_one(b, data, n_items),
                                   cl = cl)
    parallel::stopCluster(cl)
  } else {
    # Sequential with progress bar
    boot_list <- pbapply::pblapply(1:B, function(b) boot_one(b, data, n_items))
  }

  # Convert list to matrix: rows = persons, columns = bootstrap replications
  boot_mat <- do.call(cbind, boot_list)  # n_persons x B

  # ---------- 3. Individual standard errors (SD across bootstrap replications) ----------
  csem_person <- apply(boot_mat, 1, sd, na.rm = TRUE)

  # ---------- 4. Aggregate by observed total score ----------
  # Data frame with observed score and individual csem
  df <- data.frame(score = total, csem = csem_person)

  # Raw CSEM per observed score (at least 2 persons)
  raw <- aggregate(csem ~ score, data = df, FUN = function(x) {
    c(CSEM = sqrt(mean(x^2)), n = length(x))})
  raw <- do.call(data.frame, raw)  # flatten
  names(raw) <- c("score", "CSEM", "n")
  raw$CSEM <- round(raw$CSEM, digits)
  raw <- raw[raw$n >= 2, ]

  # Remove scores with n < 2 (CSEM already NA? Actually aggregate with mean gives NaN if n<2? But we only have n>=2)
  # In our aggregate, every score present has at least one person, but the csem is computed per person.
  # If a score has only one person, the csem_person exists but the mean is that value, and n=1.
  # We should keep only scores with n >= 2 to have a meaningful estimate.
  raw <- raw[raw$n >= 2, ]

  # ---------- 5. Smoothing (if requested) ----------
  if (smooth) {
    if (nrow(raw) < degree + 1)
      stop("Not enough distinct score levels (", nrow(raw), ") to fit polynomial of degree ", degree)
    fit <- lm(CSEM^2 ~ poly(score, degree, raw = TRUE), data = raw)
    if (full.range) {
      if (is.null(score.range))
        stop("full.range = TRUE requires 'score.range' to be provided.")
      all_scores <- seq(score_min_teo, score_max_teo, by = 1)
      pred_var <- predict(fit, newdata = data.frame(score = all_scores))
      pred_var <- pmax(pred_var, 0)
      csem_smooth <- sqrt(pred_var)
      n_all <- sapply(all_scores, function(s) sum(total == s))
      result_df <- data.frame(score = all_scores, n = n_all,
                              CSEM.smooth = round(csem_smooth, digits))
    } else {
      pred_var <- predict(fit, newdata = data.frame(score = raw$score))
      pred_var <- pmax(pred_var, 0)
      csem_smooth <- sqrt(pred_var)
      result_df <- data.frame(score = raw$score, n = raw$n,
                              CSEM.smooth = round(csem_smooth, digits))
    }
    csem_col <- "CSEM.smooth"
  } else {
    # No smoothing: raw CSEM
    if (full.range) {
      if (is.null(score.range))
        stop("full.range = TRUE requires 'score.range' to be provided.")
      all_scores <- seq(score_min_teo, score_max_teo, by = 1)
      n_all <- sapply(all_scores, function(s) sum(total == s))
      csem_map <- setNames(raw$CSEM, raw$score)
      csem_full <- csem_map[as.character(all_scores)]
      # Fill NA by last observation carried forward (locf)
      csem_full <- zoo::na.locf(csem_full, na.rm = FALSE)
      result_df <- data.frame(score = all_scores, n = n_all, CSEM = csem_full)
      result_df$CSEM <- round(result_df$CSEM, digits)
      csem_col <- "CSEM"
    } else {
      result_df <- raw
      csem_col <- "CSEM"
    }
  }

  # ---------- 6. Confidence intervals (if requested) ----------
  if (ci) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    if (smooth || full.range) {
      csem_vals <- result_df[[csem_col]]
    } else {
      csem_vals <- result_df$CSEM
    }
    lwr <- result_df$score - z * csem_vals
    upr <- result_df$score + z * csem_vals
    lwr <- pmax(lwr, score_min_teo)
    upr <- pmin(upr, score_max_teo)
    result_df$lwr.ci <- round(lwr, digits)
    result_df$upr.ci <- round(upr, digits)
  }

  # ---------- 7. Quantile groups (bin.score) ----------
  binned_df <- NULL
  if (!is.null(bin.score)) {
    # Need CSEM values for each person (individual csem)
    # or for each observed score? We'll use the raw aggregated CSEM per score
    # for consistency with other functions (average of CSEM per score within group)
    # But the simplest: use the individual csem and compute mean per quantile group.
    # That matches the definition of binned.CSEM in csemFSG: mean of CSEM per person.
    # We'll use the individual csem_person and group by quantiles of total.
    q <- stats::quantile(total, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
    q <- unique(q)
    groups <- cut(total, breaks = q, include.lowest = TRUE, right = TRUE)
    group_levels <- levels(groups)
    bin_list <- list()
    for (i in seq_along(group_levels)) {
      idx <- which(groups == group_levels[i])
      csem_mean <- mean(csem_person[idx], na.rm = TRUE)
      range_str <- paste0(min(total[idx]), "-", max(total[idx]))
      n_pers <- length(idx)
      mean_score <- mean(total[idx])
      bin_list[[i]] <- data.frame(group = i, range = range_str, n = n_pers,
                                  mean_score = mean_score, CSEM.mean = csem_mean,
                                  stringsAsFactors = FALSE)
    }
    binned_df <- do.call(rbind, bin_list)
    binned_df$CSEM.mean <- round(binned_df$CSEM.mean, digits)
    if (ci) {
      z <- stats::qnorm(1 - (1 - conf.level) / 2)
      lwr_bin <- binned_df$mean_score - z * binned_df$CSEM.mean
      upr_bin <- binned_df$mean_score + z * binned_df$CSEM.mean
      lwr_bin <- pmax(lwr_bin, score_min_teo)
      upr_bin <- pmin(upr_bin, score_max_teo)
      binned_df$lwr.ci <- round(lwr_bin, digits)
      binned_df$upr.ci <- round(upr_bin, digits)
    }
  }

  # ---------- 8. Output ----------
  out <- list(CSEM = result_df)
  if (!is.null(binned_df)) out$binned.CSEM <- binned_df
  return(out)
}
