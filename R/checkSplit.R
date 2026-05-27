#' Split a Test into Two Halves
#'
#' @description
#' Divides the items of a test into two halves using different methods.
#'
#' @param data A data frame or matrix of item responses (rows = persons, columns = items).
#' @param method Character: "ritc", "difficulty", "random", or "optimized".
#' @param equal.size Logical. If TRUE, forces halves to have the same number of items.
#' @param n_iter Integer. Number of random splits when method = "optimized".
#' @param top_n Integer. Number of top splits to report in summary when method = "optimized".
#' @param seed Integer for reproducibility.
#' @param na.rm Logical. If TRUE, removes rows with missing values.
#'
#' @return A list with components: half1, half2 (vectors of column names), summary (data frame),
#'   and method.
#' @export
checkSplit <- function(data,
                       method = c("ritc", "difficulty", "random", "optimized"),
                       equal.size = TRUE,
                       n_iter = 500,
                       top_n = 10,
                       seed = 123,
                       na.rm = TRUE) {

  method <- match.arg(method)
  set.seed(seed)

  # --- Data validation ---
  if (!is.data.frame(data) && !is.matrix(data)) stop("'data' must be a data frame or matrix.")
  data <- as.data.frame(data)
  if (any(!sapply(data, is.numeric))) stop("All columns must be numeric.")

  n_items <- ncol(data)
  if (n_items < 2) stop("At least 2 items required.")

  if (na.rm) {
    complete <- stats::complete.cases(data)
    if (!all(complete)) {
      data <- data[complete, , drop = FALSE]
      warning("Rows with missing values removed (na.rm = TRUE).")
    }
  }

  # Helper to compute correlation between halves from column names
  get_cor <- function(h1, h2) {
    t1 <- rowSums(data[, h1, drop = FALSE], na.rm = TRUE)
    t2 <- rowSums(data[, h2, drop = FALSE], na.rm = TRUE)
    cor(t1, t2, use = "pairwise.complete.obs")
  }

  # Helper to create summary row (without mean_diff, p_value)
  make_summary_row <- function(h1, h2) {
    data.frame(
      n1 = length(h1),
      n2 = length(h2),
      cor_halves = round(get_cor(h1, h2), 4),
      stringsAsFactors = FALSE
    )
  }

  # --- Splitting methods (except optimized) ---
  if (method != "optimized") {
    if (method == "ritc") {
      total_minus_item <- rowSums(data, na.rm = TRUE)
      ritc <- sapply(1:n_items, function(i) {
        item <- data[, i]
        total <- total_minus_item - item
        cor(item, total, use = "pairwise.complete.obs")
      })
      ritc[is.na(ritc)] <- -Inf
      ord <- order(ritc, decreasing = TRUE)
      half1 <- names(data)[ord[c(TRUE, FALSE)]]
      half2 <- names(data)[ord[c(FALSE, TRUE)]]
    } else if (method == "difficulty") {
      item_means <- colMeans(data, na.rm = TRUE)
      ord <- order(item_means)
      half1 <- names(data)[ord[c(TRUE, FALSE)]]
      half2 <- names(data)[ord[c(FALSE, TRUE)]]
    } else { # random
      idx <- sample(rep(c(1,2), length.out = n_items))
      half1 <- names(data)[idx == 1]
      half2 <- names(data)[idx == 2]
    }

    if (equal.size) {
      len1 <- length(half1); len2 <- length(half2)
      if (len1 != len2) {
        if (len1 > len2) half1 <- half1[1:len2] else half2 <- half2[1:len1]
        message("Note: An item was removed to equalize halves.")
      }
    } else {
      if (length(half1) != length(half2)) message("Note: Halves have unequal number of items.")
    }

    summary_df <- make_summary_row(half1, half2)
    return(list(half1 = half1, half2 = half2, summary = summary_df, method = method))
  }

  # --- Optimized method ---
  best_cor <- -1
  best_split <- NULL
  all_splits <- list()

  for (i in 1:n_iter) {
    idx <- sample(rep(c(1,2), length.out = n_items))
    h1 <- names(data)[idx == 1]
    h2 <- names(data)[idx == 2]
    if (length(h1) == 0 || length(h2) == 0) next
    # Optionally apply equal.size here? The user may want equal size in optimized splits.
    # We'll apply equal.size if requested, otherwise keep as is.
    if (equal.size) {
      len1 <- length(h1); len2 <- length(h2)
      if (len1 != len2) {
        if (len1 > len2) h1 <- h1[1:len2] else h2 <- h2[1:len1]
      }
    }
    ccur <- get_cor(h1, h2)
    if (is.finite(ccur)) {
      all_splits[[i]] <- list(half1 = h1, half2 = h2, cor = ccur)
      if (ccur > best_cor) {
        best_cor <- ccur
        best_split <- list(half1 = h1, half2 = h2)
      }
    }
  }

  if (is.null(best_split)) stop("No valid split found in 'optimized' method.")

  # Build summary data frame with top_n splits (sorted descending by cor)
  # Extract all non-null splits
  valid_splits <- Filter(Negate(is.null), all_splits)
  if (length(valid_splits) == 0) stop("No valid splits.")
  # Sort by cor descending
  sorted <- valid_splits[order(sapply(valid_splits, function(x) x$cor), decreasing = TRUE)]
  top <- head(sorted, top_n)

  summary_df <- data.frame(
    rank = 1:length(top),
    cor_halves = round(sapply(top, function(x) x$cor), 4),
    n1 = sapply(top, function(x) length(x$half1)),
    n2 = sapply(top, function(x) length(x$half2)),
    half1 = sapply(top, function(x) paste(x$half1, collapse = ", ")),
    half2 = sapply(top, function(x) paste(x$half2, collapse = ", ")),
    stringsAsFactors = FALSE
  )

  # The best split (first row) is also returned as vectors
  best_half1 <- top[[1]]$half1
  best_half2 <- top[[1]]$half2

  return(list(
    half1 = best_half1,
    half2 = best_half2,
    summary = summary_df,
    method = method
  ))
}
