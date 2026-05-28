#' Thorndike's method for Conditional Standard Error of Measurement (CSEM)
#'
#' @param half1, half2 Data frames/matrices con las dos mitades del test.
#' @param bin.score Entero (número de grupos cuantiles) o NULL (puntajes individuales).
#' @param min.n Tamaño mínimo por grupo (solo cuando smooth=FALSE y bin.score=NULL).
#' @param smooth Lógico: si TRUE, suaviza con polinomio.
#' @param degree Grado del polinomio (si smooth=TRUE).
#' @param full.range Lógico: si TRUE, evalúa en todo score.range (requiere smooth=TRUE y score.range).
#' @param ci Lógico: si TRUE, calcula intervalos de confianza.
#' @param conf.level Nivel de confianza (por defecto 0.95).
#' @param digits Decimales para redondeo.
#' @param score.range Vector numérico de longitud 2: mínimo y máximo teórico del puntaje total.
#'
#' @return Lista con elementos CSEM, binned.CSEM (si bin.score entero) y polynomial (si smooth=TRUE).
#' @export
csemThorndike <- function(half1, half2,
                          bin.score = NULL,
                          min.n = 20,
                          smooth = FALSE,
                          degree = 2,
                          full.range = FALSE,
                          ci = FALSE,
                          conf.level = 0.95,
                          digits = 3,
                          score.range = NULL) {

  # --- Validación de argumentos ---
  if (missing(half1) || missing(half2))
    stop("Both 'half1' and 'half2' must be provided.")
  half1 <- as.data.frame(half1)
  half2 <- as.data.frame(half2)
  if (nrow(half1) != nrow(half2))
    stop("'half1' and 'half2' must have the same number of rows.")

  total1 <- rowSums(half1, na.rm = TRUE)
  total2 <- rowSums(half2, na.rm = TRUE)
  total <- total1 + total2
  diff <- total1 - total2

  # Rango teórico
  if (is.null(score.range)) {
    theo_min <- min(total, na.rm = TRUE)
    theo_max <- max(total, na.rm = TRUE)
  } else {
    if (!is.numeric(score.range) || length(score.range) != 2)
      stop("score.range must be a numeric vector of length 2.")
    theo_min <- score.range[1]
    theo_max <- score.range[2]
  }

  # --- Función de fusión de grupos pequeños ---
  merge_small_groups <- function(scores, diffs, min_n, min_score, max_score) {
    all_scores <- seq(min_score, max_score, by = 1)
    freq <- table(scores)
    n_vec <- sapply(all_scores, function(s) if (s %in% names(freq)) freq[as.character(s)] else 0)
    csem_raw <- sapply(all_scores, function(s) {
      idx <- which(scores == s)
      if (length(idx) >= 2) sd(diffs[idx]) else NA_real_
    })
    df <- data.frame(score = all_scores, n = n_vec, CSEM.raw = csem_raw, stringsAsFactors = FALSE)

    # Puntajes que necesitan fusionarse (0 < n < min_n y no tienen CSEM válido)
    to_fuse <- which(df$n > 0 & df$n < min_n & is.na(df$CSEM.raw))
    if (length(to_fuse) == 0) return(df)

    blocks <- list()
    current_block <- c()
    current_n_sum <- 0
    i <- 1
    while (i <= nrow(df)) {
      if (df$n[i] >= min_n) {
        if (length(current_block) > 0) {
          blocks <- c(blocks, list(list(scores = current_block, total_n = current_n_sum)))
          current_block <- c()
          current_n_sum <- 0
        }
        # este puntaje se queda individual (ya está en df)
      } else if (df$n[i] > 0 && df$n[i] < min_n) {
        current_block <- c(current_block, i)
        current_n_sum <- current_n_sum + df$n[i]
        if (current_n_sum >= min_n) {
          blocks <- c(blocks, list(list(scores = current_block, total_n = current_n_sum)))
          current_block <- c()
          current_n_sum <- 0
        }
      } else if (df$n[i] == 0) {
        # los no observados no se incluyen en bloques, pero no interrumpen
      }
      i <- i + 1
    }
    if (length(current_block) > 0) {
      blocks <- c(blocks, list(list(scores = current_block, total_n = current_n_sum)))
    }

    # Asignar CSEM combinado a cada bloque
    for (blk in blocks) {
      idx_scores <- blk$scores
      all_idx <- unlist(lapply(idx_scores, function(j) which(scores == df$score[j])))
      csem_block <- sd(diffs[all_idx])
      for (j in idx_scores) {
        df$CSEM.raw[j] <- csem_block
        df$n[j] <- NA_integer_
      }
    }
    return(df)
  }

  # --- Obtener estimaciones raw ---
  if (smooth) {
    # Para suavizar, necesitamos CSEM en cada puntaje observado con n>=2
    unique_scores <- sort(unique(total))
    raw_list <- list()
    for (s in unique_scores) {
      idx <- which(total == s)
      n_s <- length(idx)
      if (n_s >= 2) {
        csem_raw <- sd(diff[idx])
      } else {
        csem_raw <- NA_real_
      }
      raw_list[[length(raw_list)+1]] <- data.frame(score = s, n = n_s, CSEM.raw = csem_raw,
                                                   stringsAsFactors = FALSE)
    }
    raw_df <- do.call(rbind, raw_list)
    raw_df <- raw_df[!is.na(raw_df$CSEM.raw), ]
  } else {
    if (is.null(bin.score)) {
      raw_df <- merge_small_groups(total, diff, min.n, theo_min, theo_max)
    } else {
      # bin.score entero: primero raw por puntaje, luego agrupar por cuantiles
      unique_scores <- sort(unique(total))
      temp_list <- list()
      for (s in unique_scores) {
        idx <- which(total == s)
        n_s <- length(idx)
        if (n_s >= 2) {
          csem_raw <- sd(diff[idx])
        } else {
          csem_raw <- NA_real_
        }
        temp_list[[length(temp_list)+1]] <- data.frame(score = s, n = n_s, CSEM.raw = csem_raw,
                                                       stringsAsFactors = FALSE)
      }
      temp_df <- do.call(rbind, temp_list)
      temp_df <- temp_df[!is.na(temp_df$CSEM.raw), ]
      # Grupos cuantiles sobre personas
      q <- stats::quantile(total, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
      q <- unique(q)
      groups <- cut(total, breaks = q, include.lowest = TRUE, right = TRUE)
      group_levels <- levels(groups)
      binned_list <- list()
      for (i in seq_along(group_levels)) {
        idx_in_group <- which(groups == group_levels[i])
        scores_in_group <- unique(total[idx_in_group])
        sub_df <- temp_df[temp_df$score %in% scores_in_group, ]
        if (nrow(sub_df) == 0) next
        csem_mean <- mean(sub_df$CSEM.raw)
        range_str <- paste0(min(scores_in_group), "-", max(scores_in_group))
        n_persons <- length(idx_in_group)
        mean_score <- mean(total[idx_in_group])
        binned_list[[i]] <- data.frame(group = i, range = range_str, n = n_persons,
                                       mean_score = mean_score, CSEM.mean = csem_mean,
                                       stringsAsFactors = FALSE)
      }
      binned_df <- do.call(rbind, binned_list)
      raw_df <- temp_df
    }
  }

  # --- Suavizamiento (si aplica) ---
  poly_out <- NULL
  csem_final <- NULL
  if (smooth) {
    fit_df <- raw_df[is.finite(raw_df$CSEM.raw) & !is.na(raw_df$CSEM.raw), ]
    if (nrow(fit_df) < degree + 1)
      stop("Not enough data points for polynomial degree ", degree)
    y <- fit_df$CSEM.raw^2
    x <- fit_df$score
    poly_terms <- paste("I(x^", 1:degree, ")", sep = "", collapse = " + ")
    form <- as.formula(paste("y ~", poly_terms))
    fit <- lm(form, data = data.frame(x = x, y = y))

    if (full.range) {
      full_scores <- seq(theo_min, theo_max, by = 1)
      pred_var <- predict(fit, newdata = data.frame(x = full_scores))
      pred_var <- pmax(pred_var, 0)
      csem_full <- sqrt(pred_var)
      full_n <- sapply(full_scores, function(s) sum(total == s))
      smooth_df <- data.frame(score = full_scores, n = full_n,
                              CSEM.smooth = round(csem_full, digits),
                              stringsAsFactors = FALSE)
    } else {
      pred_var <- predict(fit, newdata = data.frame(x = fit_df$score))
      pred_var <- pmax(pred_var, 0)
      csem_smooth <- sqrt(pred_var)
      smooth_df <- fit_df
      smooth_df$CSEM.smooth <- round(csem_smooth, digits)
      smooth_df$n <- as.integer(smooth_df$n)
    }
    # Información del polinomio
    coef_sum <- summary(fit)$coefficients
    coef_df <- data.frame(term = rownames(coef_sum), estimate = coef_sum[,1],
                          std.error = coef_sum[,2], t.value = coef_sum[,3],
                          p.value = coef_sum[,4], row.names = NULL)
    fit_stats <- list(r.squared = summary(fit)$r.squared,
                      adj.r.squared = summary(fit)$adj.r.squared,
                      AIC = AIC(fit), BIC = BIC(fit),
                      logLik = as.numeric(logLik(fit)),
                      residual.se = summary(fit)$sigma,
                      df.residual = df.residual(fit))
    poly_out <- list(coefficients = coef_df, fit.statistics = fit_stats, degree = degree)
    csem_final <- smooth_df

    # Si hay bin.score, reconstruir binned.CSEM usando valores suavizados
    if (!is.null(bin.score)) {
      q <- stats::quantile(total, probs = seq(0, 1, length.out = bin.score + 1), type = 7)
      q <- unique(q)
      groups <- cut(total, breaks = q, include.lowest = TRUE, right = TRUE)
      group_levels <- levels(groups)
      binned_list2 <- list()
      for (i in seq_along(group_levels)) {
        idx_in_group <- which(groups == group_levels[i])
        scores_in_group <- unique(total[idx_in_group])
        # Predecir CSEM.smooth para esos puntajes usando el modelo
        pred_grp <- predict(fit, newdata = data.frame(x = scores_in_group))
        pred_grp <- pmax(pred_grp, 0)
        csem_pred <- sqrt(pred_grp)
        csem_mean_smooth <- mean(csem_pred, na.rm = TRUE)
        range_str <- paste0(min(scores_in_group), "-", max(scores_in_group))
        n_persons <- length(idx_in_group)
        mean_score <- mean(total[idx_in_group])
        binned_list2[[i]] <- data.frame(group = i, range = range_str, n = n_persons,
                                        mean_score = mean_score, CSEM.mean = csem_mean_smooth,
                                        stringsAsFactors = FALSE)
      }
      binned_df <- do.call(rbind, binned_list2)
    }
  } else {
    csem_final <- raw_df
    names(csem_final)[names(csem_final) == "CSEM.raw"] <- "CSEM"
    csem_final$n <- as.integer(csem_final$n)
  }

  # --- Intervalos de confianza ---
  if (ci) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    if (smooth) {
      csem_vals <- csem_final$CSEM.smooth
    } else {
      csem_vals <- csem_final$CSEM
    }
    lwr <- csem_final$score - z * csem_vals
    upr <- csem_final$score + z * csem_vals
    csem_final$lwr.ci <- round(lwr, digits)
    csem_final$upr.ci <- round(upr, digits)

    if (exists("binned_df") && !is.null(binned_df) && nrow(binned_df) > 0) {
      lwr_bin <- binned_df$mean_score - z * binned_df$CSEM.mean
      upr_bin <- binned_df$mean_score + z * binned_df$CSEM.mean
      binned_df$lwr.ci <- round(lwr_bin, digits)
      binned_df$upr.ci <- round(upr_bin, digits)
    }
  }

  # --- Salida ---
  out <- list(CSEM = csem_final)
  if (exists("binned_df") && !is.null(binned_df) && nrow(binned_df) > 0)
    out$binned.CSEM <- binned_df
  if (smooth && !is.null(poly_out))
    out$polynomial <- poly_out
  class(out) <- "csemThorndike"
  return(out)
}
