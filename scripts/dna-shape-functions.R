denoise_signal_fft <- function(signal, N = 5, Fs = 500, plot = FALSE) {
  # Normalize signal to [-1, 1]
  norm_signal <- 2 * (signal - min(signal)) / (max(signal) - min(signal)) - 1
  L <- length(norm_signal)
  t <- seq(0, (L - 1)) / Fs
  # print(t)

  # Apply FFT
  fft_result <- fft(norm_signal)
  mags <- Mod(fft_result)

  # Frequencies corresponding to FFT bins
  freqs <- (0:(L - 1)) * Fs / L

  # Only use first half of spectrum for real-valued signals
  half_spectrum <- 2:(L %/% 2)  # Exclude DC (index 1), no Nyquist

  # Find top N frequencies (excluding DC)
  top_indices <- order(mags[half_spectrum], decreasing = TRUE)[1:N]
  top_freq_indices <- half_spectrum[top_indices]
  top_freqs <- freqs[top_freq_indices]

  # Build filtered FFT with only top N frequencies and symmetric components
  filtered_fft <- rep(0+0i, length(fft_result))
  keep_indices <- c(1, top_freq_indices, L - top_freq_indices + 2)
  filtered_fft[keep_indices] <- fft_result[keep_indices]

  # Inverse FFT to get denoised signal
  denoised_signal <- Re(fft(filtered_fft, inverse = TRUE) / L)

  # Rebuild the sine-based signal from frequencies and FFT phase/magnitude
  sine_signal <- rep(0, L)
  for (i in seq_along(top_freq_indices)) {
    idx <- top_freq_indices[i]
    amplitude <- 2 * Mod(fft_result[idx]) / L
    phase <- Arg(fft_result[idx])
    sine_signal <- sine_signal + amplitude * sin(2 * pi * freqs[idx] * t + phase)
  }

  if (plot) {
    # Plotting
    par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
    plot(t, norm_signal,
         type = "l", col = "cyan",
         ylab = "Amplitude", xlab = "Time",
         main = "Original Noisy Signal")
    grid()
    plot(t, denoised_signal,
         type = "l", col = "limegreen",
         ylab = "Amplitude", xlab = "Time",
         main = paste("Denoised Signal (Top", N, "Components)"))
    grid()
  }

  # Return everything
  invisible(list(
    top_frequencies_hz = round(top_freqs, 2),
    sine_function_signal = sine_signal,
    denoised_signal = denoised_signal,
    normalized_signal = norm_signal,
    time_vector = t
  ))
}

seqshape_df <- function(seqshape, norm_length = 100) {
  # Normalization
  norm_seqshape <-
    2 * (seqshape - min(seqshape)) / (max(seqshape) - min(seqshape)) - 1

  # Positions in 'relative' values
  rel_positions <- seq_along(seqshape) / norm_length

  return(data.frame(mgw_values = norm_seqshape,
                    positions = rel_positions))
}

trunc_fourier <- function(data, response, predictor, n_peaks, silent = FALSE,
                          bottomtop = FALSE, aem_option = 1, plot = TRUE,
                          more_sliders = TRUE, even_more_sliders = FALSE,
                          check_again = FALSE, check_again_n = 0,
                          only_coefficients = TRUE) {
  # Note: Need to add 'aem_graph' option

  if (!bottomtop) {
    i_peaks <- seq(n_peaks * 2, 2, by = -2)
    peakstart <- (n_peaks * 2)
  } else {
    i_peaks <- seq(2, n_peaks * 2, by = 2)
    peakstart <- 2
  }

  # Default model
  best_aem <- -22
  max_aem <- best_aem

  all_done <- FALSE
  cntimes <- check_again_n
  if (check_again) cntimes <- 1
  if (cntimes) i_peaks <- c(i_peaks, rep(0, cntimes))
  cstorage_size <- ifelse(cntimes,
                          length(i_peaks) - cntimes, length(i_peaks))

  # Vectors where we'll store coefficients
  dim_arr <- 2 * cstorage_size
  coeffs_storage <- array(rep(0, 2 * dim_arr), dim = c(dim_arr, 2),
                          dimnames = list(1:dim_arr, c("l_coefs", "s_coefs")))

  pcoef_tosave <- 0
  scoef_tosave <- 0

  slide_floats <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
  if (more_sliders)
    slide_floats <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
  if (even_more_sliders)
    slide_floats <- seq(0, 0.5, by = 0.025)
  slide_terms <- as.character(slide_floats[-1])
  slide_terms <- sapply(slide_terms, function(x) paste0(" - ", x, "))"))
  slide_terms <- c("))", slide_terms)
  for (ipeak in i_peaks) {
    if (!bottomtop) {
      pcoefs <- ipeak:(ipeak - 2)
      if (!(ipeak - 2)) pcoefs <- pcoefs[-length(pcoefs)]
      #  ^This removes the case of 'ipeak == 0'
    } else {
      pcoefs <- (ipeak - 2):ipeak
      if (!silent) cat("ipeak: ", ipeak, "peakstart: ", peakstart, "\n")
    }

    # Only used if 'check_again=TRUE'
    # One last run to try non-used "pi_coefs"
    if (!ipeak) {
      if (all_done) print("No further checkings necessary")
      if (all_done) break
      pcoefs <- which(coeffs_storage[, "l_coefs"] == 0)
    }

    i <- 0
    for (pcoef in pcoefs) {
      i <- i + 1
      if (pcoef == pcoef_tosave) next
      if (!ipeak) all_done <- TRUE
      if (ipeak != peakstart) {
        f_terms <- paste(last_f, "+",
                         paste0("cos(pi * ((", pcoef, " * ", predictor, ")))"))
      } else {
        f_terms <- paste0("cos(pi * ((", pcoef, " * ", predictor, ")))")
      }

      t_list <- sapply(slide_terms,
                       function(s_term) gsub("))$", s_term, f_terms))
      f_list <- paste(response, "~", t_list)
      lm_list <- lapply(f_list,
                        function(f_i) lm(as.formula(f_i), data = data))
      aem_list <- sapply(lm_list,
                         function(i_lm) adjust_error_metric(i_lm))

      newmax_aem <- max(aem_list)

      f_index <- which(aem_list == newmax_aem)

      if (!silent) {
        cat("\n==============================|==============================\n")
        print(aem_list)
        print(f_index)
        print(slide_floats[f_index])
      }

      if (newmax_aem > max_aem) {
        if (!silent) cat("New max AEM!\n")
        if (!ipeak) all_done <- FALSE
        # ^This conditional ('all_done') will help us break the loop
        # if no new terms are added while "checking again"

        pcoef_tosave <- pcoef

        scoef_tosave <- slide_floats[f_index]

        newbest_lm <- lm_list[[f_index]]
        lm_coefficients <- newbest_lm$coefficients
        lcoef_tosave <- lm_coefficients[length(lm_coefficients)]

        max_aem <- newmax_aem

        selected_t <- t_list[f_index]
      }

      if (!silent) {
        cat(t_list, sep = "\n")
        cat("Max AEM: ", max_aem, "; pcoef: ", pcoef_tosave,"; scoef: ", scoef_tosave, "\n")
      }
    }
    if (max_aem > best_aem) {
      best_aem <- max_aem
      best_lm <- newbest_lm
      last_f <- selected_t
      coeffs_storage[pcoef_tosave, "l_coefs"] <- lcoef_tosave
      coeffs_storage[pcoef_tosave, "s_coefs"] <- scoef_tosave
      if (!silent) {
        print(coeffs_storage[, "l_coefs"])
        print(coeffs_storage[, "s_coefs"])
      }
    }
  }
  coeffs_inds <- coeffs_storage[, "l_coefs"] != 0
  coeffs_storage[coeffs_inds, "l_coefs"] <- best_lm$coefficients[-1]

  # Visualize final array of pi and slide coefficients
  if (silent) {
    # print(t(coeffs_storage))
  } else if (!only_coefficients) {
    print(coeffs_storage)
  }

  if (plot) {
    fitdata <- data.frame(predictor = seq(min(data[, predictor]),
                                          max(data[, predictor]),
                                          length.out = 250))
    colnames(fitdata)[1] <- predictor
    fitdata$pred_response <- predict(best_lm, newdata = fitdata)

    plot(y = data[, response], x = data[, predictor])
    points(data[, predictor], fitted(best_lm), col = "red", pch = 20)
    lines(x = fitdata[, predictor], y = fitdata$pred_response)
  }
  if (only_coefficients) {
    return(coeffs_storage)
  } else {
    return(best_lm)
  }
}

# Note: This needs to be edited respecting its 'silent' parameter
adjust_error_metric <- function(lmodel, best_adj = FALSE,
                                silent = TRUE, brief = FALSE, aem_option = 1) {
  summ_lmodel <- summary(lmodel)

  adj_rsq <- summ_lmodel$adj.r.squared
  rsq <- summ_lmodel$r.squared
  res_err <- summ_lmodel$sigma
  p_val <- lm_pvalue(summ_lmodel)
  #: t_mean <- tvalues(summ_lmodel, mean = TRUE)
  #: t_prod <- tvalues(summ_lmodel, prod = TRUE)

  if (adj_rsq > 0.999 && rsq > 0.999) {
    if (!silent)
      cat("\t ~~>AEM:  OVER-FITTED")
    return(0.00001)
  }

  errs_names <- c("Adj. R-Squared", "Res. Std. Error",
                  "R-Squared", "P-Value")
  errs_vectr <- c(adj_rsq, res_err, rsq, p_val)
  #: errs_names <- c("Adj. R-Squared", "Res. Std. Error",
  #:                 "R-Squared", "P-Value", "Mean T", "Prod T")
  #: errs_vectr <- c(adj_rsq, res_err, rsq, p_val, t_mean, t_prod)

  aem <- switch(aem_option,
    (adj_rsq * res_err) / (rsq * p_val),
    (adj_rsq * p_val) / (rsq * res_err),
    adj_rsq / (rsq * p_val * res_err),
    res_err^2 / p_val^2,
    adj_rsq,
    rsq
  )
  names(errs_vectr) <- errs_names

  if (brief) {
    cat("\t ~~>AEM: ", aem)
    silent <- TRUE
  }
  if (!silent) {
    print(lmodel$coefficients)
    print("---  ----  ----  ----  ---  -··· ···-  ----  ----  ----  ----  ---")
    print(errs_vectr)
    cat("\n                     ~~>AEM: ", aem, "\n")
    print("=============================|||-|||==============================")
  }
  if (best_adj)
    return(adj_rsq)
  return(aem)
}

lm_pvalue <- function(lm_summ) {
  f <- lm_summ$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  attributes(p) <- NULL
  return(p)
}
