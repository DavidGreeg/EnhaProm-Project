library(stringi)
library(stringr)

scrpath <- "/home/davidfm/Projects/UBMI-IFC/EnhaProm/scripts/"
genpath <- "genome-functions.R"
genofunctions <- paste0(scrpath, genpath)
source(genofunctions)

kmer_windows <- kmer_windows
# ============================================================================ #
#						  INTERSPACE POLYNOMES
# ============================================================================ #
get_interspaces <- function(sequence, k, windows, kmer,
                            print = FALSE, rm_Nth = 0, relative = FALSE) {
  if (missing(kmer))
    stop("no 'kmer' provided to get interspaces from")
  if (missing(k) || missing(windows)) {
    if (missing(k)) {
      k <- stri_length(kmer)
      if (missing(windows))
        stop("missing arguments; needed either 'k' or (plus
             'sequence'), or 'windows' argument to continue")
      if (missing(sequence))
        lenseq <- length(windows) + k - 1
    }
    if (missing(windows)) {
      if (missing(sequence))
        stop("missing arguments; if 'k' is
             provided 'sequence' is needed")
      windows <- kmer_windows(sequence, k = k)
      lenseq <- stri_length(sequence)
    }
  }
  indexes <- which(windows == kmer)
  start_inds <- c(0, indexes + k - 1)
  final_inds <- c(indexes - 1, lenseq)

  if (print) {
    print(start_inds)
    print(final_inds)
  }

  illogic_inds <- which(final_inds <= start_inds)
  if (length(illogic_inds)) {
    start_inds <- start_inds[-illogic_inds]
    final_inds <- final_inds[-illogic_inds]
  }

  if (rm_Nth) {
    start_inds <- start_inds[-rm_Nth]
    final_inds <- final_inds[-rm_Nth]
  }

  #: spc_pos <- log(generate_n_primes(seq_len)[i_stts + 1], 16)
  #: spc_len <- 1.01^(i_ends - i_stts)
  intspace_lengths <- final_inds - start_inds

  #: spc_positions <- start_inds / lenseq
  # NOTE: Adjusted to the middle between start and final index
  intspace_positions <- start_inds + (intspace_lengths / 2)

  if (relative) {
    # Positions in terms of sequence length percentage
    intspace_positions <- intspace_positions / lenseq
    # Normalize lengths by sequence length
    intspace_lengths <- intspace_lengths / lenseq
  }

  if (print) {
    cat("InterSpaces: \n")
    cat("lengths:", intspace_lengths, "\n")
    cat("positions:", intspace_positions, "\n")
  }
  return(array(data = c(intspace_lengths, intspace_positions),
               dimnames = list(seq_along(start_inds),
                               c("lengths", "positions")),
               dim = c(length(start_inds), 2)))
}

# Better name might be 'kmer_interspace_fourier'
kmer_interspace_polynome <- function(sequence, k, windows, kmer,
                                     plot = TRUE, fixed = FALSE,
                                     case, subcase, minus = 0,
                                     kdesmooth = FALSE, bw, rel = TRUE,
                                     silent = TRUE, brief = FALSE,
                                     best_adj = FALSE, print_inds = TRUE,
                                     rm_nth = 0, acc_peaks = FALSE,
                                     aem_option = 5, aem_graph = FALSE) {
  # NOTE: AEM means 'adjust-error metric', which is the
  #       option chosen to select the best model calculated
  if (missing(kmer))
    stop("no 'kmer' provided to get interspaces from")
  if (missing(sequence) || missing(windows)) {
    if (missing(sequence)) {
      if (missing(windows))
        stop("missing arguments; needed either
             'sequence' or 'windows' argument to continue")
      intspace_arr <- get_interspaces(windows = windows, kmer = kmer,
                                      relative = rel, print = print_inds)
    }
    if (missing(windows)) {
      k <- stri_length(kmer)
      if (missing(sequence))
        stop("missing arguments; if 'windows'
             is not provided, 'sequence' is needed")
      intspace_arr <- get_interspaces(sequence, k, kmer = kmer,
                                      relative = rel, print = print_inds)
    }
  }

  spc_pos <- intspace_arr[, "positions"]
  spc_len <- intspace_arr[, "lengths"]

  # Option to accentuate peaks
  if (acc_peaks)
    spc_len <- spc_len^2

  if (aem_graph)
    plot <- FALSE

  if (missing(case))
    case <- length(spc_pos) - 1
  if (missing(subcase))
    subcase <- 1

  #: print(spc_pos)
  #: print(spc_len)
  #: model <- lm(spc_len ~ spc_pos + I(spc_pos^2) + I(spc_pos^3) + I(spc_pos^4))
  #: model <- lm(spc_len ~ cos(spc_pos))
  #: model <- lm(spc_len ~ cos(spc_pos / pi))
  #: model <- lm(spc_len ~ cos(spc_pos / pi) + sin(spc_pos / pi))
  #: print("==============================MODEL1==============================")
  #: print(summary(model))
  #: model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos))
  #: print("==============================MODEL2==============================")
  #: print(summary(model))
  #: # model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos^2) + cos(spc_pos / pi))
  #: model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos) + cos(spc_pos / pi))
  #: print("==============================MODEL3==============================")
  #: print(summary(model))
  #: # model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos^2) + cos(spc_pos / pi) + sin(spc_pos / pi))
  #: model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos) + cos(spc_pos / pi) + sin(spc_pos / pi))
  #: print("==============================MODEL4==============================")
  #: print(summary(model))
  #: print(spc_pos, pi, cos(spc_pos))
  #: model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos)
  #:             + cos(spc_pos / pi) + sin(spc_pos / pi) + cos(pi / spc_pos))
  #: print("==============================MODEL5==============================")
  #: print(summary(model))
  #: model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos)
  #:             + cos(spc_pos / pi) + sin(spc_pos / pi)
  #:             + cos(2 * pi / spc_pos) + sin(2 * pi / spc_pos))
  #: print("==============================MODEL6==============================")
  #: print(summary(model))
  #: model <- lm(spc_len
              #: ~ cos(2 * pi * spc_pos)
              #: + sin(2 * pi * spc_pos)
              #: ~ cos(spc_pos)
              #: ~ sin(spc_pos)
              #: + sin(spc_pos)
              #: + cos(spc_pos)
              #: + cos(4 * pi * spc_pos)
              #: + sin(4 * pi * spc_pos)
              #: + cos(spc_pos / (1.1^spc_pos))
              #: + sin(spc_pos / (1.1^spc_pos))
              #: + cos(6 * pi * spc_pos)
              #: + sin(6 * pi * spc_pos)
              #: + cos(spc_pos^2)
              #: + sin(spc_pos^2)
              #: + cos(8 * pi * spc_pos)
              #: + sin(8 * pi * spc_pos)
              #: + cos(spc_pos^3)
              #: + sin(spc_pos^3)
              #: + cos(spc_pos^4)
              #: + sin(spc_pos^4)
              #: + I(spc_pos)
              #: + I(spc_pos^2)
              #: + I(spc_pos^3)
              #: )
  if (kdesmooth) {
    # Kernel Density Estimation:
    #             (not really useful for interpolation, but informative)
    require("scales")
    #: z <- ksmooth(spc_pos, spc_len, kernel = "normal", bandwidth = bw)
    ## NOTE: 'bw' is obligatory and means "bandwidth" (i.e. bw = 0.1)
    print(bw)
    z <- density(spc_pos, kernel = "gaussian", bw = bw)
    print(z)
    #: plot(spc_len ~ spc_pos, pch=16, cex=0.8, col=alpha("grey", 0.9))
    #: lines(z, lwd=2, col=alpha("black", 0.9))
    plot(z, lwd = 2, col = alpha("black", 0.9))
    #: detach("package:scales", unload = TRUE)
    return("finished")
  } else {
    if (fixed) {
      model <- trig_cases(spc_len, spc_pos, case - minus, subcase)
    } else {
      model <- trig_testing(spc_len, spc_pos, case,
                            best_adj, silent, brief,
                            aem_option, aem_graph)
    }
  }
  #: print("==============================MODEL6==============================")
  #: print(summary(model))
  #: model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos)
  #:             + cos(spc_pos / pi) + sin(spc_pos / pi)
  #:             + cos(pi / spc_pos) + sin(pi / spc_pos) + cos(2 * pi / spc_pos))
  #: print("==============================MODEL7==============================")
  #: print(summary(model))
  #: model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos)
  #:             + cos(spc_pos / pi) + sin(spc_pos / pi)
  #:             + cos(pi / spc_pos) + sin(pi / spc_pos)
  #:             + cos(2 * pi / spc_pos) + sin(2 * pi / spc_pos))
  #: print("==============================MODEL8==============================")
  #: print(summary(model))
  #: model <- lm(spc_len ~ cos(spc_pos / pi) + sen(spc_pos / pi) + cos(spc_pos / pi))
  #: print(model)
  if (plot) {
    fitdata <- data.frame(spc_pos = seq(min(spc_pos),
                                        max(spc_pos),
                                        length.out = 100))
    fitdata$pred_len <- predict(model, newdata = fitdata)

    plot(y = spc_len, x = spc_pos)
    points(spc_pos, fitted(model), col = "red", pch = 20)
    lines(x = fitdata$spc_pos, y = fitdata$pred_len)
  }
  return(model)
}
trig_cases <- function(spc_len, spc_pos, case, subcase) {
  return(switch(case,
                switch(subcase,
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(3 * pi * spc_pos)
                  )
                  #: lm(spc_len
                  #:   ~ cos(4 * pi * spc_pos)
                  #: )
                ),
                switch(subcase,
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  ),
                  lm(spc_len
                    #: ~ cos(5 * pi * spc_pos) + sin(5 * pi * spc_pos)
                    ~ cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
                  ),
                  lm(spc_len
                    #: ~ cos(1 * pi * spc_pos) + sin(5 * pi * spc_pos)
                    ~ cos(2 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  ),
                  lm(spc_len
                    #: ~ cos(2 * pi * spc_pos) + sin(5 * pi * spc_pos)
                    ~ cos(2 * pi * spc_pos) + sin(6 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(4 * pi * spc_pos) + sin(6 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(1 * pi * spc_pos) + cos(3 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(1 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(3 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + cos(3 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(3 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  ),
                ),
                switch(subcase,
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                    + cos(4 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                    + cos(6 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(4 * pi * spc_pos)
                    + cos(4 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(4 * pi * spc_pos)
                    + cos(6 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(4 * pi * spc_pos)
                    + cos(8 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                    + cos(6 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                    + cos(8 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(4 * pi * spc_pos) + sin(6 * pi * spc_pos)
                    + cos(6 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(4 * pi * spc_pos) + sin(6 * pi * spc_pos)
                    + cos(8 * pi * spc_pos)
                  )
                ),
                switch(subcase,
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                    + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                    + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                    + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                    + cos(8 * pi * spc_pos) + sin(8 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(4 * pi * spc_pos)
                    + cos(4 * pi * spc_pos) + sin(6 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(6 * pi * spc_pos)
                    + cos(6 * pi * spc_pos) + sin(8 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(4 * pi * spc_pos) + sin(6 * pi * spc_pos)
                    + cos(6 * pi * spc_pos) + sin(8 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(4 * pi * spc_pos)
                    + cos(6 * pi * spc_pos) + sin(8 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
                    + cos(8 * pi * spc_pos) + sin(8 * pi * spc_pos)
                  )
                ),
                lm(spc_len
                  ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                  + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  + cos(6 * pi * spc_pos)
                ),
                lm(spc_len
                  ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                  + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
                ),
                lm(spc_len
                  ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                  + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
                  + cos(8 * pi * spc_pos)
                ),
                lm(spc_len
                  ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                  + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
                  + cos(8 * pi * spc_pos) + sin(8 * pi * spc_pos)
                ),
                lm(spc_len
                  ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                  + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
                  + cos(8 * pi * spc_pos) + sin(8 * pi * spc_pos)
                  + cos(10 * pi * spc_pos)
                ),
                lm(spc_len
                  ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                  + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
                  + cos(8 * pi * spc_pos) + sin(8 * pi * spc_pos)
                  + cos(10 * pi * spc_pos) + sin(10 * pi * spc_pos)
                )))
}

trig_testing <- function(spc_len, spc_pos, case,
                         best_adj, silent, brief,
                         aem_option, aem_graph) {
  #: lm_prev <- trig_cases(spc_len, spc_pos, 1, 1)
  lm_prev <- lm(spc_len ~ cos(2 * pi * spc_pos))
  prev_aem <- adj_error_metric(lm_prev, best_adj, silent, brief, aem_option)
  #: for (i in seq(2, sum(seq(1, case)) - 1)) {
  if (aem_graph) {
    aem_vec <- c()
    aem_cases <- c()
  }
  best_case <- c(1, 1)
  for (i in seq(1, case - 1)) {
    #: for (j in seq(1, i * 3)) {
    if (i > 3) {
      j_seq <- seq(1, 9)
    } else if (i > 1) {
      j_seq <- seq(1, 6)
    } else {
      j_seq <- seq(1, 2)
    }
    for (j in j_seq) {
      if (!silent)
        cat("\n===> Case:", i, "-", j)
      lm_post <- trig_cases(spc_len, spc_pos, i, j)
      post_aem <- adj_error_metric(lm_post, best_adj, silent, brief, aem_option)
      if (aem_graph) {
        aem_vec <- append(aem_vec, post_aem)
        case_name <- paste(i, j, sep = "-")
        aem_cases <- append(aem_cases, case_name)
      }
      #: both_rsq <- lm_post$r.squared + lm_post$adj.r.squared
      #: if (both_rsq == )
      if (post_aem > prev_aem) {
        if (!silent)
          cat("\n\t>>>> Better AEM inside Case:", i, "-", j)
        best_case <- c(i, j)
        lm_prev <- lm_post
        prev_aem <- adj_error_metric(lm_prev, best_adj,
                                     silent, brief, aem_option)
      }
    }
  }
  if (aem_graph) {
    aem_vec <- aem_vec[-1]
    aem_cases <- aem_cases[-1]
    aem_ids <- seq(1, length(aem_vec))
    aem_vec * 100
    plot(aem_ids, aem_vec)
    axis(side = 1, at = aem_ids, labels = aem_cases)
  }
  cat("\n Best Case chosen: ", best_case[1], "-", best_case[2], "\n")
  return(lm_prev)
}
adj_error_metric <- function(lmodel, best_adj = FALSE,
                             silent = FALSE, brief = FALSE, aem_option = 5) {
  summ_lmodel <- summary(lmodel)

  adj_rsq <- summ_lmodel$adj.r.squared
  rsq <- summ_lmodel$r.squared
  res_err <- summ_lmodel$sigma
  p_val <- lm_pvalue(summ_lmodel)
  t_mean <- tvalues(summ_lmodel, mean = TRUE)
  t_prod <- tvalues(summ_lmodel, prod = TRUE)

  if (adj_rsq > 0.999 && rsq > 0.999) {
    if (!silent)
      cat("\t ~~>AEM:  OVER-FITTED")
    return(0.00001)
  }

  errs_names <- c("Adj. R-Squared", "Res. Std. Error",
                  "R-Squared", "P-Value", "Mean T", "Prod T")
  errs_vectr <- c(adj_rsq, res_err, rsq, p_val, t_mean, t_prod)

  aem <- switch(aem_option,
    (adj_rsq * res_err) / p_val,
    (adj_rsq * rsq * res_err) / p_val, #******2
    (adj_rsq + res_err) * (rsq / p_val),
    (res_err / adj_rsq) * (rsq / p_val),
    #: ((adj_rsq + rsq) ** adj_rsq) / (rsq * res_err * p_val),
    (adj_rsq * res_err) / (rsq * p_val), #********5
    (adj_rsq * p_val) / (rsq * res_err), #********6
    adj_rsq / (rsq * p_val * res_err), #*******7
    res_err^2 / p_val^2, #*****8
    adj_rsq + (((rsq * res_err)^2 / p_val^2) * 100), #*****
    (adj_rsq * res_err^2) / (rsq * p_val),
    (adj_rsq * res_err^2 * (rsq - adj_rsq)^2) / (p_val * rsq), #*******11
    (adj_rsq * res_err^2) / (p_val * rsq * ((rsq * t_mean) - adj_rsq^2)), #*******11
    (adj_rsq * res_err * (rsq - adj_rsq)^2) / (p_val * rsq), #*******11
    (adj_rsq * res_err^2) / ((rsq - adj_rsq)^2 * p_val * rsq), #*******11
  )
  names(errs_vectr) <- errs_names

  if (brief) {
    cat("\t ~~>AEM: ", aem)
    silent <- TRUE
  }
  if (!silent) {
    print(lmodel$coefficients)
    print("----  ----  ----  ----  ---  -··· ···-  ----  ----  ----  ----  ----")
    print(errs_vectr)
    cat("\n                     ~~>AEM: ", aem, "\n")
    print("==============================|||-|||===============================")
  }
  if (best_adj)
    return(adj_rsq)
  return(aem)
}
lm_pvalue <- function(summ_lmodel) {
  f <- summ_lmodel$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  attributes(p) <- NULL
  return(p)
}
tvalues <- function(summ_lmodel, mean = FALSE, prod = FALSE) {
  tvals <- summ_lmodel$coefficients[, "t value"]
  if (mean) {
    return(mean(tvals))
  }
  if (prod) {
    return(prod(tvals))
  }
}

#: trigonom_cases <- function(lens, poss, case) {
#:   return(switch(case,
#:     lm(lens
#:       ~ cos(2 * pi * spc_pos)
#:     ),
#:     lm(lens
#:       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#:     ),
#:     lm(lens
#:       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#:       + cos(4 * pi * spc_pos)
#:     ),
#:     lm(lens
#:       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#:       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#:     ),
#:     lm(lens
#:       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#:       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#:       + cos(6 * pi * spc_pos)
#:     ),
#:     lm(lens
#:       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#:       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#:       + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
#:     ),
#:     lm(lens
#:       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#:       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#:       + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
#:       + cos(8 * pi * spc_pos)
#:     ),
#:     lm(lens
#:       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#:       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#:       + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
#:       + cos(8 * pi * spc_pos) + sin(8 * pi * spc_pos)
#:     ),
#:     lm(lens
#:       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#:       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#:       + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
#:       + cos(8 * pi * spc_pos) + sin(8 * pi * spc_pos)
#:       + cos(10 * pi * spc_pos)
#:     ),
#:     lm(lens
#:       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#:       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#:       + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
#:       + cos(8 * pi * spc_pos) + sin(8 * pi * spc_pos)
#:       + cos(10 * pi * spc_pos) + sin(10 * pi * spc_pos)
#:     ),
#:     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#:     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#:     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#:     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#:     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#:     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#:   ))
#: }
generate_combinations <- function(set, k) {
  # Generate all combinations of size k
  combinations <- combn(set, k, simplify = FALSE)
  # Print each combination
  for (comb in combinations) {
    print(comb)
  }
  # Return the list of combinations
  return(combinations)
}
rm(generate_combinations)
