# =============================================================================
# DEPENDENCIES -------------------------------------------------------------- #
# =============================================================================

# library(Biostrings)
# library(data.table)
# library(dplyr)
# library(mosaic)
# library(box)
library(datawizard)
library(fitdistrplus)
library(stringr)
library(stringi)
library(primes)

# =============================================================================
# BASIC UTILITIES: sequence base characterizations/transformations ---------- #
# =============================================================================

bases_count <- function(sequence, bases = c("A", "T", "C", "G")) {
  base_counts <- str_count(sequence, bases)
  names(base_counts) <- bases
  return(base_counts)
}

bases_percentage <- function(sequence,
                             bases = c("A", "T", "C", "G"), base_counts) {
  if (missing(sequence)) {
    base_percs <- base_counts / sum(base_counts)
  } else {
    base_percs <- bases_count(sequence, bases) / str_length(sequence)
  }
  names(base_percs) <- bases
  return(base_percs)
}

gc_percentage <- function(sequence, bases = c("A", "T", "C", "G")) {
  base_percs <- bases_percentage(sequence, bases)
  gc_perc <- base_percs[3] + base_percs[4]
  # names(gc_perc) <- "GC%"
  names(gc_perc) <- NULL
  return(gc_perc)
}

highlight_base <- function(sequence, base) {
  BASE <- toupper(base)
  base <- tolower(base)
  sequence <- tolower(sequence)
  return(gsub(base, BASE, sequence))
}

rev_complement <- function(sequence, bases = "ATCG", replace_bases = "TAGC") {
  return(stri_reverse(chartr(bases, replace_bases, sequence)))
}

# =============================================================================
# KMERS / WINDOWS: kmer creation, identification & utilities for them ------- #
# =============================================================================

combi_kmers <- function(bases = c("A", "T", "C", "G"), k = 2) {
  # 'sort' orders kmers list so that they follow the same order as the ones
  # transformed into factors in the (kinda homologous) function kmer_wincount()
  return(sort(do.call(paste0, expand.grid(rep(list(bases), k)))))
}

count_kmers <- function(sequence, vector_kmers, k = 2, percentage = FALSE) {
  if (missing(vector_kmers)) {
    vector_kmers <- combi_kmers(k = k)
  }
  kmer_counts <- stri_count_fixed(sequence, vector_kmers, overlap = TRUE)
  names(kmer_counts) <- vector_kmers
  if (!percentage)
    return(kmer_counts)
  if (percentage) {
    kmers_number <- stri_length(sequence) - k + 1
    return(kmer_counts / kmers_number)
  }
}

kmer_diversity <- function(kmer_counts, percentage = TRUE) {
  if (percentage)
    return(length(kmer_counts[kmer_counts > 0]) / length(kmer_counts))
  if (!percentage)
    return(length(kmer_counts[kmer_counts > 0]))
}

counts_per_window <- function(sequence, k = 2, s = 1, win_kmers, all_kmers) {
  if (missing(win_kmers)) {
    win_kmers <- kmer_windows(sequence, k = k, s = s)
  }
  if (missing(all_kmers)) {
    all_kmers <- combi_kmers(k = k)
  }
  win_kmers <- factor(win_kmers)
  all_kmers <- factor(all_kmers)
  counts_per_window <- table(all_kmers[match(win_kmers, all_kmers)])
  return(counts_per_window)
}

kmer_windows <- function(sequence, k = 2, s = 1) {
  # 'k' stands for 'kmer' (in this case its equivalent to the window size)
  # 's' stands for 'stride' (spaces taken between each window)
  seq_len <- str_length(sequence)
  return(str_sub(sequence, seq(1, seq_len + 1 - k, s), seq(k, seq_len, s)))
}

func_per_windows <- function(sequence, k = 2, s = 1, windows, kmers, func) {
  if (missing(windows)) {
    windows <- kmer_windows(sequence, k = k, s = s)
  }
  if (missing(kmers)) {
    return(sapply(windows, func, simplify = TRUE))
  } else {
    return(sapply(kmers,
                  function(x) func(kmer = x, windows = windows),
                  simplify = TRUE))
  }
}

# =============================================================================
# TM CALCULATION: basic formula for more than 13 bps or less than 14 bps ---- #
# =============================================================================

tm_calc <- function(sequence, base_counts, title = FALSE) {
  if (missing(base_counts)) {
    base_counts <- bases_count(sequence)
  }
  seq_len <- str_length(sequence)
  countA <- base_counts[1]
  countT <- base_counts[2]
  countC <- base_counts[3]
  countG <- base_counts[4]
  if (seq_len <= 13) {
    temp <- ((countA + countT) * 2) + ((countC + countG) * 4)
  } else {
    temp <- 64.9 + (41 * (countG + countC - 16.4) / (countA + countT + countC + countG))
  }
  if (title) {
    names(temp) <- "TM in Celcius"
  } else {
    names(temp) <- NULL
  }
  return(temp)
}

# tm_len_lt14 <- function(sequence, bases, base_counts, title = FALSE) {
tm_len_lt14 <- function(sequence, base_counts, title = FALSE) {
  if (missing(base_counts)) {
    base_counts <- bases_count(sequence)
  }
  countA <- base_counts[1]
  countT <- base_counts[2]
  countC <- base_counts[3]
  countG <- base_counts[4]
  temp <- ((countA + countT) * 2) + ((countC + countG) * 4)
  if (title) {
    names(temp) <- "TM in Celcius"
  } else {
    names(temp) <- NULL
  }
  return(temp)
}

# tm_len_mt13 <- function(sequence, bases, base_counts, title = FALSE) {
tm_len_mt13 <- function(sequence, base_counts, title = FALSE) {
  if (missing(base_counts)) {
    # base_counts <- bases_count(sequence, bases)
    base_counts <- bases_count(sequence)
  }
  countA <- base_counts[1]
  countT <- base_counts[2]
  countC <- base_counts[3]
  countG <- base_counts[4]
  temp <- 64.9 + (41 * (countG + countC - 16.4) / (countA + countT + countC + countG))
  if (title) {
    names(temp) <- "TM in Celcius"
  } else {
    names(temp) <- NULL
  }
  return(temp)
}

# =============================================================================
# PALINDROME UTILITIES: low-level functions for palindrome identification --- #
# =============================================================================

isRevComPalindrome <- function(sequence) {
  bases <- "ATCG"
  replace_bases <- "TAGC"
  return(sequence == chartr(bases, replace_bases, stri_reverse(sequence)))
}

isPalindrome <- function(sequence) {
  return(sequence == stri_reverse(sequence))
}

palindromes <- function(sequence, k = 2, s = 1, windows, count = FALSE) {
  if (missing(windows)) {
    windows <- kmer_windows(sequence, k = k, s = s)
  }
  pals <- windows[isPalindrome(windows)]
  if (!count) {
    return(pals)
  } else {
    return(length(pals))
  }
}

# =============================================================================
# SHANNON ENTROPY ----------------------------------------------------------- #
# =============================================================================

shannon_entropy <- function(sequence, base_counts, base_percs) {
  if (missing(base_percs)) {
    if (missing(base_counts)) {
      bpercs <- bases_percentage(sequence)
    } else {
      bpercs <- bases_percentage(base_counts = base_counts)
    }
  } else {
    bpercs <- base_percs
  }
  return(-sum(bpercs * ifelse(bpercs > 0, log(bpercs, 2), 1)))
}

# =============================================================================
# OPTIMIZED PRODUCT: kmer(count percentage in seq * GC percentage * shannon)  #
# =============================================================================

optim_prod <- function(sequence, bases_counts, base_percs)

# =============================================================================
# KMER BARCODE: transform kmer positions into numeric vector to get sum/prod  #
# =============================================================================

kmer_barcode <- function(sequence, k, windows, kmer, option = "expsum") {
  if (missing(windows)) {
    windows <- kmer_windows(sequence, k = k)
  }
  barcode_profile <- switch(option,
    "expsum" = expsum_profile(windows, kmer),
    "primes" = primes_profile(windows, kmer),
    "logprimes" = logprimes_profile(windows, kmer),
    {
      stop("No valid option selected")
    }
  )
  return(barcode_profile)
}
expsum_profile <- function(windows, kmer) {
  # : return(sum(1.001^(which(rev(windows == kmer)) - 1)))
  return(sum(1.001^(which(windows == kmer) - 1)))
}
primes_profile <- function(windows, kmer) {
  n_kmers <- length(windows)
  return(prod(generate_n_primes(n_kmers)[windows == kmer]))
}
logprimes_profile <- function(windows, kmer) {
  n_kmers <- length(windows)
  return(prod(log(generate_n_primes(n_kmers), 8)[windows == kmer]))
}

intvectodec <- function(intvec) {
  return(sum(intvec * 10^((length(intvec) - 1):0)))
}

bintodec <- function(x) {
  return(sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1)) - 1)))
}

# =============================================================================
# PARTIAL PALINDROMES: seq - rev(seq) pairs given an interspace length limit  #
# =============================================================================

# Note: maybe add option 'mark_true_palindromes' in the form of another
# boolean row in the final dataframe. To optimize this operation, we
# would check first if window is its own palindrome or differs only
# one nucleotide from its reverse match. Then each sequence would
# be applied to the function: isPalindrome()
partial_palindromes <- function(windows, restrict_length, sequence, min_size,
                                only_indexes = FALSE, only_sequences = FALSE) {
  if (missing(sequence)) {
    if (missing(windows)) {
      stop("No 'sequence' or 'windows' parameters provided")
    } else if (!missing(min_size)) {
      warning("Parameter 'min_size' unnecessary given that 'windows'
			   was provided. Size of kmers in 'windows' will be used instead")
    }
    only_indexes <- TRUE
  } else if (missing(windows)) {
    if (missing(min_size)) {
      stop("Just 'sequence' provided, 'min_size' parameter expected")
    }
    windows <- kmer_windows(sequence, k = min_size)
  }
  min_size <- str_length(windows[1])
  vector_len <- length(windows)
  if (missing(restrict_length)) {
    restrict_length <- vector_len
  }
  rev_kmers <- stri_reverse(windows)
  if (!missing(only_indexes)) {
    return(pp_chunk_onlyInds(vector_len, min_size,
                             restrict_length, rev_kmers, windows))
  } else if (!missing(only_sequences)) {
    return(pp_chunk_onlySeqs(vector_len, min_size,
                             restrict_length, rev_kmers, windows, sequence))
  }  else {
    return(pp_chunk_Both(vector_len, min_size,
                         restrict_length, rev_kmers, windows, sequence))
  }
}
pp_chunk_Both <- function(vector_len, kmer_size,
                          restrict_len, rev_kmers, win_kmers, sequence) {
  s_inds_row <- c()
  e_inds_row <- c()
  seqs_row <- c()
  for (start_ind in (1:vector_len)) {
    end_ind <- start_ind + kmer_size + restrict_len
    if (end_ind > vector_len) {
      end_ind <- vector_len
    }
    end_inds <- (start_ind :end_ind)[rev_kmers[start_ind]
                                     == win_kmers[start_ind:end_ind]]
    ei_len <- length(end_inds)
    if (ei_len) {
      start_inds <- rep(start_ind, ei_len)
      end_inds <- end_inds - 1 + kmer_size
      pal_seqs <- str_sub(sequence, start_ind, end_inds)
      s_inds_row <- c(s_inds_row, start_inds)
      e_inds_row <- c(e_inds_row, end_inds)
      seqs_row <- c(seqs_row, pal_seqs)
    }
  }
  pals_df <- data.frame(start = s_inds_row, end = e_inds_row,
                        sequences = seqs_row)
  return(pals_df)
}
pp_chunk_onlyInds <- function(vector_len, kmer_size,
                              restrict_len, rev_kmers, win_kmers) {
  s_inds_row <- c()
  e_inds_row <- c()
  for (start_ind in (1:vector_len)) {
    end_ind <- start_ind + kmer_size + restrict_len
    if (end_ind > vector_len) {
      end_ind <- vector_len
    }
    end_inds <- (start_ind :end_ind)[rev_kmers[start_ind]
                                     == win_kmers[start_ind:end_ind]]
    ei_len <- length(end_inds)
    if (ei_len) {
      start_inds <- rep(start_ind, ei_len)
      end_inds <- end_inds - 1 + kmer_size
      s_inds_row <- c(s_inds_row, start_inds)
      e_inds_row <- c(e_inds_row, end_inds)
    }
  }
  pals_df <- data.frame(start = s_inds_row, end = e_inds_row)
  return(pals_df)
}
pp_chunk_onlySeqs <- function(vector_len, kmer_size,
                              restrict_len, rev_kmers, win_kmers, sequence) {
  seqs_vec <- c()
  for (start_ind in (1:vector_len)) {
    end_ind <- start_ind + kmer_size + restrict_len
    if (end_ind > vector_len) {
      end_ind <- vector_len
    }
    end_inds <- (start_ind :end_ind)[rev_kmers[start_ind]
                                     == win_kmers[start_ind:end_ind]]
    if (length(end_inds)) {
      end_inds <- end_inds - 1 + kmer_size
      pal_seqs <- str_sub(sequence, start_ind, end_inds)
      seqs_vec <- c(seqs_vec, pal_seqs)
    }
  }
  return(seqs_vec)
}

# =============================================================================
# TABLE MAKERS: Table production (data extraction), given a set of sequences  #
# =============================================================================

sequence_characterizer <- function(sequence, k_max = 4, no_names = FALSE,
                                   optim = FALSE, as_vector = TRUE) {
  basepercs <- bases_percentage(sequence)
  seq_tms <- tm_len_mt13(sequence)
  seq_sha <- shannon_entropy(sequence)
  ## Add pairwise alignment data per sequence
  ## Local and Global alignment of sequence against its mirror (reverse)
  ## and reverse complementary

  kmers_characs <- list()
  if (!as_vector) {
    kmers_characs[[1]] <- list(perc = basepercs, temp = seq_tms, shan = seq_sha)
    names(kmers_characs)[1] <- "WS"
  } else {
    kiv_characs <- c(basepercs, temp = seq_tms, shan = seq_sha)
  }
  for (i in 2:k_max) {
    all_ki <- combi_kmers(k = i)
    ki_seq <- kmer_windows(sequence, k = i)
    ki_percs <- count_kmers(sequence, k = i,
                            percentage = TRUE)
    # ^This can be optimized to use only loaded kmers not sequence
    if (!optim)
      ki_tms <- func_per_windows(windows = all_ki,
                                 func = tm_len_lt14)
    ki_sha <- func_per_windows(windows = all_ki,
                               func = shannon_entropy)
    ki_brc <- func_per_windows(kmers = all_ki,
                               windows = ki_seq,
                               func = kmer_barcode)
    ki_pals <- func_per_windows(windows = all_ki,
                                func = isPalindrome)
        # ^Add 'isPalindrome' to 'optim' conditional in the form of
        # 'if(isPalindrome)' then "value = 2*value"
    ki_revc <- func_per_windows(windows = all_ki,
                                func = isRevComPalindrome)
        # ^Add 'isRevComPalindrome' to 'optim' conditional in the form of
        # 'if(isRevComPalindrome)' then "value = 2*value",
        #  but also add option
        # 'if((seq_length % 2) != 0)' then
        # "replace nucleotide in the middle with N" then
        # 'if(isRevComPalindrome)' then "value = 2*value"
        # Maybe later on, make more options like this for partial palindromes.
    # ^Another optimization is to only compute data from non-zero count kmers
    # and later on, organize it in the corresponding 'all_ki' order
    # Maybe I can just add an "isTrue" parameter/conditional to all functions
    ki_name <- paste("k", i, sep = "")
    if (optim) {
      ki_gcp <- func_per_windows(windows = all_ki,
                                 func = gc_percentage)
      ki_gcp <- 2^ki_gcp
      ki_sha <- 1.1^ki_sha
      ki_prod <- ki_percs * ki_gcp * ki_sha
      ki_characs <- list(prod = ki_prod, barc = ki_brc,
                         pals = ki_pals, revc = ki_revc)
    } else {
      ki_characs <- list(perc = ki_percs, temp = ki_tms, shan = ki_sha,
                         barc = ki_brc, pals = ki_pals, revc = ki_revc)
    }
    # ^This 'optim' option can be further exploited if instead of doing
    # 4 different 'func_per_windows()' and then operating over their values,
    # make a new function that multiplies all perc*gcp*sha per sequence in
    # just one go.
    if (as_vector) {
      for (kmer_i in 1:4^i) {
        # ikmer_seq <- c(all_ki[kmer_i])
        # names(ikmer_seq) <- paste("k", i, "-", kmer_i, "_", "seq", sep = "")
        ikmer_characs <- c(unlist(lapply(ki_characs, `[[`, kmer_i)))
        if (!no_names)
          names(ikmer_characs) <- paste("k", i, "-", kmer_i,
                                        "_", names(ikmer_characs), sep = "")
        # kiv_characs <- c(kiv_characs, ikmer_seq, ikmer_characs)
        kiv_characs <- c(kiv_characs, ikmer_characs)
      }
      # kmers_characs[[i]] <- kiv_characs
    } else {
      kmers_characs[[i]] <- ki_characs
      names(kmers_characs)[i] <- ki_name
    }
  }
  # return(ifelse(as_vector, kiv_characs, kmers_characs))
  if (as_vector) {
    return(kiv_characs)
  } else {
    return(kmers_characs)
  }
}

sequences_characterizer <- function(sequences, k_max, optim = FALSE, as_df = TRUE) {
  # seqs_data <- list(A_perc = c(), T_perc = c(), C_perc = c(), G_perc = c(),
  #                   seq_TM = c(), seq_SE = c())
  sc_len <- ifelse(optim,
                   cols_numb(k = k_max, n_params = 4),
                   cols_numb(k = k_max, n_params = 6))
  sc_len <- sc_len + 6
  seqs_data <- vector(mode = "list", length = sc_len)
  i <- 0
  for (sequence_i in sequences) {
    if (!i) {
      seq_i <- sequence_characterizer(sequence_i, k = k_max,
                                      optim = optim, as_vector = TRUE)
      names(seqs_data) <- names(seq_i)
      names(seq_i) <- NULL
      i <- i + 1
    } else {
      seq_i <- sequence_characterizer(sequence_i, k = k_max, no_names = TRUE,
                                      optim = optim, as_vector = TRUE)
      names(seq_i) <- NULL
    }
    for (j in 1:sc_len) {
      seqs_data[[j]] <- append(seqs_data[[j]], seq_i[j])
    }
  }
  if (as_df) {
    return(as.data.frame(seqs_data))
  } else {
    return(seqs_data)
  }
}
cols_numb <- function(k_max, n_params) {
  number <- 0
  for (i in 2:k_max) {
    number <- sum(number, 4^i)
  }
  return(number * n_params)
}

# =============================================================================
# REGRESSIONS: Polynomial, Fourier-like and Bin-like modelling for kmer data  #
# =============================================================================

# To-do: modelling in a bin-like (histogram) way
kmer_spc_polynome <- function(sequence, k, windows, kmer,
                              plot = FALSE, option = "norm") {
  if(missing(windows)) {
    if (missing(sequence)) {
      stop("No 'sequence' or 'windows' parameters provided")
    }
    windows <- kmer_windows(sequence, k = k)
    len_seq <- stri_length(sequence)
    len_kmer <- k
  } else {
    len_kmer <- stri_length(windows[1])
    len_seq <- length(windows) + len_kmer - 1
  }

  inds <- index_pairs(which(windows == kmer), len_seq, len_kmer)
  if (inds$length <= 1) {
    # The conditional filters no-match cases
    coef_names <- c("C", "x", "x^2", "x^3", "x^4", "x^5")
    errs_names <- c("Residual Std Error", "Adj R-Squared", "P-Value")
    if (inds$length) {
      if (plot)
        message("No plot can be done with just one match")
      sole <- switch(option,
        "norm" = ((inds$ends - inds$starts) / len_seq)^2,
        "log-exp" = 1.01^(inds$ends - inds$starts),
        "plain" = (inds$ends - inds$starts),
        {
          stop("No valid option selected")
        }
      )
      sole <- c(sole, 0, 0, 0, 0, 0, 0, 1, 0)
      names(sole) <- c(coef_names, errs_names)
      return(sole)
    } else {
      if (plot)
        message("No plot can be done with no match found")
      no_coefs <- c(0, 0, 0, 0, 0, 0, -1, 0, -1)
      names(no_coefs) <- c(coef_names, errs_names)
      return(no_coefs)
    }
  } else {
    return(switch(option,
      "norm" = norm_spc_profile(inds, len_seq, len_kmer, plot),
      "log-exp" = logexp_spc_profile(inds, len_seq, len_kmer, plot),
      "plain" = plain_spc_profile(inds, len_seq, len_kmer, plot),
      {
        stop("No valid option selected")
      }
    ))
  }
}
index_pairs <- function(indexes, len_seq, len_kmer) {
  index_stts <- c(0, indexes + len_kmer - 1)
  index_ends <- c(indexes - 1, len_seq)
  index_remv <- which(index_ends <= index_stts)
  if (length(index_remv)) {
    index_stts <- index_stts[-index_remv]
    index_ends <- index_ends[-index_remv]
  }
  return(list(starts = index_stts, ends = index_ends,
              length = length(index_stts)))
}
norm_spc_profile <- function(inds, len_seq, len_kmer, plot) {
  case <- inds$length - 1
  # case needs to be "length - 1" to avoid overfitting
  spc_pos <- inds$starts / len_seq
  spc_len <- ((inds$ends - inds$starts) / len_seq)^2

  lmodel <- model_case(spc_len, spc_pos, case)
  if (plot)
    model_predict_plot(lmodel$model, spc_pos, spc_len)
  return(c(lmodel$coefficients, lmodel$errors))
}
logexp_spc_profile <- function(inds, len_seq, len_kmer, plot) {
  case <- inds$length - 1

  spc_pos <- log(generate_n_primes(len_seq)[inds$starts + 1], 16)
  spc_len <- 1.01^(inds$ends - inds$starts)

  lmodel <- model_case(spc_len, spc_pos, case)
  if (plot)
    model_predict_plot(lmodel$model, spc_pos, spc_len)
  return(c(lmodel$coefficients, lmodel$errors))
}
plain_spc_profile <- function(inds, len_seq, len_kmer, plot) {
  case <- inds$length - 1

  spc_pos <- inds$starts
  spc_len <- (inds$ends - inds$starts)

  lmodel <- model_case(spc_len, spc_pos, case)
  if (plot)
    model_predict_plot(lmodel$model, spc_pos, spc_len)
  return(c(lmodel$coefficients, lmodel$errors))
}
model_case <- function(lens, poss, case) {
  # case < 0 is alreaddy filtered by the time the function is called
  coef_names <- c("C", "x", "x^2", "x^3", "x^4", "x^5")
  errs_names <- c("Residual Std Error", "Adj R-Squared", "P-Value")
  if (case == 1) {
    lmodel <- model_testing(lens, poss, case)
    coefs <- c(lmodel$coefficients, 0, 0, 0, 0)
    summ <- summary(lmodel)
    errs <- c(summ$sigma, summ$adj.r.squared, lm_pvalue(summ))
    names(coefs) <- coef_names
    names(errs) <- errs_names
  } else if (case == 2) {
    lmodel <- model_testing(lens, poss, case)
    coefs <- coef_fill(lmodel$coefficients)
    summ <- summary(lmodel)
    errs <- c(summ$sigma, summ$adj.r.squared, lm_pvalue(summ))
    names(coefs) <- coef_names
    names(errs) <- errs_names
  } else if (case == 3) {
    lmodel <- model_testing(lens, poss, case)
    coefs <- coef_fill(lmodel$coefficients)
    summ <- summary(lmodel)
    errs <- c(summ$sigma, summ$adj.r.squared, lm_pvalue(summ))
    names(coefs) <- coef_names
    names(errs) <- errs_names
  } else if (case == 4) {
    lmodel <- model_testing(lens, poss, case)
    coefs <- coef_fill(lmodel$coefficients)
    summ <- summary(lmodel)
    errs <- c(summ$sigma, summ$adj.r.squared, lm_pvalue(summ))
    names(coefs) <- coef_names
    names(errs) <- errs_names
  } else if (case >= 5) {
    lmodel <- model_testing(lens, poss, 5)
    coefs <- coef_fill(lmodel$coefficients)
    summ <- summary(lmodel)
    errs <- c(summ$sigma, summ$adj.r.squared, lm_pvalue(summ))
    names(coefs) <- coef_names
    names(errs) <- errs_names
  }
  return(list(model = lmodel, coefficients = coefs, errors = errs))
}
model_testing <- function(lens, poss, case) {
  lm_prev <- polynome_cases(lens, poss, 1)
  for (i in seq(2, sum(seq(1, case)) - 1)) {
    lm_post <- polynome_cases(lens, poss, i)
    prev_rsq <- summary(lm_prev)$adj.r.squared
    post_rsq <- summary(lm_post)$adj.r.squared
    if (post_rsq > prev_rsq)
      lm_prev <- lm_post
  }
  return(lm_prev)
}
coef_fill <- function(coefs) {
  coefs[is.na(coefs)] <- 0
  to_fill <- (6 - length(coefs))
  if (to_fill) {
    coefs <- c(coefs, rep(0, to_fill))
  } else {
    return(coefs)
  }
}
lm_pvalue <- function(lm_summ) {
  f <- lm_summ$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  attributes(p) <- NULL
  return(p)
}
model_predict_plot <- function(lmodel, poss, lens) {
  predicted <- data.frame(poss = seq(min(poss),
                                     max(poss),
                                     length.out = 100))
  print(lmodel)
  predicted$lens <- predict(lmodel, newdata = predicted)

  plot(y = lens, x = poss)
  points(poss, fitted(lmodel), col = "red", pch = 20)
  lines(x = predicted$poss, y = predicted$lens)
}

# model_peaks <- function(lens, poss, peaks) {}
polynome_cases <- function(lens, poss, case) {
  return(switch(case,
    lm(lens ~ poss),
    lm(lens ~ I(0 * poss) + I(poss^2)),
    lm(lens ~ poss + I(poss^2)),
    lm(lens ~ I(0 * poss) + I(0 * poss^2) + I(poss^3)),
    lm(lens ~ poss + I(0 * poss^2) + I(poss^3)),
    lm(lens ~ poss + I(poss^2) + I(poss^3)),
    lm(lens ~ poss + I(0 * poss^2) + I(0 * poss^3) + I(poss^4)),
    lm(lens ~ poss + I(poss^2) + I(0 * poss^3) + I(poss^4)),
    lm(lens ~ poss + I(0 * poss^2) + I(poss^3) + I(poss^4)),
    lm(lens ~ poss + I(poss^2) + I(poss^3) + I(poss^4)),
    lm(lens ~ poss + I(poss^2) + I(0 * poss^3) + I(0 * poss^4) + I(poss^5)),
    lm(lens ~ poss + I(0 * poss^2) + I(poss^3) + I(0 * poss^4) + I(poss^5)),
    lm(lens ~ poss + I(poss^2) + I(poss^3) + I(0 * poss^4) + I(poss^5)),
    lm(lens ~ poss + I(poss^2) + I(0 * poss^3) + I(poss^4) + I(poss^5)),
    lm(lens ~ poss + I(poss^2) + I(poss^3) + I(poss^4) + I(poss^5))
  ))
}

trigon_cases <- function(lens, poss, peaks, case, subcase,
                       shifter = 0, modulator = 1) {
  # : shf <- shifter
  # : mdl <- modulator
  return(switch(peaks,
                switch(case, # CASES when 1 PEAK
                  # MIN CASE POSSIBLE == 1 == 3 points
                  # MIN cos() POSSIBLE == cos(1 * pi * x)
                  switch(subcase,
                    # SUBCASES when case >= 1
                    lm(lens
                       ~ cos(1 * pi * poss)),
                    lm(lens
                       ~ sin(1 * pi * poss)),
                    lm(lens
                       ~ cos(2 * pi * poss)),
                    # Dont forget to try 'shifter' like
                    # : cos(1 * pi * poss + shf) -> cos(1 * pi * poss + 0.75)
                    # : sin(1 * pi * poss - shf) -> sin(1 * pi * poss - 0.75)
                    # Dont forget to try 'modulator' like
                    # : cos((1 + mdl) * pi * poss) -> cos((1 + 0.5) * pi * poss)
                    # : sin((1 - mdl) * pi * poss) -> sin((1 - 0.5) * pi * poss)
                  ),
                  switch(subcase,
                    # SUBCASES when case >= 2
                    lm(lens
                       ~ cos(1 * pi * poss) + sin(1 * pi * poss)),
                    lm(lens
                       ~ cos(1 * pi * poss) + cos(2 * pi * poss)),
                  ),
                ),
                switch(case, # CASES when 2 PEAKS
                  # MIN CASE POSSIBLE == 1 == 3 points
                  # MIN cos() POSSIBLE == cos(2 * pi * x)
                  switch(subcase,
                    # SUBCASES when case >= 1
                    lm(lens
                       ~ sin(1 * pi * poss)),
                    lm(lens
                       ~ cos(2 * pi * poss)),
                  ),
                  switch(subcase,
                    # SUBCASES when case >= 2
                    lm(lens
                       ~ cos(1 * pi * poss) + sin(1 * pi * poss)), # 3 points
                    lm(lens
                       ~ cos(1 * pi * poss) + cos(2 * pi * poss)), # 3 points
                    lm(lens
                       ~ sin(2 * pi * poss)),
                    lm(lens
                       ~ cos(3 * pi * poss)),
                    lm(lens
                       ~ sin(1 * pi * poss) + sin(2 * pi * poss)),
                    lm(lens
                       ~ cos(1 * pi * poss) + sin(2 * pi * poss)),
                    lm(lens
                       ~ cos(1 * pi * poss) + cos(3 * pi * poss)),
                    lm(lens
                       ~ cos(2 * pi * poss) + sin(2 * pi * poss)),
                    lm(lens
                       ~ cos(2 * pi * poss) + cos(3 * pi * poss)),
                    # Omitted, not very informative
                    # : cos(1 * pi * poss) + sin(2 * pi * poss)
                    # Dont forget to try shifter
                    # : cos(1 * pi * poss - 0.75) + sin(1 * pi * poss)
                    # : cos(1 * pi * poss) + sin(1 * pi * poss + 0.75)
                  ),
                  switch(subcase,
                    # SUBCASES when case >= 3
                    lm(lens
                       ~ sin(3 * pi * poss)),
                    lm(lens
                       ~ cos(4 * pi * poss)),
                    lm(lens
                       ~ sin(1 * pi * poss) + cos(2 * pi * poss)),
                    lm(lens
                       ~ cos(1 * pi * poss) + cos(4 * pi * poss)),
                  ),
                ),
                switch(case, # CASES when 3 PEAKS
                  # MIN CASE POSSIBLE == 3 == 5 points
                  # MIN cos() POSSIBLE == cos(4 * pi * x)
                  switch(subcase,
                    lm(lens
                       ~ cos(3 * pi * poss)),
                    lm(lens
                       ~ cos(4 * pi * poss))
                  ),
                  switch(subcase,
                    lm(lens
                       ~ cos(1 * pi * poss) + cos(4 * pi * poss)),
                    lm(lens
                       ~ cos(2 * pi * poss) + sin(1 * pi * poss)),
                    lm(lens
                       ~ cos(3 * pi * poss) + sin(2 * pi * poss)),
                  ),
                ),
                switch(case,
                  lm(lens
                     ~ cos(2 * pi * poss) + sin(2 * pi * poss)),
                  lm(lens
                     ~ cos(4 * pi * poss) + sin(4 * pi * poss)),
                  lm(lens
                     # ~ cos(5 * pi * poss) + sin(5 * pi * poss)
                     ~ cos(6 * pi * poss) + sin(6 * pi * poss)),
                  lm(lens
                     # ~ cos(1 * pi * poss) + sin(5 * pi * poss)
                     ~ cos(2 * pi * poss) + sin(4 * pi * poss)),
                  lm(lens
                     # ~ cos(2 * pi * poss) + sin(5 * pi * poss)
                     ~ cos(2 * pi * poss) + sin(6 * pi * poss)),
                  lm(lens
                     ~ cos(4 * pi * poss) + sin(6 * pi * poss)),
                  lm(lens
                     ~ cos(1 * pi * poss) + cos(3 * pi * poss)),
                  lm(lens
                     ~ cos(1 * pi * poss) + sin(4 * pi * poss)),
                  lm(lens
                     ~ cos(2 * pi * poss) + sin(3 * pi * poss)),
                  lm(lens
                     ~ cos(3 * pi * poss) + sin(3 * pi * poss)),
                  lm(lens
                     ~ cos(3 * pi * poss) + sin(4 * pi * poss)),
                ),
                switch(case,
                  lm(lens
                    ~ cos(2 * pi * poss) + sin(2 * pi * poss)
                    + cos(4 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(2 * pi * poss) + sin(2 * pi * poss)
                    + cos(6 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(2 * pi * poss) + sin(4 * pi * poss)
                    + cos(4 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(2 * pi * poss) + sin(4 * pi * poss)
                    + cos(6 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(2 * pi * poss) + sin(4 * pi * poss)
                    + cos(8 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(4 * pi * poss) + sin(4 * pi * poss)
                    + cos(6 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(4 * pi * poss) + sin(4 * pi * poss)
                    + cos(8 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(4 * pi * poss) + sin(6 * pi * poss)
                    + cos(6 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(4 * pi * poss) + sin(6 * pi * poss)
                    + cos(8 * pi * poss)
                  )
                ),
                switch(case,
                  lm(lens
                    ~ cos(2 * pi * poss) + sin(2 * pi * poss)
                    + cos(4 * pi * poss) + sin(4 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(2 * pi * poss) + sin(2 * pi * poss)
                    + cos(6 * pi * poss) + sin(6 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(4 * pi * poss) + sin(4 * pi * poss)
                    + cos(6 * pi * poss) + sin(6 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(4 * pi * poss) + sin(4 * pi * poss)
                    + cos(8 * pi * poss) + sin(8 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(2 * pi * poss) + sin(4 * pi * poss)
                    + cos(4 * pi * poss) + sin(6 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(2 * pi * poss) + sin(6 * pi * poss)
                    + cos(6 * pi * poss) + sin(8 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(4 * pi * poss) + sin(6 * pi * poss)
                    + cos(6 * pi * poss) + sin(8 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(2 * pi * poss) + sin(4 * pi * poss)
                    + cos(6 * pi * poss) + sin(8 * pi * poss)
                  ),
                  lm(lens
                    ~ cos(6 * pi * poss) + sin(6 * pi * poss)
                    + cos(8 * pi * poss) + sin(8 * pi * poss)
                  )
                ),
                lm(lens
                  ~ cos(2 * pi * poss) + sin(2 * pi * poss)
                  + cos(4 * pi * poss) + sin(4 * pi * poss)
                  + cos(6 * pi * poss)
                ),
                lm(lens
                  ~ cos(2 * pi * poss) + sin(2 * pi * poss)
                  + cos(4 * pi * poss) + sin(4 * pi * poss)
                  + cos(6 * pi * poss) + sin(6 * pi * poss)
                ),
                lm(lens
                  ~ cos(2 * pi * poss) + sin(2 * pi * poss)
                  + cos(4 * pi * poss) + sin(4 * pi * poss)
                  + cos(6 * pi * poss) + sin(6 * pi * poss)
                  + cos(8 * pi * poss)
                ),
                lm(lens
                  ~ cos(2 * pi * poss) + sin(2 * pi * poss)
                  + cos(4 * pi * poss) + sin(4 * pi * poss)
                  + cos(6 * pi * poss) + sin(6 * pi * poss)
                  + cos(8 * pi * poss) + sin(8 * pi * poss)
                ),
                lm(lens
                  ~ cos(2 * pi * poss) + sin(2 * pi * poss)
                  + cos(4 * pi * poss) + sin(4 * pi * poss)
                  + cos(6 * pi * poss) + sin(6 * pi * poss)
                  + cos(8 * pi * poss) + sin(8 * pi * poss)
                  + cos(10 * pi * poss)
                ),
                lm(lens
                  ~ cos(2 * pi * poss) + sin(2 * pi * poss)
                  + cos(4 * pi * poss) + sin(4 * pi * poss)
                  + cos(6 * pi * poss) + sin(6 * pi * poss)
                  + cos(8 * pi * poss) + sin(8 * pi * poss)
                  + cos(10 * pi * poss) + sin(10 * pi * poss)
                )))
}
#  bases_percentage(sequence)
# return(-sum(sapply(bpercs, function(x) x * log(x,2), simplify = TRUE)))

# tms_k8_sequence <- as.numeric(lapply(kmer_windows(sequence, k = 8), tm_len_lt14))
# as.data.frame(describe_distribution(tms_k8_sequence))











