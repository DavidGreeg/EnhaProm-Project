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

n_ki <- function(size) {
  return(4^size)
}

# =============================================================================
# KMERS / WINDOWS: kmer creation, identification & utilities for them ------- #
# =============================================================================

present_kmers <- function(sequence, k, all_kmers, seq_kmers,
                          kmer_counts = c(), percent = FALSE) {
  if (length(kmer_counts)) {
    return(names(kmer_counts[kmer_counts > 0]))
  } else {
    if (missing(seq_kmers)) {
      if (missing(sequence) || missing(k))
        stop("No 'seq_kmers' nor ('sequence' or 'k') arguments provided")
      seq_kmers <- kmer_windows(sequence, k = k)
    }
    if (missing(all_kmers)) {
      if (missing(k))
        stop("No 'all_kmers' nor 'k' arguments provided")
      all_kmers <- combi_kmers(k = k)
    }

    kmers_true <- all_kmers %in% seq_kmers
    if (percent) {
      return(sum(as.numeric(kmers_true)) / length(all_kmers))
    } else {
      return(all_kmers[kmers_true])
    }
  }
}

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

kmer_diversity <- function(sequence, seq_kmers, all_kmers, k,
                           kmer_counts, percentage = TRUE, relative = FALSE) {
  if (!missing(kmer_counts)) {
    rel_k <- ifelse(!relative, length(kmer_counts), sum(kmer_counts))
    n_present <- length(kmer_counts[kmer_counts > 0])
    diverse_val <- ifelse(!percentage, n_present,
                          n_present / rel_k)
  } else {
    if (missing(seq_kmers)) {
      if (missing(k) || missing(sequence))
        stop("Missing 'seq_kmers' and either 'k' or 'sequence'.
             Both needed to compute 'seq_kmers'")
      seq_kmers <- kmer_windows(sequence, k = k)
    } else if (missing(k)) {
      k <- stri_length(seq_kmers[1])
    }
    if (!relative) {
      if (missing(all_kmers))
        all_kmers <- combi_kmers(k = k)
      k_divisor <- length(all_kmers)
    } else {
      k_divisor <- length(seq_kmers)
    }
    n_present <- length(unique(seq_kmers))
    diverse_val <- ifelse(!percentage, n_present,
                          n_present / k_divisor)
  }
  return(diverse_val)
}

kmer_abs_diversity <- function(seq_kmers, all_kmers) {
  k_divisor <- length(all_kmers)
  n_present <- length(unique(seq_kmers))
  diverse_val <- n_present / k_divisor
  return(diverse_val)
}

kmer_rel_diversity <- function(seq_kmers, all_kmers) {
  k_divisor <- length(seq_kmers)
  n_present <- length(unique(seq_kmers))
  diverse_val <- n_present / k_divisor
  return(diverse_val)
}

counts_per_window <- function(sequence, k = 2, s = 1, seq_kmers, all_kmers,
                              percentage = FALSE) {
  if (missing(seq_kmers)) {
    if (missing(sequence))
      stop("No 'sequence' nor 'seq_kmers' parameters provided")
    seq_kmers <- kmer_windows(sequence, k = k, s = s)
  }
  if (missing(all_kmers)) {
    all_kmers <- combi_kmers(k = k)
  }
  seq_kmers <- factor(seq_kmers)
  all_kmers <- factor(all_kmers)
  counts_per_window <- table(all_kmers[match(seq_kmers, all_kmers)])
  if (percentage)
    counts_per_window <- counts_per_window / length(seq_kmers)
  return(counts_per_window)
}

# ↓↓↓↓ Optimized function
kmer_counts <- function(seq_kmers, all_kmers, sequence, k, percentage = FALSE) {
  if (missing(seq_kmers)) {
    if (missing(k) || missing(sequence))
      stop("Missing 'seq_kmers' and either 'k' or 'sequence'.
           Both needed to compute 'seq_kmers'")
    seq_kmers <- kmer_windows(sequence, k = k)
  } else if (missing(k)) {
    k <- stri_length(seq_kmers[1])
  }
  if (missing(all_kmers)) {
    all_kmers <- combi_kmers(k = k)
  }
  counts <- tabulate(match(seq_kmers, all_kmers), nbins = length(all_kmers))
  if (percentage) counts <- counts / length(seq_kmers)
  names(counts) <- all_kmers
  return(counts)
}

kmer_windows <- function(sequence, k = 2, s = 1) {
  # 'k' stands for 'kmer' (in this case its equivalent to the window size)
  # 's' stands for 'stride' (spaces taken between each window)
  seq_len <- str_length(sequence)
  return(str_sub(sequence, seq(1, seq_len + 1 - k, s), seq(k, seq_len, s)))
}

func_per_windows <- function(sequence, k, s, windows, kmers,
                             option = character(0), func) {
  if (missing(windows)) {
    windows <- kmer_windows(sequence, k = k, s = s)
  }
  if (missing(kmers)) {
    return(sapply(windows, func, simplify = TRUE))
  } else {
    if (missing(option)) {
      return(sapply(kmers,
                    function(x) func(kmer = x, windows = windows),
                    simplify = TRUE))
    } else {
      return(sapply(kmers,
            function(x) func(kmer = x, windows = windows, option = option),
            simplify = TRUE))
    }
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

#:  tm_len_lt14 <- function(sequence, bases, base_counts, title = FALSE) {
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

#: tm_len_mt13 <- function(sequence, bases, base_counts, title = FALSE) {
tm_len_mt13 <- function(sequence, base_counts, title = FALSE) {
  if (missing(base_counts)) {
    #: base_counts <- bases_count(sequence, bases)
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

kmers_shannon <- function(seq_kmers, all_kmers,
                          sequence, k, kmer_percs) {
  # Should add 'seq_kmers' and 'all_kmers' input validation later
  if (missing(kmer_percs)) {
    if (missing(seq_kmers)) {
      if(!missing(sequence) && !missing(k)) {
        seq_kmers <- kmer_windows(sequence, k)
      } else {
        stop("'seq_kmers' argument not provided, to compute it \n",
             "both 'sequence' and 'k' arguments are needed \n")
      }
    }
    if (missing(all_kmers)) {
      if (missing(k))
        k <- str_length(seq_kmers[1])
      all_kmers <- combi_kmers(k = k)
    }
    kmer_percs <- kmer_counts(seq_kmers = seq_kmers,
                              all_kmers = all_kmers,
                              percentage = TRUE)
    if (!missing(sequence)){
    }
  }
  return(-sum(kmer_percs * ifelse(kmer_percs > 0, log(kmer_percs, 2), 1)))
}

# =============================================================================
# EXP(GC PERCENTAGE) / SEQUENCE SHANNON ENTROPY: GC% normalized by Shannon    #
# =============================================================================

egc_ssh_ratio <- function(sequence, shannon) {
  egc_perc <- exp(gc_percentage(sequence))
  return(egc_perc / shannon)
}

# =============================================================================
# GC PERCENTAGE * SHANNON ENTROPY: Optimized function to use later ---------- #
# =============================================================================

gc_sha_prod <- function(sequence, mod_sh = 0.1, mod_gc = 0.2, version = 0) {
  bpercs <- bases_percentage(sequence)
  names(bpercs) <- NULL

  if (version) {
    shannon <- (0.2 - sum(bpercs * ifelse(bpercs > 0, log(bpercs, 2), 1))) / 2.2
    gc_perc <- (0.1 + bpercs[3] + bpercs[4]) / 1.1
  } else {
    shannon <- mod_sh^(-sum(bpercs * ifelse(bpercs > 0, log(bpercs, 2), 1)))
    gc_perc <- mod_gc^(bpercs[3] + bpercs[4])
  }

  return(shannon * gc_perc)
}

# Need to check which option is faster

# =============================================================================
# OPTIMIZED PRODUCT: kmer(count percentage in seq * GC percentage * shannon)  #
# =============================================================================

ksg_product <- function(kmer, windows, option) {
  k_perc <- kmer_counts(seq_kmers = windows,
                        all_kmers = kmer,
                        percentage = TRUE)
  if (!k_perc)
    return(0)

  k_perc <- k_perc * 100
  names(k_perc) <- NULL
  gc_sha <- gc_sha_prod(kmer)

  return(k_perc * gc_sha)
}

# =============================================================================
# KMER BARCODE: transform kmer positions into numeric vector to get sum/prod  #
# =============================================================================

kmer_barcode <- function(sequence, k, windows, kmer, option = "expsum") {
  if (missing(windows))
    windows <- kmer_windows(sequence, k = k)
  if (missing(sequence)) {
    lenseq <- length(windows) + stri_length(kmer) - 1
  } else {
    lenseq <- stri_length(sequence)
  }
  if (! kmer %in% windows)
    return(0)
  if (option != "expsum")
    n_kmers <- length(windows)
  positions <- which(windows == kmer)
  barcode_profile <- switch(option,
    "expsum" = expsum_profile(positions),
    "primes" = primes_profile(n_kmers, positions),
    "divexpsum" = divexpsum_profile(positions, lenseq),
    "logprimes" = logprimes_profile(positions, lenseq),
    # "logprimes" = logprimes_profile(n_kmers, positions, lenseq),
    {
      stop("No valid option selected")
    }
  )
  return(barcode_profile)
}
expsum_profile <- function(positions) {
  # : return(sum(1.001^(which(rev(windows == kmer)) - 1)))
  return(sum(1.001^positions - 1))
}
divexpsum_profile <- function(positions, lenseq) {
  # Minus 1 since we want to be able to represent
  # the position '1' with 'f(0)'
  return(sum(1.01^positions - 1) / 1.01^lenseq)
}
primes_profile <- function(n_kmers, positions) {
  return(prod(generate_n_primes(n_kmers)[positions]))
}
logprimes_profile <- function(positions, lenseq) {
  posit_primes <- generate_n_primes(lenseq)
  last_prime <- rev(posit_primes)[1]
  posit_primes <- posit_primes[positions]

  return(log(prod(posit_primes), last_prime))
  #: return(prod(log(posit_primes, last_prime)))
  #: return(prod(log(generate_n_primes(n_kmers), lenseq)[positions]))
}

intvectodec <- function(intvec) {
  return(sum(intvec * 10^((length(intvec) - 1):0)))
}

bintodec <- function(x) {
  return(sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1)) - 1)))
}

comps_barcode <- function(sequence, k, windows, kmer, option = "partpal") {
  if (missing(windows)) {
    windows <- kmer_windows(sequence, k = k)
  }
  if (! kmer %in% windows)
    return(0)
  primes <- generate_n_primes(length(windows) + 4)
  return(get_profile(primes[-(1:4)], option, windows, kmer))
}
get_profile <- function(primes, option, windows, kmer) {
  if (! option %in% c("partpal", "revcomp")) {
    stop("No valid bar_option selected")
  }
  pp_kmer <- stri_reverse(kmer)
  rc_kmer <- rev_complement(kmer)

  if (option == "partpal" &&
        (!isPalindrome(kmer) ||
           pp_kmer %in% windows)) {
    return(prod(2 * log(primes[windows == kmer], 7)) +
             prod(3 * log(primes[windows == pp_kmer], 7)))
  } else if (option == "revcomp" &&
               (stri_length(kmer) %% 2 ||
                  rc_kmer %in% windows)) {
    rc_pal <- ifelse(isRevComPalindrome(kmer), 5, 1)
    return(rc_pal * (prod(2 * log(primes[windows == kmer], 7)) +
                       prod(3 * log(primes[windows == rc_kmer], 7))))
  } else if (isPalindrome(kmer)) {
    return(prod(5 * log(primes[windows == kmer], 7)))
  } else {
    return(prod(log(primes[windows == kmer], 7)))
  }
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

sequence_characterizer_version0 <- function(sequence, k_max = 4, no_names = FALSE,
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
    if (!optim) {
      ki_percs <- count_kmers(sequence, k = i,
                              percentage = TRUE)
      # ^This can be optimized to use only loaded kmers not sequence
      ki_tms <- func_per_windows(windows = all_ki,
                                 func = tm_len_lt14)
      ki_sha <- func_per_windows(windows = all_ki,
                                 func = shannon_entropy)
    }
    ki_brc <- func_per_windows(kmers = all_ki,
                               windows = ki_seq,
                               func = kmer_barcode, opt = "expsum")
    ki_pals <- func_per_windows(kmers = all_ki,
                                windows = ki_seq,
                                func = comps_barcode, opt = "partpal")
    ki_revc <- func_per_windows(kmers = all_ki,
                                windows = ki_seq,
                                func = comps_barcode, opt = "revcomp")
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
      ki_prod <- func_per_windows(kmers = all_ki,
                                  windows = ki_seq,
                                  func = ksg_product)
      ki_characs <- list(prod = ki_prod, barc = ki_brc,
                         pals = ki_pals, revc = ki_revc)
    } else {
      ki_characs <- list(perc = ki_percs, temp = ki_tms, shan = ki_sha,
                         barc = ki_brc, pals = ki_pals, revc = ki_revc)
    }
    if (as_vector) {
      for (kmer_i in 1:4^i) {
        ikmer_characs <- c(unlist(lapply(ki_characs, `[[`, kmer_i)))
        if (!no_names)
          names(ikmer_characs) <- paste("k", i, "-", kmer_i,
                                        "_", names(ikmer_characs), sep = "")
        kiv_characs <- c(kiv_characs, ikmer_characs)
      }
    } else {
      kmers_characs[[i]] <- ki_characs
      names(kmers_characs)[i] <- ki_name
    }
  }
  if (as_vector) {
    return(kiv_characs)
  } else {
    return(kmers_characs)
  }
}


### IMPORTANT NOTE: Should apply rm(ls()) to all functions
###                 before their respective 'return()' in
###                 order to try to clear some RAM memory

sequence_characterizer <- function(sequence, k_max = 4, Ks,
                                   no_names = FALSE, as_vector = TRUE) {
  basepercs <- bases_percentage(sequence)
  seq_tms <- tm_len_mt13(sequence)
  seq_sha <- shannon_entropy(sequence)

  # Given that 'Biostrings' masks many functions from many packages,
  # I preferred for its' functions to be called via double colons ("::").
  seq <- Biostrings::DNAString(sequence)
  rev_seq <- Biostrings::DNAString(stri_reverse(sequence))
  revc_seq <- Biostrings::DNAString(rev_complement(sequence))

  rev_align <- Biostrings::pairwiseAlignment(seq, rev_seq, type = "global")
  ralign_id <- Biostrings::pid(rev_align)
  ralign_sc <- Biostrings::score(rev_align)

  revc_align <- Biostrings::pairwiseAlignment(seq, revc_seq, type = "global")
  rcalign_id <- Biostrings::pid(revc_align)
  rcalign_sc <- Biostrings::score(revc_align)

  if (missing(Ks))
    Ks <- 2:k_max

  kmer_diversities <- c()
  ## NOTE: Since kmer-diversities is way easier to compute, I'm tempted to
  ##       add another variable called maybe 'KDs' for kmers I only want to
  ##       have applied to kmer-diversity (or other single values obtained
  ##       kmer-wise).
  for (x in Ks) {
    # NOTE:Both kmer-shannon and absolute-diversity follow the same
    #      principle so its expected they will correlate later
    vars <- c("shan", "adiv", "rdiv")
    kd_names <- sapply(vars, function(var) paste0("k", x, "_", var))
    # Get Shannon Coefficient per kmers' percentages
    ki_sh <- kmers_shannon(sequence = sequence, k = x)
    # Get "absolut" kmer diversity
    ki_ad <- kmer_diversity(sequence = sequence, k = x)
    # Get "relative" kmer diversity
    ki_rd <- kmer_diversity(sequence = sequence, k = x, relative = TRUE)
    # ^relative to sequence length (in toher words, relative to
    #  the number of kmers possible given a sequence of n length)
    ki_diversities <- c(ki_sh, ki_ad, ki_rd)
    names(ki_diversities) <- kd_names
    kmer_diversities <- c(kmer_diversities, ki_diversities)
  }

  kmers_characs <- list()
  if (!as_vector) {
    kmers_characs[[1]] <- list(perc = basepercs, temp = seq_tms, shan = seq_sha,
                               r_id = ralign_id, r_sc = ralign_sc,
                               rc_id = rcalign_id, rc_sc = rcalign_sc,
                               k_div = kmer_diversities)
    names(kmers_characs)[1] <- "WS"
  } else {
    kiv_characs <- c(basepercs, temp = seq_tms, shan = seq_sha,
                     r_id = ralign_id, r_sc = ralign_sc,
                     rc_id = rcalign_id, rc_sc = rcalign_sc,
                     kmer_diversities)
  }

  for (i in Ks) {
    all_ki <- combi_kmers(k = i)
    ki_seq <- kmer_windows(sequence, k = i)
    ki_name <- paste("k", i, sep = "")

    ki_prod <- func_per_windows(kmers = all_ki, windows = ki_seq,
                                func = ksg_product)
    ki_bcds <- func_per_windows(kmers = all_ki, windows = ki_seq,
                                func = kmer_barcode, opt = "divexpsum")
    ki_bclp <- func_per_windows(kmers = all_ki, windows = ki_seq,
                                func = kmer_barcode, opt = "logprimes")
    ki_characs <- list(prod = ki_prod, bcds = ki_bcds, bclp = ki_bclp)
    # ^Another optimization is to only compute data from non-zero count kmers
    # and later on, organize it in the corresponding 'all_ki' order
    # Maybe I can just add an "isTrue" parameter/conditional to all functions

    if (as_vector) {
      for (kmer_i in 1:4^i) {
        ikmer_characs <- c(unlist(lapply(ki_characs, `[[`, kmer_i)))
        if (!no_names)
          names(ikmer_characs) <- paste("k", i, "-", kmer_i,
                                        "_", names(ikmer_characs), sep = "")
        kiv_characs <- c(kiv_characs, ikmer_characs)
      }
    } else {
      kmers_characs[[i]] <- ki_characs
      names(kmers_characs)[i] <- ki_name
    }
  }
  if (as_vector) {
    return(kiv_characs)
  } else {
    return(kmers_characs)
  }
}

sequences_characterizer <- function(sequences, k_max, Ks, optim = TRUE,
                                    as_df = TRUE, progressbar = TRUE) {
  # seqs_data <- list(A_perc = c(), T_perc = c(), C_perc = c(), G_perc = c(),
  #                   seq_TM = c(), seq_SE = c())
  if (progressbar)
    pb <- txtProgressBar(min = 0, max = length(sequences), initial = 0)

  if (missing(Ks))
    Ks <- 2:k_max

  # NOTE: Need to find a more optimal way to find 'ncols'
  # Number of 'kmer' characteristics' columns
  sc_len <- ifelse(optim,
                   #: cols_numb(k = k_max, n_params = 4),
                   cols_numb(Ks = Ks, n_params = 3),
                   cols_numb(Ks = Ks, n_params = 6))
  # Number of 'kmer diversity' columns
  n_kd_cols <- 3 * length(Ks)
  # Number of 'whole-sequence' characteristics' columns
  n_ws_cols <- 10 + n_kd_cols

  sc_len <- sc_len + n_ws_cols
  seqs_data <- vector(mode = "list", length = sc_len)
  i <- 0
  for (sequence_i in sequences) {
    if (!i) {
      seq_i <- sequence_characterizer(sequence_i, Ks = Ks,
                                      as_vector = TRUE)
                                      #: optim = optim, as_vector = TRUE)
      names(seqs_data) <- names(seq_i)
      names(seq_i) <- NULL
    } else {
      seq_i <- sequence_characterizer(sequence_i, Ks = Ks,
                                      no_names = TRUE, as_vector = TRUE)
                                      #: optim = optim, as_vector = TRUE)
      names(seq_i) <- NULL
    }
    for (j in 1:sc_len) {
      seqs_data[[j]] <- append(seqs_data[[j]], seq_i[j])
    }
    i <- i + 1
    if (progressbar)
      setTxtProgressBar(pb, i)
  }
  if (as_df) {
    return(as.data.frame(seqs_data))
  } else {
    return(seqs_data)
  }
}
cols_numb <- function(k_max, Ks, n_params) {
  number <- 0
  if (missing(Ks))
    Ks <- 2:k_max
  for (i in Ks) {
    number <- sum(number, 4^i)
  }
  return(number * n_params)
}

feat_abbrev_convert <- function(input) {
  fullname <- c(
    "whole sequence", "nucleotide percentages", "melting temperature",
    "shannon entropy coefficient of nucleotides", "global self-alignment",
    "reverse sequence", "reverse alignment identity", "reverse alignment score",
    "reverse complement sequence", "reverse complement alignment identity",
    "reverse complement alignment score", "kmer set", "kmer diversity",
    "shannon entropy coefficient of kmers", "absolute kmer diversity",
    "relative kmer diversity", "each kmer", "count percentage product",
    "barcode profile division-sum", "barcode profile log-primes"
  )

  abbreviations <- c(
    'WS', 'nb_perc', 'nb_temp', 'nb_shan', 'SA', 'RA', 'ralign_id', 'ralign_sc',
    'RCA', 'rcalign_id', 'rcalign_sc', 'KS', 'KD', 'kx_sh', 'kx_ad', 'kx_rd',
    'EK', 'ki_prod', 'ki_bcds', 'ki_bclp'
  )

  # Help option
  if (length(input) == 1 && (input == "help" || input == "?")) {
    help_df <- data.frame(
      FullName = fullname,
      Abbreviation = abbreviations,
      stringsAsFactors = FALSE
    )
    cat("Available features and their abbreviations:\n")
    print(help_df, row.names = FALSE)
    return(invisible(help_df))
  }

  # Create a combined lookup for fullnames and abbreviations
  combined_lookup <- c(
    setNames(abbreviations, fullname),
    setNames(abbreviations, abbreviations)
  )

  # Convert input to abbreviations
  result <- combined_lookup[input]

  # Ensure all elements are valid (per problem constraints)
  if (any(is.na(result))) {
    invalid <- input[is.na(result)]
    stop("\n Invalid input feature(s): ", paste(invalid, collapse = ", "),
         "\n To display help run: feat_abbrev_convert('?')")
  }

  return(unique(unname(result)))
}

feat_mappings <- c(
  WS = "Features per whole sequence:",
  nb_perc = "nb_perc: Nucleotide Percentages",
  nb_temp = "nb_temp: Melting Temperature",
  nb_shan = "nb_shan: Shannon entropy coefficient of Nucleotides",
  SA = "Global Self-Alignment",
  RA = "Sequence Reverse Self-Alignment",
  r_id = "r_id: reverse alignment identity",
  r_sc = "r_sc: reverse alignment score",
  RCA = "Sequence Reverse Complement Self-Alignment",
  rc_id = "rc_id: reverse complement alignment identity",
  rc_sc = "rc_sc: reverse complement alignment score",
  KS = "Features per kmer set:",
  KD = "Kmer Diversity",
  kx_shan = "kx_shan: Shannon entropy coefficient of Kmers",
  kx_adiv = "kx_adiv: Absolute kmer diversity",
  kx_rdiv = "kx_rdiv: Relative kmer diversity",
  EK = "Feature per each kmer:",
  ki_prod = "ki_prod: Count percentage product",
  ki_bcds = "ki_bcds: Barcode profile - Division of Exponent Sum",
  ki_bclp = "ki_bclp: Barcode profile - Log of Primes Product"
)

sequence_characterizer_ <- function(sequences, k_max, Ks, optim = TRUE,
                                    as_df = TRUE, progressbar = TRUE) {
  # seqs_data <- list(A_perc = c(), T_perc = c(), C_perc = c(), G_perc = c(),
  #                   seq_TM = c(), seq_SE = c())
  if (progressbar)
    pb <- txtProgressBar(min = 0, max = length(sequences), initial = 0)

  if (missing(Ks))
    Ks <- 2:k_max

  # NOTE: Need to find a more optimal way to find 'ncols'
  # Number of 'kmer' characteristics' columns
  sc_len <- ifelse(optim,
                   #: cols_numb(k = k_max, n_params = 4),
                   cols_numb(Ks = Ks, n_params = 3),
                   cols_numb(Ks = Ks, n_params = 6))
  # Number of 'kmer diversity' columns
  n_kd_cols <- 3 * length(Ks)
  # Number of 'whole-sequence' characteristics' columns
  n_ws_cols <- 10 + n_kd_cols

  sc_len <- sc_len + n_ws_cols
  seqs_data <- vector(mode = "list", length = sc_len)
  i <- 0
  for (sequence_i in sequences) {
    if (!i) {
      seq_i <- sequence_characterizer(sequence_i, Ks = Ks,
                                      as_vector = TRUE)
                                      #: optim = optim, as_vector = TRUE)
      names(seqs_data) <- names(seq_i)
      names(seq_i) <- NULL
    } else {
      seq_i <- sequence_characterizer(sequence_i, Ks = Ks,
                                      no_names = TRUE, as_vector = TRUE)
                                      #: optim = optim, as_vector = TRUE)
      names(seq_i) <- NULL
    }
    for (j in 1:sc_len) {
      seqs_data[[j]] <- c(seqs_data[[j]], seq_i[j])
    }
    i <- i + 1
    if (progressbar)
      setTxtProgressBar(pb, i)
  }
  if (as_df) {
    return(as.data.frame(seqs_data))
  } else {
    return(seqs_data)
  }
}

sequence_characterizer_by_feature <- function(sequences, k_max,
                                              Ks = c(2, 3),
                                              features_selected =
                                                c("nb_perc", "nb_temp",
                                                  "nb_shan", "RA", "RCA",
                                                  "KS", "EK"),
                                              print_feature_tree = TRUE,
                                              print_as_txt = FALSE,
                                              # progressbar = TRUE,
                                              as_df = FALSE) {
  # if (progressbar)
  #   pb <- txtProgressBar(min = 0, max = length(sequences), initial = 0)

  fs_keys <- feat_abbrev_convert(features_selected)
  feat_funcs <-
    list(WS = list(nb_perc = "bases_percentage",
                   nb_temp = "tm_len_mt13",
                   nb_shan = "shannon_entropy",
                   SA = list(RA = list(r_id = "pid r_align",
                                       r_sc = "score r_align"),
                             RCA = list(rc_id = "pid rc_align",
                                        rc_sc = "score rc_align"),
                             ID =  list(r_id = "pid r_align",
                                        rc_id = "pid rc_align"),
                             SC = list(r_sc = "score r_align",
                                       rc_sc = "score rc_align"))),
         KS = list(KD = list(kx_shan = "kmers_shannon",
                             kx_adiv = "kmer_abs_diversity",
                             kx_rdiv = "kmer_rel_diversity")),
         EK = list(ki_prod = "ksg_product no_option",
                   ki_bcds = "kmer_barcode divexpsum",
                   ki_bclp = "kmer_barcode logprimes"))
                   # ki_bclp = "kmer_barcode logprimes"),
         # RK = list(dna_shape = "dnashape_characterizer",
         #           RKS = list(ki_sp_pol = "kmer_spc_polynome",
         #                      ki_sp_lag = "kmer_spc_lagrange")))

  func_sublist <- retrieve_sublist(feat_funcs, fs_keys)
  func_vector <- retrieve_leaves_with_names(func_sublist)

  if (print_feature_tree) {
    substitute_names <- import_from_script("substitute_names",
                                           "scripts/custom-functions.R")
    print_tree_names <- import_from_script("print_tree_names",
                                           "scripts/custom-functions.R")
    # feat_tree <- print_tree_names(
    # substitute_names(func_sublist, feat_mappings))
    if (print_as_txt) {
      txt_name <- paste0("run_", format(Sys.time(), "%F_%H-%M"), ".txt")
      print_tree_names(substitute_names(func_sublist, feat_mappings),
                       file_name = txt_name)
    } else {
      print_tree_names(substitute_names(func_sublist, feat_mappings))
    }
  }

  ws_keys <- c("nb_perc", "nb_temp", "nb_shan",
               "r_id", "r_sc",
               "rc_id", "rc_sc")
  ks_keys <- c("kx_shan", "kx_adiv", "kx_rdiv")
  ek_keys <- c("ki_prod", "ki_bcds", "ki_bclp")

  ws_inds <- which(ws_keys %in% names(func_vector))
  ks_inds <- which(ks_keys %in% names(func_vector))
  ek_inds <- which(ek_keys %in% names(func_vector))

  ws_feats <- list()
  ks_feats <- list()
  ek_feats <- list()

  # Extraction of Whole-Sequence Features
  if (length(ws_inds)) {
    ws_funcs <- func_vector[names(func_vector) %in% ws_keys]
    ga_inds <- which(ws_keys %in% rev(ws_keys)[1:4])
    ga_eval <- ws_inds >= ga_inds[1]
    ga_feats <- list() ##

    if (any(ga_eval)) {
      ga_funcs <- ws_funcs[ws_keys[ws_inds[ga_eval]]]
      ws_funcs <- ws_funcs[!names(ws_funcs) %in% names(ga_funcs)]

      ga_feats <- vector(mode = "list", length = length(ga_inds))
      names(ga_feats) <- names(ga_funcs)

      for (sequence in sequences) {
        i <- 0
        if (any(ws_inds >= ga_inds[3])) {
          r_align <- Biostrings::pairwiseAlignment(sequence,
                                                   stri_reverse(sequence),
                                                   type = "global")
        }
        if (any(ws_inds[ga_eval] <= ga_inds[2])) {
          rc_align <- Biostrings::pairwiseAlignment(sequence,
                                                    rev_complement(sequence),
                                                    type = "global")
        }

        for (ga_func in ga_funcs) {
          i <- i + 1
          func_strs <- unlist(strsplit(ga_func, " "))
          ga_func <- func_strs[1]
          alignment <- func_strs[2]
          func_value <- get(ga_func,
                            envir = asNamespace("Biostrings"))(get(alignment))
          ga_feats[[i]] <- c(ga_feats[[i]], func_value)
        }
      }
    }

    if (1 %in% ws_inds) {
      ws_feats <- vector(mode = "list", length = length(ws_funcs) + 3)
      names(ws_feats) <- c("nb_perc.A", "nb_perc.T", "nb_perc.C", "nb_perc.G",
                           names(ws_funcs)[-1])
      for (sequence in sequences) {
        func_i <- 0
        for (ws_func in ws_funcs) {
          func_i <- func_i + 1
          func_value <- get(ws_func)(sequence)
          if (func_i > 1) {
            ws_feats[[func_i]] <- c(ws_feats[[func_i]], func_value)
          } else {
            for (nb in (1:4)) {
              ws_feats[[nb]] <- c(ws_feats[[nb]], func_value[nb])
            }
            func_i <- nb
          }
        }
      }
    } else {
      ws_feats <- vector(mode = "list", length = length(ws_funcs))
      names(ws_feats) <- names(ws_funcs)
      for (sequence in sequences) {
        func_i <- 0
        for (ws_func in ws_funcs) {
          func_i <- func_i + 1
          func_value <- get(ws_func)(sequence)
          ws_feats[[func_i]] <- c(ws_feats[[func_i]], func_value)
        }
      }
    }

    ws_feats <- c(ws_feats, ga_feats) ##
  }

  # When 'kmer' features are requested and Ks is missing. Set Ks
  if ((length(ks_inds) || length(ek_inds))) {
    if (!length(Ks)) {
      if (missing(k_max)) {
        kmer_keys_requested <-
          names(func_vector) %in% ks_keys | names(func_vector) %in% ek_keys
        stop("\n Requested features (", names(func_vector)[kmer_keys_requested],
             ") per 'kmer set' or per 'each kmer', but",
             "\n no 'Ks' or 'k_max' arguments provided.")
      }
      Ks <- 2:k_max
    }

    if (length(ks_inds)) {
      ks_funcs <- func_vector[names(func_vector) %in% ks_keys]
      ks_feats <- vector(mode = "list",
                         length = length(ks_inds) * length(Ks))
      k1_i <- 0
    }

    if (length(ek_inds)) {
      ek_funcs <- func_vector[names(func_vector) %in% ek_keys]
      ek_feats <- vector(mode = "list",
                         length = length(ek_inds) * sum(4^Ks))
      k2_i <- 0
    }

    for (k in Ks) {
      all_ki <- combi_kmers(k = k)
      seqs_ki <- lapply(sequences, function(sequence) kmer_windows(sequence, k))

      if (length(ks_inds)) {
        for (seq_ki in seqs_ki) {
          func_i <- 0
          for (ks_func in ks_funcs) {
            func_i <- func_i + 1
            ks_i <- k1_i + func_i
            ks_feats[[ks_i]] <- c(ks_feats[[ks_i]],
                                  get(ks_func)(seq_kmers = seq_ki,
                                               all_kmers = all_ki))
            names(ks_feats)[ks_i] <- sub("x", k, names(ks_funcs)[func_i])
          }
        }
        k1_i <- k1_i + length(ks_funcs)
      }

      if (length(ek_inds)) {
        for (seq_ki in seqs_ki) {
          func_i <- 0
          kset_i <- 0
          for (ek_func in ek_funcs) {
            func_strs <- unlist(strsplit(ek_func, " "))
            func_option <- func_strs[2]
            ek_func <- func_strs[1]
            func_i <- func_i + 1
            kmer_i <- 0
            for (each_ki in all_ki) {
              kmer_i <- kmer_i + 1
              ek_i <- k2_i + kset_i + kmer_i
              ek_feats[[ek_i]] <- c(ek_feats[[ek_i]],
                                    get(ek_func)(windows = seq_ki,
                                                 kmer = each_ki,
                                                 option = func_option))
              names(ek_feats)[ek_i] <- sub("i", paste0(k, ".", kmer_i),
                                           names(ek_funcs)[func_i])
            }
            kset_i <- kset_i + kmer_i
          }
        }
        k2_i <- k2_i + (length(ek_funcs) * 4^k)
      }
    }
  }
  seqs_data <- c(ws_feats, ks_feats, ek_feats)
  if (as_df) return(as.data.frame(seqs_data)) else return(seqs_data)
}

dnashape_characterizer <- function(fasta_path, shape_feats, as_df = TRUE,
                                   feat_ncols = 10, round_number = 3) {
  # All possible shape features:
  # -(defalut) MGW
  # -(defalut) HelT
  # -(defalut) ProT
  # -(defalut) Roll
  # -(defalut) EP
  # - Stretch
  # - Tilt
  # - Buckle
  # - Shear
  if (missing(shape_feats)) {
    seqs_shapes <- DNAshapeR::getShape(fasta_path)
  } else {
    seqs_shapes <- DNAshapeR::getShape(fasta_path, shapeType = shape_feats)
  }
  sf_list <- list()
  for (i_sf in seq_along(seqs_shapes)) {
    featname_i <- names(seqs_shapes)[i_sf]
    feat_i <- asplit(apply(seqs_shapes[[i_sf]], 1,
                           function(shape_vector) {
                             denoise_signal_fft(na.omit(shape_vector),
                                                N = feat_ncols,
                                                round_int = round_number)
                           }), MARGIN = 1)
    i_name <- paste0("shape", featname_i)
    names(feat_i) <- paste(i_name,
                           c("main", paste0("sub", (1:(feat_ncols - 1)))),
                           sep = ".")
    sf_list <- c(sf_list, feat_i)
  }

  if (as_df) return(as.data.frame(sf_list)) else return(sf_list)
}

# =============================================================================
# REGRESSIONS: Polynomial, Fourier-like and Bin-like modelling for kmer data  #
# =============================================================================

# Need to-do: modelling in a bin-like (histogram) way
kmer_spc_polynome <- function(sequence, k, windows, kmer,
                              plot = FALSE, option = "norm") {
  if (missing(windows)) {
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
tvalues <- function(summ_lmodel, mean = FALSE, prod = FALSE) {
  # Taken from 'trash-code.R'
  tvals <- summ_lmodel$coefficients[, "t value"]
  if (mean) {
    return(mean(tvals))
  }
  if (prod) {
    return(prod(tvals))
  }
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

#: model_peaks <- function(lens, poss, peaks) {}
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

# To do later:
# - function that gets the relative number of peaks by
#   looking non-increasing values in the response variable
#   (after removing ~80% of lower-end points).
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
    print(t(coeffs_storage))
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

denoise_signal_fft <- function(signal, N = 5, Fs = 500, relative_Fs = FALSE,
                               predict = FALSE, plot = FALSE, round_int = 3) {
  # Normalize signal to [-1, 1]
  norm_signal <- 2 * (signal - min(signal)) / (max(signal) - min(signal)) - 1
  if (relative_Fs) Fs <- length(norm_signal)
  L <- length(norm_signal)
  t <- seq(0, (L - 1)) / Fs

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

  if (plot) predict <- TRUE
  if (predict) {
    # Build filtered FFT with only top N frequencies and symmetric components
    filtered_fft <- rep(0+0i, length(fft_result))
    keep_indices <- c(1, top_freq_indices, L - top_freq_indices + 2)
    filtered_fft[keep_indices] <- fft_result[keep_indices]

    # Inverse FFT to get denoised signal
    denoised_signal <- Re(fft(filtered_fft, inverse = TRUE) / L)

    if (plot) {
      # Plotting
      plot(t, norm_signal,
           type = "l", col = "deepskyblue",
           ylab = "Amplitude", xlab = "Time",
           main = "Original Noisy Signal vs Denoised Signal", lwd = 2)
      lines(t, denoised_signal, type = "l", col = "magenta", lwd = 2, lty = 2)
      legend("topright", legend = c("Noisy Signal", "Denoised Signal"),
             col = c("deepskyblue", "magenta"), bg = "white", pch = 16,
             horiz = TRUE, cex = 0.9)
    }

    return(invisible(list(
      top_frequencies_hz = round(top_freqs, round_int),
      denoised_signal = denoised_signal,
      normalized_signal = norm_signal,
      time_vector = t
    )))
  }

  # Return everything
  top_frequencies <- round(top_freqs, round_int)
  return(c(top_frequencies[1], sort(top_frequencies[-1], decreasing = TRUE)))
}

# =============================================================================
# UTILITIES: Tools used in some functions along this file                     #
# =============================================================================

required_libs <- function(libs) {
  for (lib in libs)
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

import_from_script <- function(function_name, script_path) {
  script_env <- new.env()
  script_content <- readLines(script_path)
  script_expr <- parse(text = script_content)
  eval(expr = script_expr, envir = script_env)
  return(get(function_name, envir = script_env))
}


# Function to retrieve a subset of the nested list based on acceptable keys.
retrieve_sublist <- function(x, keys) {
  # If x is not a list, there's nothing to prune.
  if (!is.list(x)) return(NULL)

  pruned <- list()
  for (nm in names(x)) {
    element <- x[[nm]]
    # If the current key is in the accepted keys, include it entirely.
    if (nm %in% keys) {
      pruned[[nm]] <- element
    } else if (is.list(element)) {
      # Otherwise, recursively prune the element.
      sub <- retrieve_sublist(element, keys)
      # Only include the element if the recursion returned a non-empty list.
      if (length(sub) > 0) {
        pruned[[nm]] <- sub
      }
    }
    # For atomic elements (non-list) not in keys, we do nothing.
  }
  pruned
}

# Function to recursively retrieve all the "leaves" from a nested list.
retrieve_leaves_with_names <- function(x, current_name = NULL) {
  if (!is.list(x)) {
    # x is atomic; if a current name exists, assign it to all elements
    n <- length(x)
    if (!is.null(current_name)) {
      names(x) <- rep(current_name, n)
    }
    return(x)
  }

  result <- c()
  nms <- names(x)

  for (i in seq_along(x)) {
    # If the current element has a name, use it; otherwise, use the parent's name.
    key <- if (!is.null(nms) && nms[i] != "") nms[i] else current_name
    result <- c(result, retrieve_leaves_with_names(x[[i]], key))
  }

  result
}

