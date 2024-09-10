setwd("/home/davidfm/Projects/UBMI-IFC/EnhaProm/")
source("scripts/genome-functions.R")

# ============================================================================ #
# ============================================================================ #
partial_palindromes_test1 <- function(windows, restrict_length, sequence, min_size,
                                only_indexes = FALSE, only_sequences = FALSE) {
  if (missing(sequence)) {
    if (missing(windows)) {
      stop("No 'sequence' or 'windows' parameters provided")
    }
    only_indexes <- TRUE
  } else if (missing(windows)) {
    windows <- kmer_windows(sequence, k = min_size)
  }
  k <- str_length(windows[1])
  vector_len <- length(windows)
  if (missing(restrict_length)) {
    restrict_length <- vector_len
  }
  rev_kmers <- stri_reverse(windows)
  if (!missing(only_indexes)) {
    return(pp_chunk_onlyInds(vector_len, k,
                             restrict_length, rev_kmers, windows))
  } else if (!missing(only_sequences)) {
    return(pp_chunk_onlySeqs(vector_len, k,
                             restrict_length, rev_kmers, windows, sequence))
  }  else {
    return(pp_chunk_Both_test1(vector_len, k,
                         restrict_length, rev_kmers, windows, sequence))
  }
}
pp_chunk_Both_test1 <- function(vector_len, kmer_size,
                          restrict_len, rev_kmers, win_kmers, sequence) {
  # pals_df <- data.frame(start = integer(), end = integer(),
                        # sequences = character())
  # loop_counter <- 0
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
      # df_new_rows <- data.frame(start = start_inds, end = end_inds,
                                # sequences = pal_seqs)
      # pals_df <- bind_rows(pals_df, df_new_rows)
      s_inds_row <- c(s_inds_row, start_inds)
      e_inds_row <- c(e_inds_row, end_inds)
      seqs_row <- c(seqs_row, pal_seqs)
      # loop_counter <- loop_counter + 1
    }
  }
  pals_df <- list(start = s_inds_row, end = e_inds_row,
                        sequences = seqs_row)
  return(pals_df)
}

# ============================================================================ #
# ============================================================================ #
partial_palindromes_test2 <- function(windows, restrict_length, sequence, min_size,
                                only_indexes = FALSE, only_sequences = FALSE) {
  if (missing(sequence)) {
    if (missing(windows)) {
      stop("No 'sequence' or 'windows' parameters provided")
    }
    only_indexes <- TRUE
  } else if (missing(windows)) {
    windows <- kmer_windows(sequence, k = min_size)
  }
  k <- str_length(windows[1])
  vector_len <- length(windows)
  if (missing(restrict_length)) {
    restrict_length <- vector_len
  }
  rev_kmers <- stri_reverse(windows)
  if (!missing(only_indexes)) {
    return(pp_chunk_onlyInds(vector_len, k,
                             restrict_length, rev_kmers, windows))
  } else if (!missing(only_sequences)) {
    return(pp_chunk_onlySeqs(vector_len, k,
                             restrict_length, rev_kmers, windows, sequence))
  }  else {
    return(pp_chunk_Both_test2(vector_len, k,
                         restrict_length, rev_kmers, windows, sequence))
  }
}
pp_chunk_Both_test2 <- function(vector_len, kmer_size,
                          restrict_len, rev_kmers, win_kmers, sequence) {
  # pals_df <- data.frame(start = integer(), end = integer(),
                        # sequences = character())
  # loop_counter <- 0
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
      # df_new_rows <- data.frame(start = start_inds, end = end_inds,
                                # sequences = pal_seqs)
      # pals_df <- bind_rows(pals_df, df_new_rows)
      s_inds_row <- c(s_inds_row, start_inds)
      e_inds_row <- c(e_inds_row, end_inds)
      seqs_row <- c(seqs_row, pal_seqs)
      # loop_counter <- loop_counter + 1
    }
  }
  pals_df <- data.frame(start = s_inds_row, end = e_inds_row,
                        sequences = seqs_row)
  return(pals_df)
}

# ============================================================================ #
# ============================================================================ #
partial_palindromes_test3 <- function(windows, restrict_length, sequence, min_size,
                                only_indexes = FALSE, only_sequences = FALSE) {
  if (missing(sequence)) {
    if (missing(windows)) {
      stop("No 'sequence' or 'windows' parameters provided")
    }
    only_indexes <- TRUE
  } else if (missing(windows)) {
    windows <- kmer_windows(sequence, k = min_size)
  }
  k <- str_length(windows[1])
  vector_len <- length(windows)
  if (missing(restrict_length)) {
    restrict_length <- vector_len
  }
  rev_kmers <- stri_reverse(windows)
  if (!missing(only_indexes)) {
    return(pp_chunk_onlyInds(vector_len, k,
                             restrict_length, rev_kmers, windows))
  } else if (!missing(only_sequences)) {
    return(pp_chunk_onlySeqs(vector_len, k,
                             restrict_length, rev_kmers, windows, sequence))
  }  else {
    return(pp_chunk_Both_test3(vector_len, k,
                         restrict_length, rev_kmers, windows, sequence))
  }
}
pp_chunk_Both_test3 <- function(vector_len, kmer_size,
                          restrict_len, rev_kmers, win_kmers, sequence) {
  # pals_df <- data.frame(start = integer(), end = integer(),
                        # sequences = character())
  # loop_counter <- 0
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
      # df_new_rows <- data.frame(start = start_inds, end = end_inds,
                                # sequences = pal_seqs)
      # pals_df <- bind_rows(pals_df, df_new_rows)
      s_inds_row <- c(s_inds_row, start_inds)
      e_inds_row <- c(e_inds_row, end_inds)
      seqs_row <- c(seqs_row, pal_seqs)
      # loop_counter <- loop_counter + 1
    }
  }
  pals_df <- data.frame(start = s_inds_row, end = e_inds_row,
                        sequences = seqs_row)
  return(pals_df)
}
