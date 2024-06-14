library(stringi)
library(stringr)
library(primes)

# ============================================================================ #
#						  PARTIAL PALINDROMES TESTS
# ============================================================================ #

proto_partial_palindromes <- function(windows, restrict_in_kmers) {
  # This one just returns dataframe of indexes and is 30x times slower than
  # first implementation ('partial_pals1')
  k <- str_length(windows[1])
  vector_len <- length(windows)
  if (missing(restrict_in_kmers)) {
    restrict_len <- vector_len
  } else {
    restrict_len <- k * restrict_in_kmers
  }
  rev_kmers <- stri_reverse(windows)
  pal_inds_df <- data.frame(start = integer(), end = integer())
  for (i_start in (1:vector_len)) {
    i_end <- i_start - 1 + k + restrict_len
    if (i_end > vector_len) {
      i_end <- vector_len
    }
    pal_inds <- (i_start:i_end)[rev_kmers[i_start] == windows[i_start:i_end]]
    if (length(pal_inds)) {
      pal_inds <- pal_inds - 1 + k
      df_new_rows <- data.frame(start = i_start, end = pal_inds)
      pal_inds_df <- rbind(pal_inds_df, df_new_rows)
    }
  }
  return(pal_inds_df)
}

# Lets see if 'rbindlist()' makes the above code faster
# Note: 'rbindlist()' comes from the 'data.table' library
# which when loaded says:
# 'data.table 1.15.2 using 8 threads (see ?getDTthreads)'
# Later I'll have to check how and when to use more threads.
# Since 'rbindlist()' didn't perform as well as expected,
# I tried with 'bind_rows()' from the 'dplyr' library which
# performed better: now it performs just 10 TIMES SLOWER xD
# UPDATE: Now I tried with 'add_row()' function, also from
# 'dplyr', I had my hopes up, howver it doesnt perform as
# well as the previous approach with 'bind_rows()', not
# even using 'tibble()' instead of 'data.frame()' to create
# the dataframe
partial_palindromes <- function(windows, restrict_length, sequence, min_size,
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
    return(pp_chunk_Both(vector_len, k,
                         restrict_length, rev_kmers, windows, sequence))
  }
}
pp_chunk_Both <- function(vector_len, kmer_size,
                          restrict_len, rev_kmers, win_kmers, sequence) {
  pals_df <- data.frame(start = integer(), end = integer(),
                        sequences = character())
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
      df_new_rows <- data.frame(start = start_ind, end = end_inds,
                                sequences = pal_seqs)
      pals_df <- bind_rows(pals_df, df_new_rows)
    }
  }
  return(pals_df)
}
pp_chunk_onlyInds <- function(vector_len, kmer_size,
                              restrict_len, rev_kmers, win_kmers) {
  pals_df <- data.frame(start = integer(), end = integer())
  for (start_ind in (1:vector_len)) {
    end_ind <- start_ind + kmer_size + restrict_len
    if (end_ind > vector_len) {
      end_ind <- vector_len
    }
    end_inds <- (start_ind :end_ind)[rev_kmers[start_ind]
                                     == win_kmers[start_ind:end_ind]]
    if (length(end_inds)) {
      end_inds <- end_inds - 1 + kmer_size
      df_new_rows <- data.frame(start = start_ind, end = end_inds)
      pals_df <- bind_rows(pals_df, df_new_rows)
    }
  }
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

partial_palindromes_test <- function(windows, restrict_length, sequence, min_size,
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
    return(pp_chunk_Both_test(vector_len, k,
                         restrict_length, rev_kmers, windows, sequence))
  }
}
pp_chunk_Both_test <- function(vector_len, kmer_size,
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

# Originally called 'partialPalindromes_ALLindexes1()'
pp_df_all1 <- function(sequence, k = 2, s = 1, windows) {
  # Implementation of 'partial_pals1()' using dataframe
  # like in 'proto_partial_palindromes()'
  if (missing(windows)) {
    windows <- kmer_windows(sequence, k = k, s = s)
  } else if (missing(sequence)) {
    k <- str_length(windows[1])
  }
  len_vec <- length(windows)
  rev_kmers <- stri_reverse(windows)
  pal_inds_df <- data.frame(start = integer(), end = integer())
  for (i in (1:len_vec)) {
    pal_inds <- (i:len_vec)[rev_kmers[i] == windows[i:len_vec]]
    if (length(pal_inds)) {
      pal_inds <- pal_inds - 1 + k
      df_new_rows <- data.frame(start = i, end = pal_inds)
      pal_inds_df <- rbind(pal_inds_df, df_new_rows)
    }
  }
  return(pal_inds_df)
}

# Originally called 'partialPalindromes_ALLindexes2()'
pp_df_all2 <- function(sequence, k = 2, s = 1, win_kmers) {
  if (missing(win_kmers)) {
    win_kmers <- kmer_windows(sequence, k = k, s = s)
  } else if (missing(sequence)) {
    k <- str_length(win_kmers[1])
  }
  len_vec <- length(win_kmers)
  rem_inds <- (1:len_vec)
  rev_kmers <- stri_reverse(win_kmers)
  pal_inds_df <- data.frame(start = integer(), end = integer())
  for (i in (1:len_vec)) {
    pal_inds <- rem_inds[rev_kmers[i] == win_kmers]
    if (length(pal_inds)) {
      pal_inds <- pal_inds - 1 + k
      df_new_rows <- data.frame(start = i, end = pal_inds)
      pal_inds_df <- rbind(pal_inds_df, df_new_rows)
    }
    win_kmers <- win_kmers[-1]
    rem_inds <- rem_inds[-1]
  }
  return(pal_inds_df)
}
# Both functions seem to perform in almost the same time
# ('partialPalindromes_ALLindexes1()' and 'partialPalindromes_ALLindexes2()')
# although the first one seems a little faster. By practical means, differences
# between them are negligible but first is preferable since uses one less vector

# Originally called 'partial_pals1()'
pp_list_all1 <- function(sequence, win_kmers) {
  partial_pals <- list()
  len_vec <- length(win_kmers)
  rev_kmers <- stri_reverse(win_kmers)
  for (i in (1:len_vec)) {
    pal_index <- (i:len_vec)[rev_kmers[i] == win_kmers[i:len_vec]]
    pals_len <- length(partial_pals)
    if (length(pal_index)) {
      partial_pals[[pals_len + 1]] <- pal_index
    } else {
      partial_pals <- append(partial_pals, NA)
    }
  }
  if (missing(sequence)) {
    return(partial_pals)
  } else {
    len_pp <- length(partial_pals)
    partial_pals_seqs <- list()
    for (i in (1:len_pp)) {
      if (!is.na(partial_pals[[i]][1])) {
        temp_vect <- c()
        for (indx in partial_pals[[i]]) {
          temp_vect <- append(temp_vect, str_sub(sequence, i, indx + 3))
        }
        partial_pals_seqs[[i]] <- temp_vect
      } else {
        partial_pals_seqs <- append(partial_pals_seqs, NA)
      }
    }
    return(partial_pals_seqs)
  }
}

# Originally called 'partial_pals2()'
pp_list_all2 <- function(sequence, win_kmers) {
  partial_pals <- list()
  len_vec <- length(win_kmers)
  rev_kmers <- stri_reverse(win_kmers)
  rem_kmers <- win_kmers
  rem_inds <- (1:len_vec)
  for (i in (1:len_vec)) {
    pal_index <- rem_inds[rev_kmers[i] == rem_kmers]
    pals_len <- length(partial_pals)
    if (length(pal_index)) {
      partial_pals[[pals_len + 1]] <- pal_index
    } else {
      partial_pals <- append(partial_pals, NA)
    }
    rem_kmers <- rem_kmers[-1]
    rem_inds <- rem_inds[-1]
  }
  if (missing(sequence)) {
    return(partial_pals)
  } else {
    len_pp <- length(partial_pals)
    partial_pals_seqs <- list()
    for (i in (1:len_pp)) {
      if (!is.na(partial_pals[[i]][1])) {
        temp_vect <- c()
        for (indx in partial_pals[[i]]) {
          temp_vect <- append(temp_vect, str_sub(sequence, i, indx + 3))
        }
        partial_pals_seqs[[i]] <- temp_vect
      } else {
        partial_pals_seqs <- append(partial_pals_seqs, NA)
      }
    }
    return(partial_pals_seqs)
  }
}
# However, 'partial_pals' (which also have similar performance speeds respective
# to each other while performing the same task as partialPalindromes_ALLindexes
# pair) gives an 'uglier' output, they perform around 30 TIMES FASTER than
# 'partialPalindromes_ALLindexes'


# ============================================================================ #
#						  INTERSPACE POLYNOMES
# ============================================================================ #
kmer_interspace_polynome <- function(sequence, k, windows, kmer,
                                     plot = TRUE, fixed = FALSE,
                                     case, subcase, minus = 0,
                                     kdesmooth = FALSE, bw,
                                     silent = FALSE, brief = FALSE,
                                     best_adj = FALSE, which_inds = TRUE,
                                     rm_nth = 0,
                                     aem_option=5, aem_graph=FALSE) {
  k_len <- stri_length(windows[1])
  seq_len <- length(windows) + k_len - 1

  indexes <- which(windows == kmer)
  i_stts <- c(0, indexes + k_len - 1)
  i_ends <- c(indexes - 1, seq_len)

  if (which_inds) {
    print(i_stts)
    print(i_ends)
  }

  i_remv <- which(i_ends <= i_stts)
  if (length(i_remv)) {
    i_stts <- i_stts[-i_remv]
    i_ends <- i_ends[-i_remv]
  }
  if (rm_nth) {
    i_stts <- i_stts[-rm_nth]
    i_ends <- i_ends[-rm_nth]
  }

  # spc_pos <- log(generate_n_primes(seq_len)[i_stts + 1], 16)
  # spc_len <- 1.01^(i_ends - i_stts)
  spc_pos <- i_stts / seq_len
  # print(spc_pos)
  spc_len <- ((i_ends - i_stts) / seq_len)^2
  # print(spc_len)
  if (which_inds) {
    print(spc_pos)
    print(spc_len)
  }
  if (missing(case))
    case <- length(spc_pos) - 1
  if (missing(subcase))
    subcase <- 1

  # print(spc_pos)
  # print(spc_len)
  # model <- lm(spc_len ~ spc_pos + I(spc_pos^2) + I(spc_pos^3) + I(spc_pos^4))
  # model <- lm(spc_len ~ cos(spc_pos))
  # model <- lm(spc_len ~ cos(spc_pos / pi))
  # model <- lm(spc_len ~ cos(spc_pos / pi) + sin(spc_pos / pi))
  # print("==============================MODEL1==============================")
  # print(summary(model))
  # model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos))
  # print("==============================MODEL2==============================")
  # print(summary(model))
  # # model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos^2) + cos(spc_pos / pi))
  # model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos) + cos(spc_pos / pi))
  # print("==============================MODEL3==============================")
  # print(summary(model))
  # # model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos^2) + cos(spc_pos / pi) + sin(spc_pos / pi))
  # model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos) + cos(spc_pos / pi) + sin(spc_pos / pi))
  # print("==============================MODEL4==============================")
  # print(summary(model))
  # print(spc_pos, pi, cos(spc_pos))
  # model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos)
  #             + cos(spc_pos / pi) + sin(spc_pos / pi) + cos(pi / spc_pos))
  # print("==============================MODEL5==============================")
  # print(summary(model))
  # model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos)
  #             + cos(spc_pos / pi) + sin(spc_pos / pi)
  #             + cos(2 * pi / spc_pos) + sin(2 * pi / spc_pos))
  # print("==============================MODEL6==============================")
  # print(summary(model))
  # model <- lm(spc_len
              # ~ cos(2 * pi * spc_pos)
              # + sin(2 * pi * spc_pos)
              # ~ cos(spc_pos)
              # ~ sin(spc_pos)
              # + sin(spc_pos)
              # + cos(spc_pos)
              # + cos(4 * pi * spc_pos)
              # + sin(4 * pi * spc_pos)
              # + cos(spc_pos / (1.1^spc_pos))
              # + sin(spc_pos / (1.1^spc_pos))
              # + cos(6 * pi * spc_pos)
              # + sin(6 * pi * spc_pos)
              # + cos(spc_pos^2)
              # + sin(spc_pos^2)
              # + cos(8 * pi * spc_pos)
              # + sin(8 * pi * spc_pos)
              # + cos(spc_pos^3)
              # + sin(spc_pos^3)
              # + cos(spc_pos^4)
              # + sin(spc_pos^4)
              # + I(spc_pos)
              # + I(spc_pos^2)
              # + I(spc_pos^3)
              # )
  if (kdesmooth) {
    require('scales')
    # z <- ksmooth(spc_pos, spc_len, kernel = "normal", bandwidth = bw)
    print(bw)
    z <- density(spc_pos, kernel = "gaussian", bw = bw)
    print(z)
    # plot(spc_len ~ spc_pos, pch=16, cex=0.8, col=alpha("grey", 0.9))
    # lines(z, lwd=2, col=alpha("black", 0.9))
    plot(z, lwd=2, col=alpha("black", 0.9))
    detach("package:scales",unload=TRUE)
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
  # print("==============================MODEL6==============================")
  # print(summary(model))
  # model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos)
  #             + cos(spc_pos / pi) + sin(spc_pos / pi)
  #             + cos(pi / spc_pos) + sin(pi / spc_pos) + cos(2 * pi / spc_pos))
  # print("==============================MODEL7==============================")
  # print(summary(model))
  # model <- lm(spc_len ~ cos(spc_pos) + sin(spc_pos)
  #             + cos(spc_pos / pi) + sin(spc_pos / pi)
  #             + cos(pi / spc_pos) + sin(pi / spc_pos)
  #             + cos(2 * pi / spc_pos) + sin(2 * pi / spc_pos))
  # print("==============================MODEL8==============================")
  # print(summary(model))
  # model <- lm(spc_len ~ cos(spc_pos / pi) + sen(spc_pos / pi) + cos(spc_pos / pi))
  # print(model)
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
                  # lm(spc_len
                  #   ~ cos(4 * pi * spc_pos)
                  # )
                ),
                switch(subcase,
                  lm(spc_len
                    ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
                  ),
                  lm(spc_len
                    ~ cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  ),
                  lm(spc_len
                    # ~ cos(5 * pi * spc_pos) + sin(5 * pi * spc_pos)
                    ~ cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
                  ),
                  lm(spc_len
                    # ~ cos(1 * pi * spc_pos) + sin(5 * pi * spc_pos)
                    ~ cos(2 * pi * spc_pos) + sin(4 * pi * spc_pos)
                  ),
                  lm(spc_len
                    # ~ cos(2 * pi * spc_pos) + sin(5 * pi * spc_pos)
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
						 best_adj, 
						 silent, brief, 
						 aem_option, aem_graph) {
  # lm_prev <- trig_cases(spc_len, spc_pos, 1, 1)
  lm_prev <- lm(spc_len ~ cos(2 * pi * spc_pos))
  prev_aem <- adj_error_metric(lm_prev, best_adj, silent, brief, aem_option)
  # for (i in seq(2, sum(seq(1, case)) - 1)) {
  if (aem_graph) {
    aem_vec <- c()
    aem_cases <- c()
  }
  best_case <- c(1,1)
  for (i in seq(1, case - 1)) {
    # for (j in seq(1, i * 3)) {
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
        case_name <- paste(i, j, sep="-")
        aem_cases <- append(aem_cases, case_name)
      }
      # both_rsq <- lm_post$r.squared + lm_post$adj.r.squared
      # if (both_rsq == )
      if (post_aem > prev_aem) {
        if (!silent)
          cat("\n\t>>>> Better AEM inside Case:", i, "-", j)
        best_case <- c(i,j)
        lm_prev <- lm_post
        prev_aem <- adj_error_metric(lm_prev, best_adj, silent, brief, aem_option)
      }
    }
  }
  if (aem_graph) {
    aem_vec <- aem_vec[-1]
    aem_cases <- aem_cases[-1]
    aem_ids <- seq(1,length(aem_vec))
    aem_vec*100
    plot(aem_ids, aem_vec)
    axis(side=1, at=aem_ids, labels=aem_cases)
  }
  cat("\n Best Case chosen: ", best_case[1], "-", best_case[2], "\n")
  return(lm_prev)
}
adj_error_metric <- function(lmodel, best_adj=FALSE,
                             silent=FALSE, brief=FALSE, aem_option=5) {
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

  errs_names <- c("Adj. R-Squared", "Res. Std. Error", "R-Squared", "P-Value", "Mean T", "Prod T")
  errs_vectr <- c(adj_rsq, res_err, rsq, p_val, t_mean, t_prod)

  aem <- switch(aem_option,
    (adj_rsq * res_err) / p_val,
    (adj_rsq * rsq * res_err) / p_val, #******2
    (adj_rsq + res_err) * (rsq / p_val),
    (res_err / adj_rsq) * (rsq / p_val),
    # ((adj_rsq + rsq) ** adj_rsq) / (rsq * res_err * p_val),
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

# trigonom_cases <- function(lens, poss, case) {
#   return(switch(case,
#     lm(lens
#       ~ cos(2 * pi * spc_pos)
#     ),
#     lm(lens
#       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#     ),
#     lm(lens
#       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#       + cos(4 * pi * spc_pos)
#     ),
#     lm(lens
#       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#     ),
#     lm(lens
#       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#       + cos(6 * pi * spc_pos)
#     ),
#     lm(lens
#       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#       + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
#     ),
#     lm(lens
#       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#       + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
#       + cos(8 * pi * spc_pos)
#     ),
#     lm(lens
#       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#       + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
#       + cos(8 * pi * spc_pos) + sin(8 * pi * spc_pos)
#     ),
#     lm(lens
#       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#       + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
#       + cos(8 * pi * spc_pos) + sin(8 * pi * spc_pos)
#       + cos(10 * pi * spc_pos)
#     ),
#     lm(lens
#       ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)
#       + cos(4 * pi * spc_pos) + sin(4 * pi * spc_pos)
#       + cos(6 * pi * spc_pos) + sin(6 * pi * spc_pos)
#       + cos(8 * pi * spc_pos) + sin(8 * pi * spc_pos)
#       + cos(10 * pi * spc_pos) + sin(10 * pi * spc_pos)
#     ),
#     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#     lm(lens ~ cos(2 * pi * spc_pos) + sin(2 * pi * spc_pos)),
#   ))
# }
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
