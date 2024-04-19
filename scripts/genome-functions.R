library(stringr)
library(stringi)

# seq <- "gtatgggaatcagccgggtctcactatgtgcaaaggagattcggtcgtgtggtacttattcag"
seq <- "gtatgggaatcagccgggtctcactatgtgcaaa"
# seq <- "gtatgggaat"
sequence <- toupper(seq)
bases_vector <- c("A", "T", "C", "G")

bases_count <- function(sequence, bases = c("A", "T", "C", "G")) {
  base_counts <- str_count(sequence, bases)
  names(base_counts) <- bases
  return(base_counts)
}

bases_percentage <- function(sequence, bases = c("A", "T", "C", "G")) {
  base_percs <- bases_count(sequence, bases) / str_length(sequence)
  names(base_percs) <- bases
  return(base_percs)
}

gc_percentage <- function(sequence, bases = c("A", "T", "C", "G")) {
  base_percs <- bases_percentage(sequence, bases)
  gc_perc <- base_percs[3] + base_percs[4]
  names(gc_perc) <- "GC%"
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

combi_kmers <- function(bases = c("A", "T", "C", "G"), k = 2) {
  # 'sort' orders kmers list so that they follow the same order as the ones
  # transformed into factors in the (kinda homologous) function kmer_wincount()
  return(sort(do.call(paste0, expand.grid(rep(list(bases), k)))))
}

count_kmers <- function(sequence, vector_kmers, k = 2) {
  if (missing(vector_kmers)) {
    vector_kmers <- combi_kmers(k = k)
  }
  kmer_counts <- stri_count_fixed(sequence, vector_kmers, overlap = TRUE)
  names(kmer_counts) <- vector_kmers
  return(kmer_counts)
}

kmer_windws <- function(sequence, k = 2, s = 1) {
  # 'k' stands for 'kmer' (in this case its equivalent to the window size)
  # 's' stands for 'stride' (spaces taken between each window)
  seq_len <- str_length(sequence)
  return(str_sub(sequence, seq(1, seq_len + 1 - k, s), seq(k, seq_len, s)))
}

kmer_wincount <- function(sequence, k = 2, s = 1, win_kmers, all_kmers) {
  if (missing(win_kmers)) {
    win_kmers <- kmer_windws(sequence, k = k, s = s)
  }
  if (missing(all_kmers)) {
    all_kmers <- combi_kmers(k = k)
  }
  win_kmers <- factor(win_kmers)
  all_kmers <- factor(all_kmers)
  #=====================================================================#
  counts_per_window <- table(all_kmers[match(win_kmers, all_kmers)])
  return(counts_per_window)
}

tm_calc <- function(sequence, bases = c("A", "T", "C", "G"), base_counts) {
  if (missing(base_counts)) {
    base_counts <- bases_count(sequence, bases)
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
  names(temp) <- "TM in Celcius"
  return(temp)
}

tm_len_lt14 <- function(sequence, bases = c("A", "T", "C", "G"), base_counts) {
  if (missing(base_counts)) {
    base_counts <- bases_count(sequence, bases)
  }
  countA <- base_counts[1]
  countT <- base_counts[2]
  countC <- base_counts[3]
  countG <- base_counts[4]
  temp <- ((countA + countT) * 2) + ((countC + countG) * 4)
  names(temp) <- "TM in Celcius"
  return(temp)
}

tm_len_mt13 <- function(sequence, bases = c("A", "T", "C", "G"), base_counts) {
  if (missing(base_counts)) {
    base_counts <- bases_count(sequence, bases)
  }
  countA <- base_counts[1]
  countT <- base_counts[2]
  countC <- base_counts[3]
  countG <- base_counts[4]
  temp <- 64.9 + (41 * (countG + countC - 16.4) / (countA + countT + countC + countG))
  names(temp) <- "TM in Celcius"
  return(temp)
}




