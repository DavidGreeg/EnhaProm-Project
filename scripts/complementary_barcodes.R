pals_barcode <- function(sequence, k, windows, kmer,
						 bar_option = "simple", calc_option = "expsum") {
  if (missing(windows)) {
    windows <- kmer_windows(sequence, k = k)
  }
  if (! kmer %in% windows)
    return(0)

  n_kmers <- length(windows)

  positions <- which(windows == kmer)
  barcode_profile <- switch(calc_option,
    "expsum" = expsum_profile(n_kmers, positions),
    "primes" = primes_profile(n_kmers, positions),
    "logprimes" = logprimes_profile(n_kmers, positions),
    {
      stop("No valid calc_option selected")
    }
  )
  return(barcode_profile)
}
get_boolpositions <- function(bar_option, windows, kmer) {
  if (! bar_option %in% c("simple", "partialpal", "revcomp", "pp_rc"))
    stop("No valid bar_option selected")

  if (bar_option == "partialpal" && !isPalindrome(kmer)) {
      rev_kmer <- stri_reverse(kmer)
      if (rev_kmer %in% windows)) {
        return(which(windows == kmer && windows == rev_kmer))
      } else {
        return(which(windows == kmer))
      }
        
  } else if (bar_option == "revcomp" && !isRevComPalindrome(kmer)) {
      if (rev_complement(kmer) %in% windows) {  
    
  } else {
    return(which(windows == kmer))
  }
}
expsum_profile <- function(n_kmer, positions) {
  # : return(sum(1.001^(which(rev(windows == kmer)) - 1)))
  return(sum(1.001^positions - 1))
}
primes_profile <- function(n_kmers, positions) {
  return(prod(generate_n_primes(n_kmers)[positions]))
}
logprimes_profile <- function(n_kmers, positions) {
  return(prod(log(generate_n_primes(n_kmers), 8)[positions]))
}
