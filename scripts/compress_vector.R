compress_vector <- function(vec, target_size, normalized = FALSE, tst = FALSE) {
  bin_size <- ceiling(length(vec) / target_size)
  if (!tst) {
    compressed_vec <- sapply(seq(1, length(vec), by = bin_size),
                             function(i) {
                               mean(vec[i:min(i + bin_size - 1, length(vec))])
                             })
  } else {
    compressed_vec <- sapply(seq(1, length(vec), by = bin_size),
                             function(i) {
                               sum(vec[i:min(i + bin_size - 1, length(vec))])
                             })
  }

  compressed_vec <- compressed_vec[1:target_size]
  if (normalized) {
    normalization_factor <- sum(vec) / sum(compressed_vec)
    compressed_vec_normal <- compressed_vec * normalization_factor
    return(compressed_vec_normal)
  } else {
    return(compressed_vec)
  }
}
