source("/home/davidfm/Projects/UBMI-IFC/EnhaProm/scripts/genome-functions.R")
cre <- "enhancers"
cre_fasta <- paste("datasets/GenomicBenchmarks/", cre,
                   "_train_positive.fasta", sep = "")
seqs <- scan(cre_fasta, character(), quote = "")[seq(2, 29484, 2)]
library(doParallel)
library(foreach)
n <- 14742 / 819 #18
corescluster <- makeCluster(detectCores() - 6)
registerDoParallel(corescluster)
foreach(i = 1:n) %dopar% {
  library(fitdistrplus)
  library(stringr)
  library(stringi)
  library(primes)
  i_start <- ((i - 1) * 819) + 1
  i_final <- i * 819
  if (i > 1) {
    write.table(sequences_characterizer(promseqs[i_start:i_final],
                                        k_max = 6, optim = TRUE),
                paste("datasets/GB-Testing/", cre, "-training_",
                      i, ".csv", sep = ""), sep = ",",
                row.names = FALSE, col.names = FALSE)
  } else {
    write.csv(sequences_characterizer(promseqs[i_start:i_final],
                                      k_max = 6, optim = TRUE),
              paste("datasets/GB-Testing/", cre, "-training_",
                    i, ".csv", sep = ""), row.names = FALSE)
  }
}
stopCluster(corescluster)

# TEST
#: corescluster <- makeCluster(4)
#: registerDoParallel(corescluster)
#: foreach(i = 1:4) %dopar% {
#:   library(fitdistrplus)
#:   library(stringr)
#:   library(stringi)
#:   library(primes)
#:   i_start <- ((i - 1) * 819) + 1
#:   i_final <- i * 819
#:   if (i > 1)
#:     write.table(sequences_characterizer(promseqs[i_start:i_final], k_max = 4, optim = TRUE), 
#:                 paste("datasets/GB-Testing/promoters-training_", i, ".csv", sep = ""), 
#:                 sep = ",", row.names = FALSE, col.names = FALSE)
#:   else
#:     write.csv(sequences_characterizer(promseqs[i_start:i_final], k_max = 4, optim = TRUE),
#:               paste("datasets/GB-Testing/promoters-training_", i, ".csv", sep = ""), 
#:               row.names = FALSE)
