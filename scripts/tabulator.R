# Execution time: Start
start_time <- Sys.time()

# Project path and additional functions
projpath <- "/home/davidfm/Projects/UBMI-IFC/EnhaProm/"
source(paste0(projpath, "scripts/genome-functions.R"))

# Required libraries
library(optparse)
library(doParallel)
library(foreach)

# Define command-line args
option_list <- list(
  make_option(c("--cre"), type = "character", default = "test",
              help = "Type of CRE: enhancers, promoters, OCRs",
              metavar = "character"),
  make_option(c("--n_sequences"), type = "integer", default = 60,
              help = "Number of sequences in the fasta file",
              metavar = "integer"),
  make_option(c("--n_elements_per_csv"), type = "integer", default = 5,
              help = "Number of elements per CSV file",
              metavar = "integer"),
  make_option(c("--n_cores"), type = "integer", default = 6,
              help = "Number of cores to use for parallel processing",
              metavar = "integer"))

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Assign variables from arguments
cre <- opt$cre
n_sequences <- opt$n_sequences
n_elements_per_csv <- opt$n_elements_per_csv
n_cores <- opt$n_cores

# Load sequences
fasta_path <- paste0(projpath, "datasets/GB-Test/", cre, "_training.fasta")
seqs <- scan(fasta_path, character(), quote = "")[seq(2, n_sequences * 2, 2)]

# Number of iterations
n <- n_sequences / n_elements_per_csv

# Set up parallel computing
corescluster <- makeCluster(n_cores)
registerDoParallel(corescluster)

required_libs <- required_libs ##
libs <- c("stringr", "stringi", "primes")

# Process sequences in parallel
foreach(i = 1:n) %dopar% {
  required_libs(libs)
  i_start <- ((i - 1) * n_elements_per_csv) + 1
  i_final <- i * n_elements_per_csv
  if (i > 1) {
    write.table(sequences_characterizer(seqs[i_start:i_final], Ks = c(2,3,5)),
                paste0(projpath, "datasets/GB-Test/", cre, "-training_", i, ".csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
  } else {
    write.csv(sequences_characterizer(seqs[i_start:i_final], Ks = c(2,3,5)),
              paste0(projpath, "datasets/GB-Test/", cre, "-training_", i, ".csv"),
              row.names = FALSE)
  }
}
stopCluster(corescluster)

# Execution time: end
end_time <- Sys.time()
exec_time <- start_time - end_time
timeFile <- file(paste0("exec_time-", cre, ".txt"))
writeLines(exec_time, timeFile)
close(timeFile)

# Print completion message
cat("Processing complete. ", cre, "'s CSV files saved. \n")


