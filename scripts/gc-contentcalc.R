library(stringr)

# seq <- "gtatgggaatcagccgggtctcactatgtgcaaaggagattcggtcgtgtggtacttattcag"
seq <- "gtatgggaat"
sequence <- toupper(seq)
bases <- c("A", "T", "C", "G")
base_counts <- c()
base_percs <- c()

for (base in bases) {
	base_counts <- c(base_counts, str_count(sequence, base))
	base_percs <- c(base_percs, str_count(sequence, base) / str_length(sequence))
}

print(paste(sequence, "GC content: ", base_percs[3] + base_percs[4]))
# base_counts
# base_percs

