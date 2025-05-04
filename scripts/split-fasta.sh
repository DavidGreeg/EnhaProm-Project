#!/bin/bash

# Usage: ./split_fasta.sh input.fasta N
input="$1"
N="$2"

# Create a temporary directory
tempdir=$(mktemp -d)

# Count the number of sequences
num_records=$(grep -c '^>' "$input")

# Calculate records per file (ceiling division)
records_per_file=$(( (num_records + N - 1) / N ))

# Split the FASTA into N files
awk -v rpf="$records_per_file" -v tempdir="$tempdir" '
  BEGIN { part=1; count=0; outfile=sprintf("%s/part_%d.fasta", tempdir, part) }
  /^>/ {
    if (count >= rpf) {
      close(outfile)
      part++
      outfile=sprintf("%s/part_%d.fasta", tempdir, part)
      count=0
    }
    count++
  }
  { print > outfile }
' "$input"

echo "Temporary files created in: $tempdir"
# Pass $tempdir to your R script here (e.g., Rscript script.R "$tempdir")
