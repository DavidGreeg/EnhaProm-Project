#!/bin/bash

# Function to generate .txt files with random DNA sequences
generate_dna_files() {
  local num_files=$1  # Number of files to generate
  local seq_length=$2 # Length of each DNA sequence
  local output_dir=$3 # Directory to store the generated files

  # Create the output directory if it doesn't exist
  mkdir -p "$output_dir"

  # DNA bases
  bases=("A" "T" "C" "G")

  # Generate files
  for i in $(seq 1 "$num_files"); do
    # Generate a random DNA sequence
    dna_sequence=""
    for _ in $(seq 1 "$seq_length"); do
      dna_sequence+="${bases[RANDOM % 4]}"
    done

    # Write the sequence to a .txt file
    echo "$dna_sequence" > "$output_dir/dna_sequence_${i}.txt"
  done

  echo "$num_files DNA sequence files generated in '$output_dir'."
}


# Example usage:
# generate_dna_files <number_of_files> <sequence_length> <output_directory>
#generate_dna_files 10 50 "./random_dna_files"
