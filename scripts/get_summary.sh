#!/bin/bash

# Directory containing the CSV files
input_dir="GB-Testing"
output_file="summaries.csv"

# Initialize the output file with headers
echo "File,Type,Column,Mean,StDev" > $output_file

# Loop through CSV files
for csv_file in "$input_dir"/*.csv; do
  # Extract the base name of the file
  base_name=$(basename -s .csv "$csv_file")
  type_name=$(echo "$base_name" | sed 's/test-1638-//; s/s-6mers//')  # Extract type (e.g., Promoter, Enhancer, OCR)
  
  # Process each CSV file with awk
  awk -v type="$type_name" '
    BEGIN { OFS="," }
    NR>1 {
      for (i = 1; i <= NF; i++) {
        sum[i] += $i
        sumsq[i] += $i * $i
        count[i]++
      }
    }
    END {
      for (i = 1; i <= NF; i++) {
        mean = sum[i] / count[i]
        stddev = sqrt((sumsq[i] / count[i]) - (mean * mean))
        print type, "Column" i, mean, stddev
      }
    }
  ' FS="," "$csv_file" > $output_file
done
