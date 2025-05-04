#!/usr/bin/awk -f

BEGIN {
    FS = ","      # Input field separator (CSV format)
    OFS = ","     # Output field separator (CSV format)
	OFMT = "%.7g" # Same number of decimals as R
}

FNR == 1 {
    if (ARGIND == 1) {
        # Capture column headers from the first file
        for (i = 1; i <= NF; i++) {
            headers[i] = $i
        }
        # Print output headers
        print "Type", "ColName", "Mean", "StDev"
    }
    next  # Skip the header row of each file
}

{
    # Process data rows
    for (i = 1; i <= NF; i++) {
        sum[i] += $i         # Sum of each column
        sumsq[i] += $i * $i  # Sum of squares for each column
    }
}

ENDFILE {
    # Determine the type based on the filename
    type = (FILENAME ~ /promoter/) ? "Promoter" : \
           (FILENAME ~ /enhancer/) ? "Enhancer" : \
           (FILENAME ~ /ocr/) ? "OCR" : "Unknown"

	count = FNR - 1
    # Output summaries for each column
    for (i = 1; i <= NF; i++) {
        mean = sum[i] / count
        stddev = sqrt((sumsq[i] / count) - (mean * mean))
        print type, headers[i], mean, stddev
    }

    # Reset data for the next file
    delete sum
    delete sumsq
}
