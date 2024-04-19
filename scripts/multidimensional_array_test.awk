#!/bin/awk
#
# Simple AWK script made in order to practice matrices' manipulation.
# It rearranges a previously imputed matrix in a 90 degree angle.
#
# Example use:
# (suppose that you are located in:
#	~/Epigenomica/UBMI-IFC/EnhaProm-Project/scripts)
# awk -f multidimensional_array_test.awk test-files/input_matrix.txt

{
	if (max_nf < NF)
		max_nf = NF
	max_nr = NR
	for (x=1; x <= NF; x++)
		vector[x, NR] = $x
	}

END{
	for (x=1; x <= max_nf; x++){
		for (y=max_nr; y >= 1; --y){
			printf("%s ", vector[x,y])
			}
		printf("\n")
		}
	}
