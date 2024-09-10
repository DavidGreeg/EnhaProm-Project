#!/bin/awk
#
# =Column Filter Tool=
#
# Example use:
# awk -v columns="6,7,8" -f scripts/column-filter.awk datasets/GeneHancer_RegElems_chr1.tsv
{
	split(columns,cols,",")
	i=1
	for(icol=1;icol<=NF;icol++){
		if (icol != cols[i]){
			printf "%s\t",$icol
		}
		else{
			i++
		}
	}
	printf "\n"
}
