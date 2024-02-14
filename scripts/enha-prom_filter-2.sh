#!/bin/bash
#
# =Filter Tool=
# Made to retrieve enhancer and promoter fields from tables
# 
# Example of its use is:
# bash scripts/filter_enhpro.sh 10 RefSeq_FuncElems_chr1.tsv | column -s";" -t

# chrom="$1";
# start="$2";
# final="$3";

type="$1";
file="$2";

enha_prom_starts() {
awk 'BEGIN {FS="\t"}
{
	if ($'$type' ~ /nhancer/ || $'$type' ~ /romoter/)
	{
		print $2;
	}
}
' $file
}

function enha_prom_finals {
awk 'BEGIN {FS="\t"}
{
	if ($'$type' ~ /nhancer/ || $'$type' ~ /romoter/)
	{
		print $3;
	}
}
' $file
}

locations_starts="$(enha_prom_starts)"
locations_finals="$(enha_prom_finals)"

# size_list="$(echo $locations_starts | wc -l)"
# echo $size_list
echo $locations_starts
#| wc -l

# for start in $list_locations1
# do
# 	echo $start;
# done
