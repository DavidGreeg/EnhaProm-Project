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
etype="$1";

awk 'BEGIN {FS="\t"}
{
	if ($'$etype' ~ /nhancer/ || $'$etype' ~ /romoter/)
	{
		print $1,$2,$3,$'$etype';
	}
}
' "$2"
