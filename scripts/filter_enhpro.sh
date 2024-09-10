#!/bin/bash
#
# =Filter Tool=
# Made to retrieve enhancer and promoter fields from tables
# 
# Example of its use is:
# bash scripts/filter_enhpro.sh 10 RefSeq_FuncElems_chr1.tsv | column -s";" -t


nfield="$1"
awk -F "\t" '
{
	if ($'$nfield' ~ /enhancer/ || $'$nfield' ~ /promoter/)
	{
		if (!$11)
		{
			print $1,";",$10,";","no-note",";",$12;
		}
		else
		{
			print $1,";", $10,";","has-note",";",$12;
		}
	}
}
' "$2"
