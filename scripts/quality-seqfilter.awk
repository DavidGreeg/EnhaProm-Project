#!/bin/awk
#
# =Phred Quality Sequence Filter Tool=
# 
# Example use:
# (suppose that you are located in:
#	~/Epigenomica/UBMI-IFC/EnhaProm-Project/scripts)
# awk -f scripts/quality-seqfilter.awk ~/XX01tiles_DEE1.fastq

function char_to_val(char){
	phred_chars="!\"#$%&\x27()*+,-./0123456789:;<=>?@ABCDEFGHIJK"
	phred_val=index(phred_chars, char)-1
	return phred_val
}
{
	if ( NR == 400001 ) #1001 ) #13 )
	{
		exit "terminated"
	}
	if( $1 ~ /^@HS1/ ) #( $1 == "@HS1" )
	{
		header=$0;
		header_line=NR; #-------------------------------
		getline sequence;
	}
	else if ( $0 ~ /^\+$/ && header_line == ((NR-2)) ) #( $0 == "+" && header_line == ((NR-2)) )
	{
		plus=$0;
		getline quality;

		q_sum=0;
		seq_len=length(quality);
		split(quality, q_chars, "");
		for (i=1; i<=seq_len; i++)
		{
			q_val=char_to_val(q_chars[i]);
			q_sum += q_val;
		}
		q_mean= q_sum/seq_len;

		print "=Header line: "header_line" =Header: ",header;
		print (q_mean >= 33)? "Line "NR": "quality"----STATUS: ACCEPTABLE" : "Line "NR": "quality"----STATUS: NOT ACCEPTABLE";
		print "\t-Mean: "q_mean"\n" #"\t-StDev: "line_sd;
	}
	else
	{
		next;
	}
}
