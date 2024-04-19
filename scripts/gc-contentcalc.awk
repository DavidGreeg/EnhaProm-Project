#!/bin/awk
#
# =GC Content Calculator Tool=

function base_perc(b_count,s_length){
	percent= b_count/seq_len
	return percent
}
BEGIN{
	split("A T C G",bases," ")

	# for (i=1; i<=4; i++){
	# 	print "base: ", bases[i]
	# }
}
{
	if ( NR == 400001 ) #1001 ) #13 )
	{
		exit "terminated"
	}
	if( $1 ~ /^@HS1/ ) #( $1 == "@HS1" )
	{
		getline sequence


		seq_len = length(sequence)
		split(sequence, seq_bases, "")
		
		for (i=1; i<=seq_len; i++){
			base_counts[seq_bases[i]]++
		}

		print "=Line:", NR+1, "\tSequence:", sequence, "\nPercentage per Base: "
		for (i=1; i<=4; i++){
			b_perc[bases[i]]= base_perc(base_counts[bases[i]],seq_len)
			printf("%s%%: %g \t", bases[i], b_perc[bases[i]])
			base_counts[bases[i]]=0
		}
		print "\tGC Percentage", b_perc["C"]+b_perc["G"]

	}
	else
	{
		next;
	}
}
