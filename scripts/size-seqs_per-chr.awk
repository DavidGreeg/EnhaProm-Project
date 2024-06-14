#!/bin/awk
#
# =Overlapping Sequences Filter Tool=
#
# Example use:
# (supposing you are located in:
#	~/Epigenomica/UBMI-IFC/EnhaProm-Project/)
# tail -n500 datasets/FEbiolregions_RS_2023_10_T2T-CHM13v2.0.bed | awk -f scripts/overlap-filter.awk

{
	post_start = $2;
	post_end = $3;

	if ($1 == "chr1")
	{
		print post_end - post_start
	}
	# if ( NR == 1 )
	# {
	# 	prev_start = post_start;
	# 	prev_end = post_end;
	# }
	# else if ( post_start < prev_start )
	# {
	# 	print "posterior element starts at:", post_start, " which before previous element which starts at:", prev_start
	# }
	# else if ( post_start == prev_start || post_end == prev_end)
	# {
	# 	#coincidental ends
	# 	print "subset", NR, $1, prev_start, post_end
	# }
	# else if ( post_start < prev_end )
	# {
	# 	#overlap 
	# 	print "overlap", NR, $1, prev_start, post_end
	# }
	#
	# prev_start = post_start;
	# prev_end = post_end;
}
