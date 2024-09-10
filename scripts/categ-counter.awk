#!/bin/awk
#
# =Category Counter Tool=
#
# Example use:
# (supposing you are located in:
#	~/Epigenomica/UBMI-IFC/EnhaProm-Project/)
# awk -f scripts/categ-counter.awk datasets/GeneHancer_RegElems_chr1.tsv

{ 
	elemType[$11]++ 
} 

END{ 
	for(type in elemType){ 
		print "type: ", type, "\tcounted:", elemType[type],"times" 
	} 
}
