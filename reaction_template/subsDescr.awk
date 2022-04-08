#! /usr/bin/awk -f

# This script takes all reaction templates produced by getRxnTemps.awk and, for each substance, constructs the associated set of reaction templates.

# Takes as arguments: a grouping element (4 digit code), and the file with all reaction templates.
# Aim is to read in parallel for all possible grouping elements.

BEGIN{
	FS="\t";
	OFS="\t"}

# First, read grouping element (4 digit number given as first argument)
NR==FNR{ge=$0; next}

# Using this grouping element, filter to files with such element only

$1~ge{
# Collect and concatenate information
if(count[$2]<50000){	# Keep increasing string only if ammount of data is reasonable
	templ[$2]=templ[$2] "-" $3
	ids[$2]=ids[$2] "-" $4
	dates[$2]=dates[$2] "-" $5
	count[$2]+=1
}
else if(overcount[$2]==0){
	# Activate flag and delete all information for this substance
	overcount[$2]=1
	templ[$2]=0
	ids[$2]=0
	dates[$2]=0
	print "Ignoring substance with ID" , $2
}
}

END{
for(s in templ){
	print s, count[s], templ[s], ids[s], dates[s]
}
}
