#! /usr/bin/awk -f

# Calculate density, defined as #ssRxns(y)/#ssSubs(y)

BEGIN{
	OFS="\t"
	FS="\t"
}

# Process first file: get all values of #ssRxns(y)
NR==FNR{
ssrx[$1]=$4
next
}

# Process second file: get #ssSubs(y)
{
sssub[$1]=$8
}

END{
for(y=1800;y<=2022;y++){
	print y , ssrx[y] / (sssub[y]+1e-6)
}

}

