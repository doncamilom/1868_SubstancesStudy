#!/usr/bin/gawk -f

function mean(arr){
	sum=0
	for(i in arr){
		sum+=arr[i]
	}
	return sum/length(arr)

}

function median(arr){
	n=asort(arr,sarr)
	med=sarr[int(length(arr)/2)]
	return med
}

BEGIN{
	OFS="\t"
	FS=OFS
	rxns_year[4000][1]=0
}

{
# This array contains, for each year, an array of the new reactions in a publication of that year
rxns_year[$1][1]=0
i=length(rxns_year[$1])+1
rxns_year[$1][i]=$3
}

END{
for(y in rxns_year){
	cmean=mean(rxns_year[y])
	cmedian=median(rxns_year[y])
	cnt=length(rxns_year[y])

	# For every year, print:
	# year, mean # rxns per doc, median # rxns per doc, total docs.
	print y, cmean, cmedian, cnt
	
}

}
