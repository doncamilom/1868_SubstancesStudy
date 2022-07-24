#!/usr/bin/gawk -f

# Return, for each year,
# a distribution of the number of new reactions published by authors each year

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
}

{
split($5,auths,";")
for(i in auths){
	all_rxns[$1][auths[i]]+=1	# Count num rxns in which author participates this year
}
}

END{
for(y in all_rxns){
	print "Year = " y
	# year, num authors, mean new rxns per auth, median
	for(a in all_rxns[y]){
		print "", a, all_rxns[y][a]
	}
}

}
