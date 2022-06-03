#!/usr/bin/gawk -f

# Count number of authors involved in publishing new reactions each year
# Also return mean and median number of new reactions each author is associated with, each year

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
	# year, num authors, mean new rxns per auth, median
	print y, length(all_rxns[y]), mean(all_rxns[y]), median(all_rxns[y])
}

}
