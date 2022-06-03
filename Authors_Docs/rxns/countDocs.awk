#!/usr/bin/gawk -f

# Count number of docs involved in publishing new reactions each year
# Also return mean and median number of new reactions published per doc

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
all_rxns[$1][$3]+=1	# Count new rxns in doc $3
if($4>0)	docs[$1][$3]=$4		# Count auths in doc $3
}

END{
for(y in all_rxns){
	# year, total docs, mean new rxns per doc, median, mean num authors per doc, median
	if(length(docs[y])>0)
		print y, length(all_rxns[y]), mean(all_rxns[y]), median(all_rxns[y]), mean(docs[y]), median(docs[y])
}

}
