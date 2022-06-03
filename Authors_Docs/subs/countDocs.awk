#!/usr/bin/gawk -f

# Count number of docs involved in publishing new substances each year
# Also return mean and median number of new substances published per doc

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
all_subs[$2][$3]+=1	# Count new substs in doc $3
if($4>0)	docs[$2][$3]=$4		# Count auths in doc $3
}

END{
for(y in all_subs){
	# year, total docs, mean new substs per doc, median, mean num authors per doc, median
	print y, length(all_subs[y]), mean(all_subs[y]), median(all_subs[y]), mean(docs[y]), median(docs[y])
}

}
