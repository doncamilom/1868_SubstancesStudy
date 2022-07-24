#!/usr/bin/gawk -f

# Return, for each year,
# a distribution of the number of new reactions published in each document

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
	print "Year = " y
	if(length(docs[y])>0)
	for(d in docs[y]){
		print "",d, docs[y][d]
	}
}

}
