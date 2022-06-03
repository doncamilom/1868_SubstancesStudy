#!/usr/bin/gawk -f

# Count number of authors involved in publishing new substances each year
# Also return mean and median number of new substances each author is associated with, each year

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
	all_subs[$2][auths[i]]+=1	# Count num subs in which author participates this year
	
}
}

END{
for(y in all_subs){
	# year, num authors, mean new substs per auth, median
	print y, length(all_subs[y]), mean(all_subs[y]), median(all_subs[y])
}

}
