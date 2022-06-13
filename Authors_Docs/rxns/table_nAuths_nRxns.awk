#!/usr/bin/gawk -f

# Calculate table of frequencies in documents, one for each year:
# row: number of new reactions
# col: number of authors


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
	FS="\t"
	OFS=FS
}

{
all_rxns[$1][$3]+=1	# Count new rxns in doc $3
if($4>0)	docs[$1][$3]=$4		# Count auths in doc $3
}

END{
# First, count frequencies, that is merge `all_rxns` and `docs`
for(y in all_rxns){	# For each year
	for(d in all_rxns[y]){	# for each document
		# 3D array with indexes: year, #rxns, #auths. Values are freqs of combinations.
		T[y][all_rxns[y][d]][docs[y][d]] += 1		
	}
}

# Now plot these tables
max_na=100
for(y in T){
	print "Year = " y
	# print header
	head=""
	for(na=0; na<=max_na; na++){	
		head=head OFS na	
	}
	print head


	for(nr in T[y]){	# num new reactions
		line = nr
		for(na=0; na<=max_na; na++){	# num authors (suppose max value is 100)

			line=line OFS T[y][nr][na]
		}
		print line
	}
}



}
