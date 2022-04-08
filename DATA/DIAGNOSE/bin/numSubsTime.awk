#! /bin/awk -f

### This script finds the date of publication for every substance in single step reactions.
# Fields:
# ID	RX.RXRN	RXD.SXRN	RXD.RGTXRN	RXD.CATXRN	RXD.SOLXRN	RXD.STP	CIT.PREPY	CIT.PPY	CIT.PY	RX.PXRN

# TODO: Produce counts for substances first reported as reactants, reagents, ..., products.

function min(arr){
	if(arr[1]!=""){
		m=arr[1]
		for(i in arr){
			if(arr[i]<m) m=arr[i]
		}
		return m
	}
	else return 4000
}

BEGIN{
	FS="\t"; 
	OFS="\t"
	split("2,3,4,5,6,11",is,",")
}

{
if($7==1){
# Obtain dates. Work only with academic year (py) as ppy has problems.
split($10,py,":")
mpy=min(py)

## Extract all substances in line
# Substances are in $2, $3, $4, $5, $6, $11 (products are in $11)

for(i in is){ # Columns where substances are
	delete subs
	split($is[i],subs,":") # Create arrays of substances for each of $2, $3, ...
	for(s in subs){	# each substance in reaction role $is[i]
		if(y[subs[s] ":" i]){	# If an entry exists for this substance in this role
			if(mpy<y[subs[s] ":" i])	y[subs[s] ":" i]=mpy # Update only if date is older
			# If already exists but date is larger, do nothing
		}
		else{
			y[subs[s] ":" i]=mpy
		}
		g[subs[s]] # Add subs to list of substances
	}
}
}
}



END{
# Loop through all substances checking what they've been reported as, and in what year

for(s in g){
	for(i=1;i<=6;i++)	arr[i]=y[s ":" i]

	# Find which of these is minimum
	m=4000
	am=1
	for(i in arr){
		if(arr[i] && arr[i]<m){	m=arr[i];	am=i	}
	}

	# Increase this year's count, for this role.
	c[m ":" am]+=1

	yrs[m]
}

print "#####\n"
for(yr in yrs){
	line=yr
	sum=0
	for(i=1;i<=6;i++){
		line=line "\t" c[yr ":" i]
		sum+=c[yr ":" i]
	}
	print line, sum
	# Fields:
	#	year	RX.RXRN	RXD.SXRN	RXD.RGTXRN	RXD.CATXRN	RXD.SOLXRN	RX.PXRN	Sum
}
}
	
