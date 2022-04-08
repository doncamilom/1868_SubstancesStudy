#! /bin/awk -f

### Canonicalize all reactant sets (reactants, reagents, etc...): Write ordered as a single field
# Fields:
# ID	RX.RXRN	RXD.SXRN	RXD.RGTXRN	RXD.CATXRN	RXD.SOLXRN	RXD.STP	CIT.PREPY	CIT.PPY	CIT.PY	RX.PXRN

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
	split("2,3,4,5,6",is,",")
}

{
split($9,ppy,":")
split($10,py,":")
mpy=min(py)
mppy=min(ppy)
if(mpy<mppy)	date=mpy
else		date=mppy

delete set

## Extract all substances in line
# Relevant substances are in $2, $3, $4, $5, $6 (don't use products).
if(date!=4000){
for(i in is){
	split($is[i],subs,":") # Create arrays inside array subs
	for(s in subs){
		if(subs[s]){
			c=subs[s]
			if(!(c in set)){	# Only add new entry if this compound is not yet in array
				set[c]
			}
		}
	}
}
}

# Add reaction set to dataset
if($7==1){	# Only if reaction is single step


asorti(set)
line=""
for(i in set)	line=line ":" set[i]

if(line in rxs){
	if(date<rxs[line]){
		rxs[line]=date	# Update date for this reaction only if new date lower
		id[line]=$1
		cs[line]=length(set)
	}
}
else{
	rxs[line]=date	# Update date for this reaction only if new date lower
	id[line]=$1
	cs[line]=length(set)
}
}
}

END{
for(r in rxs){
	print id[r], cs[r], r, rxs[r]	# Print canonicalized rxn set and date
}

}
