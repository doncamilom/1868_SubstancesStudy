#! /bin/awk -f

### This program extracts the absolute minimum date associated with each substance, from reactions
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
	split("2,3,4,5,6,11",is,",")
}

# If we only care about PY, use this:
{
split($10,py,":")
mpy=min(py)

## Extract all substances in line
# Substances are in $2, $3, $4, $5, $6, $11
for(i in is){
	split($is[i],subs,":") # Create arrays inside array subs
	for(s in subs){
		if(mpy!=4000){
			if(all[subs[s]]){ # If an entry exists for this substance
				if(mpy<all[subs[s]])		all[subs[s]]=mpy # Update only if mpy is older
			}
			else all[subs[s]]=mpy
		}
	}
}

}

# Use this in case substances in patents are also considered
#{
## Extract minimum amongst all py and ppy.
## They're in $9, $10
#split($9,ppy,":")
#split($10,py,":")
#mppy=min(ppy)
#mpy=min(py)
#
#if(mppy==1800 || mpy==1800) print $9,$10; c+=1
#
#if(xor(mppy!=4000,mpy!=4000)){ # Add date only if at least one date exists
#	if(mppy<mpy) date=mppy
#	else date=mpy
#}
#else date=9999
#
### Extract all substances in line
## Substances are in $2, $3, $4, $5, $6, $11
#for(i in is){
#	split($is[i],subs,":") # Create arrays inside array subs
#	for(s in subs){
#		if(all[subs[s]]){ # If an entry exists for this substance
#			if(date<all[subs[s]])		all[subs[s]]=date # Update only if date is older
#		}
#		else all[subs[s]]=date
#	}
#}
#}

	
END{
for(s in all){
	print s, all[s]
}
}
	
