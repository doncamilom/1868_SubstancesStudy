#! /bin/awk -f

# Extract number of single step reactions, both in academic journals and patents yearly from output of getRxnDet.awk.
# Fields:
# ID	RX.RXRN	RXD.SXRN	RXD.RGTXRN	RXD.CATXRN	RXD.SOLXRN	RXD.STP	CIT.PREPY	CIT.PPY	CIT.PY

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
	mppy=4000
	mpy=4000
}

{
split($1,id_rxd,":")
if(FNR==1)	id=id_rxd[1]

# We have RXID:RXD. Count only reactions on one side, and reactions+rxds in parallel
if($7==1){  # Count only single step
	# Get years
	split($9,ppy,":")
	split($10,py,":")
		
	if(id_rxd[1]!=id){	# New reaction scheme. 
		# Write years for previous RXID
		if(mppy<mpy && mppy!=4000)	ypp[mppy]+=1
		else if(mpy!=4000)		yp[mpy]+=1
		
		# Extract year from this rxid:rxd
		mppy=min(ppy)
		mpy=min(py)
		id=id_rxd[1]
	}
	else{	# If this is an additional rxd to the same reaction
		# Compare this entrie's years against previously collected
		currppy=min(ppy)
		currpy=min(py)
		if(currppy<mppy)	mppy=currppy
		if(currpy<mpy)		mpy=currpy
	}

	# Write rxn+rxd dates
	ppyrxd[min(ppy)]+=1
	pyrxd[min(py)]+=1
}
}

END{
	for(y=1771;y<=2022;y++){
		print y, ypp[y], yp[y], ypp[y]+yp[y], ppyrxd[y], pyrxd[y], ppyrxd[y]+pyrxd[y]
	}
}
	
