#! /bin/awk -f

# For each reaction id (reaction scheme) count frequency of number of references and number of variants
# Fields:
# ID	RX.RXRN	RXD.SXRN	RXD.RGTXRN	RXD.CATXRN	RXD.SOLXRN	RXD.STP	CIT.PREPY	CIT.PPY	CIT.PY	RX.PXRN

BEGIN{
	FS="\t"; 
	OFS="\t"
}

{
# Compute number of references as total number of entries in CIT.PPY and CIT.PY together.
split($1,id_rxd,":")
if(FNR==1)	id=id_rxd[1]

if($7==1){  # Count only single step
	n=split($9,_,":")+split($10,_,":")	# Total number of references for this rxn detail

	if(id_rxd[1]!=id){	# New reaction scheme. 

		if(idlist[id]){	# Check if this (previous) RXID had been treated before
			if(refs_id>idlist[id]){	# If current entry is higher
				freq[idlist[id] ":100"]-=1	# Reduce freq of previous entry by 1 (only happens with more than 100 vars, RXD problem)
				freq[refs_id ":" vars]+=1	# If entry existed before, write only if this one is larger
			}
		}
		else{	# If entry didn't exist, count it normally
			freq[refs_id ":" vars]+=1
			# And add it to list of known ids
			idlist[id]=refs_id
		}
		# Reset this rxid's count, starting from n
		refs_id=n
		# Update id
		id=id_rxd[1]
		# Restart count of num of variants
		vars=1
	}
	else{	# If this is an additional rxd to the same reaction
		# Compare this entrie's years against previously collected
		refs_id+=n	# Increment this ID's count
		vars+=1
	}
}
}

END{
for(i in freq){
	split(i,is,":")
	# Num refs, num vars, freq
	print is[1],is[2],freq[i]
}
}
	
