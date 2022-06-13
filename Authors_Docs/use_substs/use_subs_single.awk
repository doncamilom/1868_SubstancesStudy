#! /usr/bin/awk -f

# Get the number of times some substance was used every year.

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
	SID="385737"

	OFS="\t"
	FS=OFS
	a="RXD.SXRN,RXD.RGTXRN,RXD.CATXRN,RXD.SOLXRN"	# Sources of substances
	split(a,rxdets,",")

	yq="CIT.PREPY,CIT.PPY,CIT.PY"	# Sources of years
	split(yq,year_q,",")
}

{
# for each RXD:
#	- check if SID is mentioned in some role.
#	- If so, count +1 to this role, this year


if(match($0,/<reaction index/)){
	rxd_num=0 # Reaction detail number. For creating unique ID for each RXN,RXD pair
	delete rctn
	delete prds
	# Delete here cause reactives and prods are not reported in each RXD, but only in main reaction.
	}

if(match($0,/<RX.RXRN>/)){ 
		if(match($0,SID)){
			rctn[gensub(/.*<RX.RXRN>(.+)<\/RX.RXRN>/, "\\1","g",$0)]=1
		}
}

if(match($0,/<RX.PXRN>/)){
		if(match($0,SID)){
			prds[gensub(/.*<RX.PXRN>(.+)<\/RX.PXRN>/, "\\1","g",$0)]=1
		}
}

# Delete every previous reaction detail, whenever <RXD> is found.
if(match($0,/<RXD>/)){
	delete subs
	delete years_rxd   
	yr=""
}

for(i in rxdets){	# Collect all other substances
	field=rxdets[i]
	if(match($0,field)){
		rgx=".*<" field ">(.+)</" field ">"
			if(match($0,SID)){	# Only append if SID matches line
				subs[field][gensub(rgx, "\\1","g",$0)]=1
			}
	}

}

for(i in year_q){  	# Collect years of pub for this RXD
	field=year_q[i]
	if(match($0,field)){
		rgx=".*<" field ">(.+)</" field ">"
		ind=length(years_rxd)+1
		years_rxd[ind]=gensub(rgx, "\\1","g",$0)  # Concatenate with previous entries found for this reaction detail
	}
}


if(match($0,/<\/RXD>/)){	# Last step: decide which substances are new and update data for them
	# Count +1 in every field where SID is in array
	yr=min(years_rxd)
		if(SID in rctn){
			count[yr]["RX.RXRN"]+=1
		}
		if(SID in prds){
			count[yr]["RX.PXRN"]+=1
		}
		for(field in subs){
			if(SID in subs[field]){
				count[yr][field]+=1
			}
		}
	}
}

END{
		for(yr in count){
			for(field in count[yr]){
				print yr, field, count[yr][field]
			}
		}
}

