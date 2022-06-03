#! /usr/bin/awk -f

# Get all NEW SUBSTANCES published, with year of publication and authors associated.
# TODO maybe later: count mean and median number of new substs published per document.

# Extract data from raw XML files for reactions, donwloaded from Reaxys,
# The script collects reactants, reagents, solvents, catalysts.

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

function update(subs, new_yr, auts_rxd, cnr_rxd){		# Update this substance's data.


	
	if(subs in all_subs){
		if(new_yr<all_subs[subs]["year"]){	# Update all data for this substance
			all_subs[subs]["year"]=new_yr
			all_subs[subs]["auth"]=auts_rxd
			all_subs[subs]["cnr"]=cnr_rxd
		}
	}
	else{
		all_subs[subs]["year"]=new_yr
		all_subs[subs]["auth"]=auts_rxd
		all_subs[subs]["cnr"]=cnr_rxd
	}
	return 0

}

BEGIN{
	OFS="\t"
	FS=OFS
	a="RXD.SXRN,RXD.RGTXRN,RXD.CATXRN,RXD.SOLXRN"	# Sources of substances
	split(a,rxdets,",")

	yq="CIT.PREPY,CIT.PPY,CIT.PY"	# Sources of years
	split(yq,year_q,",")
}

{
# for each RXD: 
# 	- get all substances in an array
#	- get authors
# 	- get date of pub
#	- get CNR (Id for doc)
#	- When encounter </RXD>, check every substance and update a big array of substances:
#		- big_arr[subst_i][year]
#		- big_arr[subst_i][auths]
#		(obvs only update if old year is higher than year of this RXD)


if(match($0,/<reaction index/)){
	rxd_num=0 # Reaction detail number. For creating unique ID for each RXN,RXD pair
	delete subs_rxd		# Delete here cause reactives and prods are not reported in each RXD, but only in main reaction.
	}

if(match($0,/<RX.RXRN>/)){ # Add all reactants to subs_rxd array
	ind=length(subs_rxd)+1
	subs_rxd[ind]=gensub(/.*<RX.RXRN>(.+)<\/RX.RXRN>/, "\\1","g",$0)
}

if(match($0,/<RX.PXRN>/)){ # Add all products to subs_rxd array
	ind=length(subs_rxd)+1
	subs_rxd[ind]=gensub(/.*<RX.PXRN>(.+)<\/RX.PXRN>/, "\\1","g",$0)
}

# Delete every previous reaction detail, whenever <RXD> is found.
if(match($0,/<RXD>/)){
	delete subs_rxd_other
	delete years_rxd   
	new_yr=""
	auts_rxd=""
	cnr_rxd=""
}

for(i in rxdets){	# Collect all other substances
	field=rxdets[i]
	if(match($0,field)){
		rgx=".*<" field ">(.+)</" field ">"
		ind=length(subs_rxd_other)+1
		subs_rxd_other[ind]=gensub(rgx, "\\1","g",$0)  # Concatenate with previous entries found for this reaction detail
		#print field, subs_rxd_other[ind]
		
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

# Collect author list
if(match($0,/<CIT.AU>/)){
	auts_rxd=gensub(/.*<CIT.AU>(.+)<\/CIT.AU>/, "\\1","g",$0)
}

# Collect cnr
if(match($0,/<CNR.CNR>/)){
	# New citation for this RXD
	cnr_rxd=gensub(/.*<CNR.CNR>(.+)<\/CNR.CNR>/, "\\1","g",$0)
}
	
if(match($0,/<\/RXD>/)){	# Last step: decide which substances are new and update data for them


# Do for reactants and products
for(i in subs_rxd){
	subs=subs_rxd[i]
	new_yr=min(years_rxd)
	update(subs, new_yr, auts_rxd, cnr_rxd)
}

# Do for all the rest
for(i in subs_rxd_other){
	subs=subs_rxd_other[i]
	new_yr=min(years_rxd)
	update(subs, new_yr, auts_rxd, cnr_rxd)
}


}
}

END{
for(subs in all_subs){
	# For each substance, print SUBS_ID, pub_year, cit_pub_year, author_list
	print subs, all_subs[subs]["year"], all_subs[subs]["cnr"], all_subs[subs]["auth"]
}

}

