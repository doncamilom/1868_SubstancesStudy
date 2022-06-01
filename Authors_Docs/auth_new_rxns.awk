#! /usr/bin/gawk -f

# Get number of authors publishing NEW REACTIONS each year

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
	OFS="\t"
}

{
if(match($0,/<reaction index/)){
	rxn_year=4000
	}

if(match($0,/<RX.ID high/)){
	rxid=gensub(/.*<RX.ID highlight="true">(.+)<\/RX.ID>/, "\\1","g",$0)
}

if(match($0,/<RXD>/)){
	# New RXD, meaning collect new set of citations
	yi=1
	delete year
}

if(match($0,/<citation index/)){
	# New citation for this RXD
	cit_rxd=gensub(/.*<citation index="(.+)">/, "\\1","g",$0)
}

if(match($0,/<CIT.AU>/)){
	auts_rxd=gensub(/.*<CIT.AU>(.+)<\/CIT.AU>/, "\\1","g",$0)
	#print auts_rxd
}

########### Collect year for this RXD
if(match($0,/<CIT.PREPY>/)){	
	year[yi]=gensub(/.*<CIT.PREPY>(.+)<\/CIT.PREPY>/, "\\1","g",$0)
	yi+=1
}
if(match($0,/<CIT.PY>/)){
	year[yi]=gensub(/.*<CIT.PY>(.+)<\/CIT.PY>/, "\\1","g",$0)
	yi+=1
}
if(match($0,/<CIT.PPY>/)){
	year[yi]=gensub(/.*<CIT.PPY>(.+)<\/CIT.PPY>/, "\\1","g",$0)
	yi+=1
}
########### Finish collecting year


if(match($0,/<\/RXD>/)){	# If done with this RXD
	# Decide whether this RXD is the first for this reaction (oldest year, as we only want new reactions now)

	rxd_year=min(year) # is defined as minimum among prepy, py and ppy

	if(rxn_year>=rxd_year){
		rxn_year=rxd_year
		cit=cit_rxd
		auts=auts_rxd
	}
}

if(match($0,/<\/reaction>/)){
	# Get list of authors
	auts=gensub("et al.",";et al.","g",auts)	# Treat 'et al.' as a different author. ACHTUNG
	auts=gensub(" ","","g",auts)		# Remove spaces
	n=split(auts,aut_l,";")

	for(i in aut_l){
		auts_total[rxn_year][aut_l[i]]+=1	# Add 1 to the number of works published by this author, this year.
	}

	#print rxid, cit, rxn_year, n, auts
}

}

END{

#	##### Print Number of authors that published this year
#	for(y in auts_total){
#		print y, length(auts_total[y])
#	}




#### This prints year, num. of contributions of this author this year, author
for(y in auts_total){

	for(i in auts_total[y]){
		print y,auts_total[y][i], i
	}

}

}



# Delete every previous reaction detail, whenever <RXD> is found.
#if(match($0,/<RXD>/)){
#	delete data   
#	rxd_num+=1
#}
#
#for(i in rxdets){  # Collect data for every field
#	field=rxdets[i]
#	if(match($0,field)){
#		rgx=".*<" field ">(.+)</" field ">"
#		data[field]=data[field] ":" gensub(rgx, "\\1","g",$0)  # Concatenate with previous entries found for this reaction detail
#	}
#}
#
#	
## Dump all the collected data into a single line
#if(match($0,/<\/RXD>/)){
#	line=rxid ":" rxd_num "\t" rctn
#
#	# Add reagents, catalysts, etc...
#	for(i in rxdets){ 
#		field=rxdets[i]
#		sub(":","",data[field]) # Remove first ":" 
#		line=line "\t" data[field]
#	}
#	print line, prds
#
#}


