#! /usr/bin/gawk -f

# Get number of NEW REACTIONS published, per document

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

if(match($0,/<CNR.CNR>/)){
	# New citation for this RXD
	cit_rxd=gensub(/.*<CNR.CNR>(.+)<\/CNR.CNR>/, "\\1","g",$0)
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
	cits_total[rxn_year][cit]+=1	# Add 1 to the number of reactions published in this doc.
}

}

END{

#### This prints year, num. of contributions of this author this year, author
for(y in cits_total){
	for(i in cits_total[y]){	# Dump Year, CNR, Num New Reactions Published in this doc.
		print y, i, cits_total[y][i]
	}

}
}


