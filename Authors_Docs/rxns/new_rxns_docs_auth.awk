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

if(match($0,/<citation index/)){
	# New citation for this reaction
	yi=1
	delete year
	cit_rxd=""
	auts_rxd=""
	dt_rxd=""
}

if(match($0,/<CNR.CNR>/)){
	# New citation for this RXD
	cit_rxd=gensub(/.*<CNR.CNR>(.+)<\/CNR.CNR>/, "\\1","g",$0)
}

if(match($0,/<CIT.AU>/)){
	auts_rxd=gensub(/.*<CIT.AU>(.+)<\/CIT.AU>/, "\\1","g",$0)
}

if(match($0,/<CIT.DT>/)){
	dt_rxd=gensub(/.*<CIT.DT>(.+)<\/CIT.DT>/, "\\1","g",$0)
}

########### Collect year for this citation
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


if(match($0,/<\/citation>/)){	# If done with this RXD
	# Decide whether this citation is the first for this reaction (oldest year, as we only want new reactions now)

	rxd_year=min(year) # is defined as minimum among prepy, py and ppy

	if(rxn_year>=rxd_year){
		rxn_year=rxd_year
		cit=cit_rxd
		auts=auts_rxd
		dt=dt_rxd
	}
}

if(match($0,/<\/reaction>/)){
	# End of reaction statement

	# Get list of authors
	auts=gensub("et al.",";et al.","g",auts)	# Treat 'et al.' as a different author. ACHTUNG
	auts=gensub(" ","","g",auts)		# Remove spaces
	n=split(auts,aut_l,";")
	
	if(dt=="Article")	dtn=1
	else			dtn=0

	#year, rxid, cnr, num auts, auth list, doc type (article=1, patent=0)
	print rxn_year, rxid, cit, n, auts, dtn
}

}

