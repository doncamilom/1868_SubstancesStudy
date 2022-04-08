#! /usr/bin/awk -f

# Extract relevant data from raw XML files for reactions, donwloaded from Reaxys,
# The script collects reactants, reagents, solvents, catalysts, and years of publication for every reaction:reaction detail pair


BEGIN{
	OFS="\t"
	a="RXD.SXRN,RXD.RGTXRN,RXD.CATXRN,RXD.SOLXRN,RXD.STP,CIT.PREPY,CIT.PPY,CIT.PY"
	split(a,rxdets,",")
	
	# Header can be printed and then placed in some file. Just for reference
	gsub(",","\t",a)
	print "ID\tRX.RXRN\t" a "\tRX.PXRN"  # Header
}

{

if(match($0,/<reaction index/)){
	rxd_num=0 # Reaction detail number. For creating unique ID for each RXN,RXD pair
	rctn=""
	prds=""
	}

if(match($0,/<RX.ID.+>/)){
	# Extract RX ID
	rxid=gensub(/.*<RX.ID.+>(.+)<\/RX.ID>/, "\\1","g",$0)
	}

if(match($0,/<RX.RXRN>/)){ # Concatenate all reactants
	rctn=rctn ":" gensub(/.*<RX.RXRN>(.+)<\/RX.RXRN>/, "\\1","g",$0)
}

if(match($0,/<RX.PXRN>/)){ # Concatenate all products
	prds=prds ":" gensub(/.*<RX.PXRN>(.+)<\/RX.PXRN>/, "\\1","g",$0)
}

# Delete every previous reaction detail, whenever <RXD> is found.
if(match($0,/<RXD>/)){
	delete data   
	rxd_num+=1
}

for(i in rxdets){  # Collect data for every field
	field=rxdets[i]
	if(match($0,field)){
		rgx=".*<" field ">(.+)</" field ">"
		data[field]=data[field] ":" gensub(rgx, "\\1","g",$0)  # Concatenate with previous entries found for this reaction detail
	}
}

	
# Dump all the collected data into a single line
if(match($0,/<\/RXD>/)){
	line=rxid ":" rxd_num "\t" rctn

	# Add reagents, catalysts, etc...
	for(i in rxdets){ 
		field=rxdets[i]
		sub(":","",data[field]) # Remove first ":" 
		line=line "\t" data[field]
	}
	print line, prds

}
# Done :)
}

