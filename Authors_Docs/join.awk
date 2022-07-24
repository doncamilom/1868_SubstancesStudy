#!/bin/gawk -f

# Use as:
# ./join.awk subs/detailedDistribNsubs.txt rxns/detailedDistribRxns.txt > detailedDistribAll.txt

BEGIN{
	FS="\t"
	OFS=FS
	#print header
	print "Year", "CNR", "Author", "# new substs", "# new reactions"
}

FNR==NR{	# Read substances file

# Store all fields in array
#	4000	1599111	Wizinger;Renckhoff	1
data[$1][$2]["subs"]=$4
data[$1][$2]["auts"]=$3
next
}

{	# Now store data for rxns
data[$1][$2]["rxns"]=$4
data[$1][$2]["auts"]=$3
}

END{

for(yr in data){
	for(cnr in data[yr]){
		n=split(data[yr][cnr]["auts"], auts, ";")

		if(n>0){
			for(i in auts){
				au=auts[i]
				print yr, cnr, au, data[yr][cnr]["subs"],data[yr][cnr]["rxns"]
			}
		}
		else	
				print yr, cnr, data[yr][cnr]["auts"], data[yr][cnr]["subs"],data[yr][cnr]["rxns"]
	}
}
}
