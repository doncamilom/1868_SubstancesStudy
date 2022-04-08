#!/usr/bin/gawk -f

# Produce distributions of stoichiometric subindex for each element

BEGIN{
FS="\t"
OFS=FS
}

{
split($2,s,"-")
for(e in s){
	split(s[e],en,":")
	d[en[1]][en[2]]+=1	# Distribution array
	d[en[1]]["total"]+=1
}
}

END{
for(e in d){
	print "============ Distribution for " e
	for(n in d[e]){
		if(n!="total")	print e,n,d[e][n]/d[e]["total"]
	}
}
}
