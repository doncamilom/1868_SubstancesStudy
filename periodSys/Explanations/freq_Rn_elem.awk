#!/usr/bin/gawk -f

BEGIN{OFS="\t";FS=OFS}


{

split($1,rcomp,/:[0-9]+-/)	# Composition of R
n=substr(rcomp[length(rcomp)],3)

split($2, elem, /:/)
split($4, year, /:/)

for(e in elem){
	for(i in rcomp){
		all_el[elem[e]][year[e]][rcomp[i]]+=1	# element:year:element_in_R --> freq
		#print elem[e], year[e], rcomp[i]
	}
}

}

END{

for(e in all_el){	# For each element
	print "Element X = " e

	for(y in all_el[e]){	# For every year
		print "\tYear " y

		for(re in all_el[e][y]){
			if(!match(re,/X:[0-9]+/))
				print re, ":", all_el[e][y][re]
		}
	}
}

}
