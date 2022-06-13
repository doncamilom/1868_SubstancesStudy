#!/bin/gawk -f

BEGIN{
	OFS="\t"
	FS=OFS
}

{
# 2326046	1977	76415	Schaumann,E.; Grabley,F.-F.


auts=gensub("et al.",";et al.","g",$4)	# Treat 'et al.' as a different author. ACHTUNG
auts=gensub(" ","","g",auts)		# Remove spaces

data[$2][$3 FS auts]+=1
}

END{
for(yr in data){
	for(cnr in data[yr]){
			print yr, cnr, data[yr][cnr]
	}
}

}
