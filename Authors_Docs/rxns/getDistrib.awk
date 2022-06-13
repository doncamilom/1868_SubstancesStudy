#!/bin/gawk -f

BEGIN{
	OFS="\t"
	FS=OFS
}

{
#1923	14003	1223922	2	Robinson,G.;Robinson,R.	1

data[$1][$3 FS $5]+=1
}

END{
for(yr in data){
	for(cnr in data[yr]){
			print yr, cnr, data[yr][cnr]
	}
}

}
