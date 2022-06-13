#!/bin/bash

# Parallel process all raw XML files, and dump results into scr
for y in $(seq 0 19)
do
	parallel -j60 --max-args 500 "./use_subs.awk {} >> scr/count_use{%}.txt" ::: ../../DATA/RXN/p$y/*
	echo "Done with data from p$y"
done

# Collect results

gawk '
BEGIN{
	FS="\t"
	OFS=FS
	fields="RXD.SXRN,RXD.RGTXRN,RXD.CATXRN,RXD.SOLXRN,RX.RXRN,RX.PRXN"	# Sources of substances
	split(fields,flds,",")

	# print header
	line="SUBID" FS "Year"
	for(i in flds){
		line=line FS flds[i]
	}
	print line
}

{
# 605461	2019	RXD.SXRN	12
# count[SID][yr][field]+=num
count[$1][$2][$3]+=$4
}

END{
for(s in count)
	for(yr in count[s]){
		line=s FS yr
		for(i in flds){
			field=flds[i]
			line=line FS count[s][yr][field]
		}
		print line
	}

}' scr/* > use_toolkit_subs.txt



