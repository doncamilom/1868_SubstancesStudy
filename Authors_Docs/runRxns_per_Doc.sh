#!/bin/bash

# Parallel process all raw XML files, and dump results into scr
for y in $(seq 0 19)
do
	parallel -j60 --max-args 500 "./rxns_per_doc.awk {} >> scr/rxns_per_doc_{%}.txt" ::: ../DATA/RXN/p$y/*
	echo "Done with data from p$y"
done

echo "Done with raw data. Now collecting results."


# Get number of authors publishing NEW REACTIONS each year
gawk '
BEGIN{OFS="\t"}

{
# Sum the results for this document, this year. This number is number of published new reactions in this document
rxns_total[$1][$2]+=$3	
}

END{	##### Print Number of authors that published this year
	for(y in rxns_total){
		for(cnr in rxns_total[y]){
			print y, cnr, rxns_total[y][cnr]
		}
	}
}' scr/rxns_per_doc* > rxns_per_doc.txt

