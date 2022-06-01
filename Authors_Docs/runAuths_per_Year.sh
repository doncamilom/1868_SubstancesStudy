#!/bin/bash

# Parallel process all raw XML files, and dump results into scr
#for y in $(seq 0 19)
#do
#	parallel -j60 --max-args 500 "./auth_new_rxns.awk {} >> scr/yr_num_auth_{%}.txt" ::: ../DATA/RXN/p$y/*
#	echo "Done with data from p$y"
#done
#
#echo "Done with raw data. Now collecting results."
#
## Collect all results in scr, with one processor.
#./collect_auth_new_rxns.awk scr/yr_num_auth* > yr_num_auth.txt

# Get number of authors publishing NEW REACTIONS each year
gawk '
BEGIN{OFS="\t"}

{
# Sum the results for this author, this year. This number is their number of published new reactions in year $1
auts_total[$1][$3]+=$2	
}

END{	##### Print Number of authors that published this year
	for(y in auts_total){
		print y, length(auts_total[y])
	}
}' scr/yr_num_auth* > yr_num_auth.txt

