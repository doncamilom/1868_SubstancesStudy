#!/bin/bash

# Parallel process all raw XML files, and dump results into scr
for y in $(seq 0 19)
do
	parallel -j60 --max-args 500 "./new_subs_docs_auth.awk {} >> scr/subs_yr_cnr_auth_{%}.txt" ::: ../../DATA/RXN/p$y/*
	echo "Done with data from p$y"
done

echo "Done with raw data. Now collecting results."

# Just do the same process once again: Get substances + substance data, and update every time.
gawk '
BEGIN{OFS="\t"}
{
# For each substance, add to array, and update array data depending on year.

if($1 in all_subs){
	if($2<all_subs[$1]["y"]){	# Update only if year is less than previous
		all_subs[$1]["y"]=$2
		all_subs[$1]["cnr"]=$3
		all_subs[$1]["aut"]=$4
	}
}
else{				# Or if substance not yet in list
	all_subs[$1]["y"]=$2
	all_subs[$1]["cnr"]=$3
	all_subs[$1]["aut"]=$4
}

}

END{	##### Print Number of authors that published this year
	for(subs in all_subs){
		auth=all_subs[subs]["aut"]
		auth=gensub("et al.",";et al.","g",auth)	# Treat "et al." as a different author. ACHTUNG
		auth=gensub(" ","","g",auth)		# Remove spaces

		# Check: if last character of auth is ;, then n-=1
		lt=substr(auth , length(auth))

		if(lt==";")
			auth=substr(auth, 1, length(auth)-1)

		n=split(auth,aut_l,";")

		# SUB_ID, YEAR, CNR, NUM AUTHORS, AUTHOR LIST
		print subs, all_subs[subs]["y"],all_subs[subs]["cnr"], n, auth
	}
}' scr/subs_yr_cnr_auth_* > subs_yr_cnr_auth.txt

