#! /usr/bin/env bash

STARTTIME=$(date +%s.%N)

# First run perl code that changes data format
perl format.pl $1 > ./Data/perlFormatted.tsv

# Then run awk code to find all possible substitution formulas
awk -f getAllRs.awk ./Data/perlFormatted.tsv > ./Data/allRs_dirty.tsv

# Finally, run awk code to find all non unique Rs + relationships between elements.
awk -f nonUniqRs.awk ./Data/allRs_dirty.tsv > ./Data/allRs_clean.tsv

ENDTIME=$(date +%s.%N)
echo "Elapsed time: $( echo "$ENDTIME - $STARTTIME" | bc -l ) seconds"
