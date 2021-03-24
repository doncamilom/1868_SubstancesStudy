#! /usr/bin/env bash

STARTTIME=$(date +%s.%N)

python3.7 getTabsI.py $1 $2 $3
source nonUniqRs.sh
python3.7 getTabsII.py $2

ENDTIME=$(date +%s.%N)
echo "Time elapsed: $( echo "$ENDTIME - $STARTTIME" | bc -l ) seconds"

