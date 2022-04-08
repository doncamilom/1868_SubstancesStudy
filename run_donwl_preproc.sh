##! /usr/bin/bash
#
###############
## This code executes all steps for download of raw data from Reaxys, to produce of usable csv files.
###############
#
#date=$(date +%d_%b_%y)	# Store date as a variable
date="21_Feb_22"
#
#
## Download reactions from database using 20 processes (check periodically that 60M is large enough to contain all IDs)
#n=60000000
N=20
#errorStatus=$(./DATA/rxs_download.py -t R -n $n -N $N --name log --check_rxd) 
#while [ $errorStatus -eq -1 ]  # Loop until error status is positive :)
#do
#	sleep 300  # Sleep for 5 minutes in case it was a database reset error
#	errorStatus=$(./DATA/rxs_download.py -t R -n $n -N $N --name log --check_rxd -k) 
#done
#echo "Done downloading reactions"
#
#
## Download substances, using 20 processes (check 40M is large enough)
#n=40000000
#N=20
#errorStatus=$(./DATA/rxs_download.py -t S -n $n -N $N --name log) 
#while [ $errorStatus -eq -1 ]  # Loop until error status is positive :)
#do
#	sleep 300  # Sleep for 5 minutes in case it was a database reset error
#	errorStatus=$(./DATA/rxs_download.py -t S -n $n -N $N --name log -k) 
#done
#echo "Done downloading substances"


### Run preprocessing and routinary results on these data

# Extract rxn ID, reactants, reagents, etc..., dates. from rxn data
# Run getRxnDet.awk in parallel for all files in RXNs/p*
#mkdir "PREP/scr_$date"
#for i in $(seq 1 $N) # Iterate over directories (and parallel process each) so that argument list for parallel is short enough
#do
#	d=$(echo $i-1 | bc -l)	# compute $i - 1
#	parallel -j40 "PREP/getRxnDet.awk {} >> PREP/scr_$date/rxds{%}.tsv" ::: $(ls DATA/RXN/p$d/n*xml)
#done
#
### Produce some routinary plots (sanity check)
#
#diag_bin="DATA/DIAGNOSE/bin"		# Dir with code to produce diagnostic results
#diag_results="DATA/DIAGNOSE/Results"	# Dir to dump diagnostic results
#
## Concatenate results of previous script
#cat PREP/scr_$date/* > PREP/scr_$date/all_rxds.tsv
#
## Count number of single step reactions for every year
#./$diag_bin/numRxnsTime.awk PREP/scr_$date/all_rxds.tsv > $diag_results/countRxnsSS_$date.txt
## Count number of substances in single step reactions for every year (separated by role in reaction)
#./$diag_bin/numSubsTime.awk PREP/scr_$date/all_rxds.tsv > $diag_results/countSubsSS_$date.txt
## Compute density 
#./$diag_bin/calcDensity.awk $diag_results/countRxnsSS_$date.txt $diag_results/countSubsSS_$date.txt > $diag_results/density_$date.txt
## Count frequency of number of references per reaction and reaction details.
#./$diag_bin/numRefsFreq.awk PREP/scr_$date/all_rxds.tsv > $diag_results/countRefsFreqSS_$date.txt
## Count frequency of number of references and number of details per reaction.
#./$diag_bin/numVarsRefs.awk PREP/scr_$date/all_rxds.tsv > $diag_results/countVarsRefsSS_$date.txt
## Compute frequency of number of references per RXID (integrate over number of variations)
#awk 'BEGIN{FS="\t";OFS="\t"}{f[$1]+=$3}END{for(i in f) print i,f[i]}' $diag_results/countVarsRefsSS_$date.txt > $diag_results/numRefsRXID_$date.txt
## Compute frequency of number of variations per RXID (integrate over number of references)
#awk 'BEGIN{FS="\t";OFS="\t"}{f[$2]+=$3}END{for(i in f) print i,f[i]}' $diag_results/countVarsRefsSS_$date.txt > $diag_results/numVariationsRXID_$date.txt
#
## Make plots for the above calculations
#gnuplot -e "date='$date'" "$diag_bin/plotRxns.gpi" > _gnuplot_
#rm _gnuplot_
#
## Extract date from RXN for each substance in SS reactions (from file all_rxds.tsv)
#./PREP/getSubsDates.awk PREP/scr_$date/all_rxds.tsv > PREP/Data/subs_dates.tsv

##########
### Extracting data for PeriodSys project
##########

# Extract MFs from substance data, and merge with Dates found above
# This time, run parallel on the directories, as otherwise we'd have to load the dates file for every small file.
parallel "./PREP/getMFs.awk PREP/Data/subs_dates.tsv DATA/SUB/{}/* > PREP/scr_$date/id_mf_py{%}.tsv" ::: $(ls DATA/SUB/ -l | grep "^d" | awk '{print $NF}')

# Concatenate all these results
cat PREP/scr_$date/id_mf_py* > PREP/scr_$date/all_id_mf_py.tsv

# Further process this data to produce final usable file for project
./PREP/cleanMFs_final.awk PREP/scr_$date/all_id_mf_py.tsv > PREP/Data/MFs_ids_year.tsv

# Execute extraction of similarities for periodSys project
cd periodSys
./run_pipe.sh


# To run everything for a custom CS. e.g. CS of only organic substances, etc. Use pip_customCS.sh
# First argument is datafile (already formated with Perl script)
# Second argument is output dir. Here, a whole file system will be created to host results, scr, etc.
cd periodSys
./pipe_customCS.sh SplitOrgInorg/org_subs.tsv SplitOrgInorg/RunOrg &




##########
### Extracting data for reaction_template project
##########

#./PREP/getSetsRxn.awk PREP/scr_$date/all_rxds.tsv > PREP/Data/reaction_sets.tsv

