#!/usr/bin/bash

#####################
# Execute the whole pipeline, from raw data formating to production of optimal permutations and similarity classes.
#####################

### Copy datafile to this dir
#cp ../PREP/Data/MFs_ids_year.tsv ./Data/.
##
### Format data into tractable format for awk
#./similarity/format.pl ./Data/MFs_ids_year.tsv > Data/format_MFs.tsv
##
### Decompose MFs into Rs
#./similarity/decomposeInRs.awk Data/format_MFs.tsv > scr/all_decomp.tsv
#cd scr
#csplit all_decomp.tsv /^==*/ {*} -z --suppress-matched
#cd ..
#
## Get similarity relationships. Compute in parallel for each of the produced files above.
## sed is for cleaning trailing colon
#parallel -j30 "./similarity/getSimRel.awk {} | sed 's/\t:/\t/g' > ./scr/grouped_rs{%}" ::: scr/xx*
#cat scr/grouped_rs* > Data/allRs.tsv

# Convert all this data into python usable format
./similarity/dataToPy.py 

# Produce similarity matrices for years between 1800 and 2022, (takes ~ 4h using 50 processors)
./similarity/simMat.py -n 50 > ./log/sim.log

# Run optimization of ordering of elements with genetic algorithms
rm ./log/genetic1D.log -f
./Genetic1D/genetic1D.py -N 50	# Run 50 optimizations each year, so we can pick the best 20.

# Explore historical evolution of PS. run comparison of the optimized permutations over time
# Produces the file ./Results/mathistory.npy
./Genetic1D/comparePS.py -N 50

# For this step, we have to switch to k70 as the module we need for that is only here.
# Run code for finding families of elements using computer vision algorithms
./familiesElem/findGroups.py -N 20



