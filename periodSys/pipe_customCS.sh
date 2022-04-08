#!/usr/bin/bash

#####################
# Execute the whole pipeline, from raw data formating to production of optimal permutations and similarity classes.
#####################

inp_file=$1
out_dir=$2
home=$(pwd)
mkdir $out_dir/scr
mkdir $out_dir/Data
mkdir $out_dir/Results
mkdir $out_dir/log
cp $home/Data/ElementList.txt $out_dir/Data/
#
#### Decompose MFs into Rs
#./similarity/decomposeInRs.awk $inp_file > $out_dir/scr/all_decomp.tsv
#cd $out_dir/scr
#csplit all_decomp.tsv /^==*/ {*} -z --suppress-matched
#cd $home
#
## Get similarity relationships. Compute in parallel for each of the produced files above.
## sed is for cleaning trailing colon
#parallel -j30 "./similarity/getSimRel.awk {} | sed 's/\t:/\t/g' > $out_dir/scr/grouped_rs{%}" ::: $out_dir/scr/xx*
#cat $out_dir/scr/grouped_rs* > $out_dir/Data/allRs.tsv

## Convert all this data into python usable format
#./similarity/dataToPy.py -i $out_dir/

# Produce similarity matrices for years between 1800 and 2022, (takes ~ 4h using 50 processors)
#./similarity/simMat.py -n 50 -i $out_dir/ > $out_dir/log/sim.log

# Run optimization of ordering of elements with genetic algorithms
#rm $out_dir/log/genetic1D.log -f
#./Genetic1D/genetic1D.py -N 50 -i $out_dir/	# Run 50 optimizations each year, so we can pick the best 20.

# Explore historical evolution of PS. run comparison of the optimized permutations over time
# Produces the file ./Results/mathistory.npy
./Genetic1D/comparePSs.py -N 50 -i $out_dir/




