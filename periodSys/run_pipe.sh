#!/bin/bash

#####################
# Execute the whole pipeline, from raw data formating to production of optimal permutations and similarity classes.
#####################

inp_preproc=''
inp_file='./Data/MFs_ids_year.tsv'
out_dir='.'

get_rs=''
get_sim=''
opt_perm=''
comparePS=''
findGrps=''

print_usage() {
  printf "\n\nUsage %s:\n\n" $0 
  printf "\n	[-i] Input file.\n"
  printf "		Default is ./Data/MFs_ids_year.tsv\n"
  printf "	[-p] Control what the nature of input file is.\n"
  printf "		If -p: File should be already preprocessed using format.pl\n"
  printf "		Else: File is file with raw MFs, index and date.\n"
  printf "	[-o] Directory to output all files. Default is '.'\n"
  printf "		Choose an empty directory, the whole filesystem will be produced there.\n\n"
  printf "	[-R] Compute all possible sub-molecular formula/subindex pair in dataset\n"
  printf "	[-S] Compute Similarity Matrices for every year.\n"
  printf "	[-P] Optimize permutations of sequence of elements, given similarity matrices.\n"
  printf "	[-H] Compare the optimized permutations in time.\n"
  printf "	[-G] Compute families of elements. ***Need to change to k70*** \n"
}

# Parse arguments
while getopts 'i:o:pRSPHG' flag; do
  case "${flag}" in
    i) inp_file="$OPTARG" ;;
    p) inp_preproc=1 ;;
    o) out_dir="$OPTARG" ;;
    R) get_rs=1 ;;
    S) get_sim=1 ;;
    P) opt_perm=1 ;;
    H) comparePS=1 ;;
    G) findGrps=1 ;;
    ?) print_usage
       exit 1 ;;
  esac
done


home=$(pwd)
mkdir $out_dir/scr
mkdir $out_dir/Data
mkdir $out_dir/Results
mkdir $out_dir/log
cp $home/Data/ElementList.txt $out_dir/Data/


if [ -z "$inp_preproc" ]; then	# If data is reported to NOT be preprocessed with format.pl
	## Format data into tractable format for awk
	echo "Formating $inp_file and saving in $out_dir/Data/format_MFs.tsv"
	#./similarity/format.pl $inp_file > $out_dir/Data/format_MFs.tsv
	inp_file="$out_dir/Data/format_MFs.tsv"
fi

if [ -n "$get_rs" ]; then
	## Decompose MFs into Rs
	./similarity/decomposeInRs.awk $inp_file > $oit_dir/scr/all_decomp.tsv

	cd $out_dir/scr
	csplit all_decomp.tsv /^==*/ {*} -z --suppress-matched
	cd $home
fi

if [ -n "$git_sim" ]; then
	## Get similarity relationships. Compute in parallel for each of the produced files above.
	## sed is for cleaning trailing colon
	parallel -j30 "./similarity/getSimRel.awk {} | sed 's/\t:/\t/g' > $out_dir/scr/grouped_rs{%}" ::: $out_dir/scr/xx*
	cat $out_dir/scr/grouped_rs* > $out_dir/Data/allRs.tsv
	
	## Convert all this data into python usable format
	python3 similarity/dataToPy.py -i $out_dir/
	
	# Produce similarity matrices for years between 1800 and 2022, (takes ~ 4h using 50 processors)
	python3 similarity/simMat.py -n 50 -i $out_dir/ > $out_dir/log/sim.log
fi

if [ -n "$opt_perm" ]; then
	# Run optimization of ordering of elements with genetic algorithms
	rm $out_dir/log/genetic1D.log -f
	python3 Genetic1D/genetic1D.py -N 50 -i $out_dir/	# Run 50 optimizations each year, so we can pick the best 20.
fi

if [ -n "$comparePS" ]; then
	# Explore historical evolution of PS. run comparison of the optimized permutations over time
	# Produces the file ./Results/mathistory.npy
	python3 Genetic1D/comparePSs.py -i $out_dir/
fi

if [ -n "$findGrps" ]; then
	# For this step, we have to switch to k70 as the module we need for that is only here.
	# Run code for finding families of elements using computer vision algorithms
	python3 familiesElem/findGroups.py -N 20	#TODO modify script to accept $out_dir as an arg
fi


