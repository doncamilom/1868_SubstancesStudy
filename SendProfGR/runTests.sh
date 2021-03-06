#! /usr/bin/env bash

# This script runs the code for variable data size and number of cores. Run as:
# ./runTests.sh {path/to/DataFile}.csv {NP}
# NP being the number of processors to use on this run

for i in 50000 100000 300000 500000 1000000
do
    ./MPIgetTabs_v2.py $1 $2 $i >> ./logs/NProc"$2"_NMaxData_"$i".log
done
