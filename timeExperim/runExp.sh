#! /usr/bin/env bash

for i in 100 500 1000 2000 2500 3000 3500 4000 
do
    mpiexec -n 2 python timeExperiment.py ../sample.csv $i >> ../Data/resultsTimeExperim_2P.txt
done

for i in 100 500 1000 2000 2500 3000 3500 4000 
do
    mpiexec -n 3 python timeExperiment.py ../sample.csv $i >> ../Data/resultsTimeExperim_3P.txt
done

for i in 100 500 1000 2000 2500 3000 3500 4000 
do
    mpiexec -n 4 python timeExperiment.py ../sample.csv $i >> ../Data/resultsTimeExperim_4P.txt
done
