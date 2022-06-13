#!/bin/bash

# Parallel process all raw XML files, and dump results into scr
for y in $(seq 0 19)
do
	parallel -j60 --max-args 500 "./use_subs_single.awk {} >> scr/count_use{%}.txt" ::: ../../DATA/RXN/p$y/*
	echo "Done with data from p$y"
done

