#!/bin/bash

# Parallel process all raw XML files, and dump results into scr
for y in $(seq 0 19)
do
	parallel -j60 --max-args 500 "./new_rxns_docs_auth.awk {} >> scr/yr_rx_cnr_auts_dt_{%}.txt" ::: ../../DATA/RXN/p$y/*
	echo "Done with data from p$y"
done

echo "Done with raw data. Now collecting results."

cat scr/* > rxns_data.txt

