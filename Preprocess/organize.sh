#! /usr/bin/env bash

mkdir Data
mv *tsv Data/
cp ../Data/ElementList.txt Data/
python3 dataToPy.py
