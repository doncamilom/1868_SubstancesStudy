#!/usr/bin/bash

# Filter between organic and inorganic
./split.awk ../Data/format_MFs.tsv > filtered.tsv

# Split into two different files
csplit filtered.tsv /^==*/ {*} -z --suppress-matched

# Rename files
mv xx00 inorg_subs.tsv
mv xx01 org_subs.tsv
