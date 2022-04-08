#!/usr/bin/bash

# Compute all reaction templates
./getRxnTemp.awk Data/reaction_sets.tsv > Data/subs_rxn_templates.tsv

# Group substances and compute set of reaction templates for each, in parallel
rm scr/*
parallel -j40 "echo {} > scr/f{}.tmp; ./subsDescr.awk scr/f{}.tmp Data/subs_rxn_templates.tsv >> scr/res{%}.tsv; rm scr/f{}.tmp" ::: $(seq -f "%04g" 0 9999)
cat scr/res*.tsv > Data/subst_descr.tsv
rm scr/*


