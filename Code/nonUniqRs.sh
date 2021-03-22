#! /usr/bin/env bash

# Concatenate all resulting files

cat ./Data/strs_*_.txt > ./Data/AllRs_strs.txt

awk '
BEGIN {}
{if ($1 in uniq)
    nonuniq[$1];
else
    uniq[$1];
}
END{ for (c in nonuniq) print c   }
' ./Data/AllRs_strs.txt > ./Data/AllRs_clean_strs.txt
