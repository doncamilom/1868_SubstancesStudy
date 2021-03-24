#! /usr/bin/env bash

# Concatenate all resulting files

cat ./Data/strs_*_.txt > ./Data/AllRs_strs.txt
rm ./Data/strs_*_.txt

if [ -f "./Data/AllRs_clean_strs.txt" ]; then
    echo "File of clean Rs already exists. Deleting it...";
    rm "./Data/AllRs_clean_strs.txt";
fi

awk '
BEGIN {}
{if ($1 in uniq)
    nonuniq[$1];
else
    uniq[$1];
}
END{ for (c in nonuniq) print c   }
' ./Data/AllRs_strs.txt > ./Data/AllRs_clean_strs.txt
