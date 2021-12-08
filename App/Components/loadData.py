#! /usr/bin/env python3

import pickle

## Load all groups generated with script
grps={}
for yr in range(1800,2018,2):
    try:
        fh = open(f'../CoreCalcs/foundGroups/grp{yr}.bin','rb')
        grps[yr] = pickle.load(fh)[1]
    except:        pass#print(f"{yr} not found")

grps[1864] = grps[1862] # Just copy results from previous year as there was an error  
