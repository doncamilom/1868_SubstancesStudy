#! /usr/bin/env python3

# Instructions: Run as
# ./III_getTabs.py NP
# NP=Number of processors to run in parallel
# This script continues the work of the second part, where all non-unique Rs were found.

import multiprocessing as mp
import numpy as np
from scipy import sparse as sp
from time import time
from itertools import chain
import pickle
import bz2
import sys
import os
import glob


def main():
    global size

    size = int(sys.argv[1])

    t0=time()

    # Load compounds
    cmpnds = sp.load_npz('./Data/Cmpnds_sparse.npz')
    with bz2.BZ2File('./Data/year_ID_elems_nmax.bin', 'r') as f:
        years,subsID, FullElemntList , NMax = pickle.load(f)

    #Construct a dict to go from element symbol to an index
    elemDict = {i:elem for i,elem in enumerate(FullElemntList)}

    ###############
    # From part II, we got a bunch of files named `` in Data/, which contain all Rs of interest
    # So we have to load all of them.
    ###############
    file_list = glob.glob('Data/scr/CleanRs*.npz')
    Rs_list = [sp.load_npz(f) for f in file_list]

    writeLogs(f"\n * Everything loaded in: {time()-t0:.3f} s")

    t0 = time()
    sz = sum([i.shape[0] for i in Rs_list])     # Calculate total amount of Rs.
    Rs_list = rejoin_chunks(Rs_list,sz)         # Rejoin some chunks

    t = time()-t0
    writeLogs(f"\n * Found a total of {sz} non-unique Rs in {t:.3f} s")
    writeLogs(f"\nStarting finding matchings.\n")
    t = time()
    
    with mp.Pool(processes=size) as pool:
        R_results = [pool.apply_async(get_matches,args=(rs,cmpnds,years,subsID,elemDict,))
                     for rs in Rs_list ]
        R_get = [r.get() for r in R_results]
    
    Matches = list(chain(*R_get))
    Rs = sp.vstack(Rs_list)

    writeLogs("\nDone!")
    writeLogs(f"\nTime: {time()-t0:.3f} s\t NProc: {size}\t NMax: {NMax} ")


    ##### Calculate stats
    lensR = [len(x[0]) for x in Matches]
    meanLen = sum(lensR)/len(Matches)
    writeLogs(f"\n\tMean number of compounds per R: {meanLen}. Total R(n)s found: {len(Matches)}")  
    writeLogs(f"\n\tMax number of compounds per R (first 3 max) : {sorted(set(lensR))[-3:]}")
        
    return Matches,Rs


#######################
# Auxiliary functions #
#######################
    
def writeLogs(string):
    with open('./Data/logs_III_getTabs.log','a') as f:
        f.write(string)
         
def rejoin_chunks(Rs_uniq,sz):
    """ Rejoin some of the chunks produced.
    Some are very small and it's inefficient to call one single process for such small chunk"""
    final_chunk_sz = sz//size
    
    Rs_lists = []    
    curr_pointer = 0
    exhausted = False
    for i in range(size):
        if i+1<size:
            curr_len = 0
            curr_list = []

            for j,l in enumerate(Rs_uniq[curr_pointer:]):
                if curr_len < final_chunk_sz:
                    curr_list.append(l)
                    curr_len += l.shape[0]
                    exhausted = True
                else:
                    curr_pointer += j
                    exhausted = False
                    break
            Rs_lists.append(sp.vstack(curr_list))  # Make a new concatentated chunk of approx size of totalRs/size 
            if exhausted: break

        elif not exhausted:    Rs_lists.append(sp.vstack(Rs_uniq[curr_pointer:]))  # Append a chunk of anything left

    return Rs_lists

def get_matches(Rs,cmpnds,years,subsID,elemDict):
    """Get matches. For each R(n) find all elements X such that compound R-Xn exists in dataset. 
    Build an element set for each R(n).
    """
    totalRs = Rs.shape[0]
    ns = Rs[:,-1].data
    R  = Rs[:,:-1]
    
    sumCmpnds = cmpnds.sum(axis=1)
    sumRaxis1 = np.array(    R.sum(axis=1).flatten() + ns    ).flatten()
    
    Matches = []
    for i,n in enumerate(ns):
        if i%1000==0:       writeLogs( f"\n\t{i}/{totalRs}\t R(n)s evaluated..." )
            
        r = R[i] #The actual R
        
        """Encode a condition to search only within a subset of compounds
        fulfulling certain conditions based on R"""
        # 1. R is contained in compound        
        cond1 = ((cmpnds - r.toarray())>=0).all(axis=1)
        # 2. sum of atoms in cmpnd == sum of atoms in R_ (sum(R) + n)
        cond2 = (sumCmpnds == sumRaxis1[i])
        
        cond = np.array(cond1 & cond2).flatten()  # Combine conditions
        subsetCmpnds = cmpnds[cond]  # Select subset of cmpnds
        curr_years = years[cond]
        curr_subsID = subsID[cond]
        
        cmpnds_no_R = (subsetCmpnds - r.toarray())
        
        # Now select only those cmpnds where residual is due to one element only (X_n)
        # Only useful for n!=1
        if n!=1:
            cond = np.array((cmpnds_no_R!=0).sum(axis=1)==1).flatten()
            subsetCmpnds = subsetCmpnds[cond]
            curr_years = curr_years[cond]
            curr_subsID = curr_subsID[cond]
            
        # At this point, subsetCmpnds contains all compounds that match with R(n).
        elemIndex = (subsetCmpnds - r.toarray()).nonzero()[1]
        curr_list = list(map(lambda x: elemDict[x]  , elemIndex))  # Map dict to above list of elems

        Matches.append( [curr_list, curr_years, curr_subsID] )

        ###########
        ## Deal with this after you get all data created above
        #Table_list.append(getTable(curr_list,curr_years,curr_subsID,useID=useID))      
        ###########
    return Matches      # [list_elems,list_years, list_ids] for each R(n)


if __name__ == '__main__':
    t0 = time()

    # Remove log file if it existed before this execution
    if 'logs_III_getTabs.log' in os.listdir('./Data/'):    os.remove('./Data/logs_III_getTabs.log')
    
    Matches,Rs = main()
    writeLogs(f"\n\nTotal runtime: {time()-t0:.3f} s")

    sfile = bz2.BZ2File('./Data/AllMatches.bin', 'w')
    pickle.dump(Matches, sfile)

    sp.save_npz('./Data/Rs_clean.npz',Rs)
