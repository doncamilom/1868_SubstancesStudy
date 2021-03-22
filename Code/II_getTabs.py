#! /usr/bin/env python3

# Instructions: Run as
# ./II_getTabs.py NP LCP
# NP=Number of processors to run in parallel
# LCP=Last chunk proccessed. Restart non-unique R finding starting from this chunk.
# This script continues the work of the first part, which calculates all initial sparse mats.

import multiprocessing as mp
import numpy as np
from scipy import sparse as sp
from time import time
from itertools import chain
import pickle
import bz2
import sys
import os

def main():
    global size,last_chunk_processed

    size = int(sys.argv[1])
    if len(sys.argv)<3: last_chunk_processed = 0
    else : last_chunk_processed = int(sys.argv[2])

    t0=time()

    # Load the list of R chunks    : IF splitting was successful, aka.  './Data/List_Splitted_Rs_dirty.bin' does exist
    with bz2.BZ2File('./Data/Chunks_Rs_dirty.bin', 'r') as f:
        Rs_list = pickle.load(f)

    writeLogs(f"\n * Everything loaded in: {time()-t0:.3f} s")

    t0 = time()
    writeLogs("\nFinding unique Rs...\n")

    # Find non-unique Rs starting from chunk `last_chunk_processed`
    if not os.path.isdir('./Data/scr/'): os.mkdir('./Data/scr/') # If scr dir hasn't been created, create it
    Rs_list = unique_mp(Rs_list,size,time(),last_chunk_processed) 


#######################
# Auxiliary functions #
#######################
    
def writeLogs(string):
    with open(f'./Data/logs_II_getTabs_{last_chunk_processed}_.log','a') as f:
        f.write(string)

def validRs(Rs,ids,start_from_id=0):
    """Get all non-unique R-n vectors out of a chunk Rs."""
    if ids>=start_from_id:
    
        shap1 = Rs.shape[0]
        new_rs , c = np.unique(Rs.toarray(),axis=0, return_counts=True)
        
        if (c>1).sum()!=0:    
            new_rs = sp.csr_matrix(new_rs[c > 1],dtype=np.short) 
            writeLogs(f"\t{shap1}\t{new_rs.shape[0]}\n")

            ## Save as results are produced.
            sp.save_npz(f"./Data/scr/CleanRs_{ids}_.npz",new_rs)
            return new_rs

def unique_mp(dist_list,size,t0,last_chunk_processed):
    """Create data chunks for finding non-unique Rs in parallel.

    unique_from_id = last chunk id processed in a previous run. 
    Allows to restart computation from that id.
    """

    writeLogs(f"\n Shape1\t Shape2\n")

    with mp.Pool(processes=size) as pool:
        R_results = [pool.apply_async(validRs,args=(r,ids,last_chunk_processed,))
                     for ids,r in enumerate(dist_list)]        

        Rs_get = [r.get() for r in R_results]

    Rs_get = [r for r in Rs_get if r is not None]

    return Rs_get

if __name__ == '__main__':
    t0 = time()

    main()

    writeLogs(f"\n\nTotal runtime: {time()-t0:.3f} s")
