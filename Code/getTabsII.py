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
import re

def main():
    global size
    size = int(sys.argv[1])

    t0=time()

    # Load compounds
    cmpnds = sp.load_npz('./Data/Cmpnds_sparse.npz')
    with bz2.BZ2File('./Data/year_ID_elems_nmax.bin', 'r') as f:
        years,subsID, FullElemntList , NMax = pickle.load(f)

    # Load clean Rs and build sparse matrix (do this in parallel, so actually get a list of Rs chunks)
    Rs_list = load_Rs('./Data/AllRs_clean_strs.txt', FullElemntList)
    ## Save concat version of this
    sp.save_npz('./Data/AllRs_clean_sparse.npz',sp.vstack(Rs_list))
    #with bz2.BZ2File(, 'w') as f:
    #    pickle.dump(Rs_list,f)

    writeLogs(f"\nRs loaded and saved as sparse matrix. Time: {time()-t0:.3f} s")

    #Construct a dict to go from element symbol to an index
    elemDict = {i:elem for i,elem in enumerate(FullElemntList)}

    writeLogs(f"\n * Everything loaded in: {time()-t0:.3f} s")

    t0 = time()
    sz = sum([i.shape[0] for i in Rs_list])     # Calculate total amount of Rs.

    writeLogs(f"\n * Found a total of {sz} non-unique Rs in {time()-t0:.3f} s")
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

def get_sparse_vecs(tx, elemList):

    Li = re.split(r"(?<!^)(?=[A-Z])",tx)  #Split as ['H2','O']
    
    # Adds 1 if no subindex. Result is ['H2','O1']. 
    # Right after, split chem symbol from subindex as [['H',2],['O',1]]
    
    li = [re.split(r"([A-z]+)(([0-9]*[.])?[0-9]+)",i)
          if bool(re.match(r'[A-z]*([0-9]*[.])?[0-9]+',i))
          else re.split(r"([A-z]+)(([0-9]*[.])?[0-9]+)",i+'1') for i in Li]  
    
    # Construct two lists: input for sparse matrix construction
    col  = [elemList.index(i[1]) for i in li]  # Index of element i to put correspondent data
    data = [int(i[2]) for i in li]           # Num. atoms of element i
    return col,data


def Rs_str_to_sparse(rs_str,elemList):
    
    # List of lists [col,data]
    colXdata = list(map(lambda x: get_sparse_vecs(x,elemList) , rs_str))

    # See docs for scipy.sparse.csr_matrix to understand the syntaxis
    indptr = np.cumsum([0]+list(map(lambda x: len(x[0]) , colXdata)))
    indices = np.array(list(chain(*[l[0] for l in colXdata])))
    data = np.array(list(chain(*[l[1] for l in colXdata])))

    return sp.csr_matrix((data, indices, indptr), 
                         shape=(len(colXdata), len(elemList)),
                         dtype=np.short)

         
def load_Rs(path,elemList):
    
    with open(path,'r') as f:
        rs_strs = f.read().splitlines()

    chunk_sz = len(rs_strs)//size
    Rs_lists_strs = [rs_strs[chunk_sz*i:chunk_sz*(i+1)] if i<size-1 else rs_strs[chunk_sz*i:] for i in range(size)  ]
    
    # Process these strs in parallel to get sparse matrices
    _elemList = elemList + ['X']
    with mp.Pool(processes=size) as pool:
        R_results = [pool.apply_async(Rs_str_to_sparse,args=(rs,_elemList,))
                     for rs in Rs_lists_strs ]
        R_get = [r.get() for r in R_results]
    
    return R_get


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
