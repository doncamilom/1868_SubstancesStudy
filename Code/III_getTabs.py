#! /usr/bin/env python3

# Instructions: Run as
# ./II_getTabs.py NP MaxNRowsMP LCP
# NP=Number of processors to run in parallel
# NC=Number of compunds to use from dataset
# LCP=Last chunk proccessed. Restart non-unique R finding starting from this chunk.
# This script continues the work of the first part, which calculates all initial sparse mats.

import multiprocessing as mp
import numpy as np
from scipy import sparse as sp
from time import time
from itertools import chain
from makeVecs_sparse import allVecs_sparse
import pickle
import bz2
import sys
import os

def main():
    global NMax, size, maxLenArray

    size = int(sys.argv[1])
#    maxLenArray = int(sys.argv[2])  # If Rs list was splitted already, this arg is useless.
    last_chunk_processed = int(sys.argv[2])

    t0=time()

    # Load compounds
    cmpnds = sp.load_npz('./Data/Cmpnds_sparse.npz')
    sfile = bz2.BZ2File('./Data/year_ID_elems_nmax.bin', 'r')
    years,subsID, FullElemntList , NMax = pickle.load(sfile)

    #Construct a dict to go from element symbol to an index
    elemDict = {}
    for i,elem in enumerate(FullElemntList):
        elemDict[i] = elem
    

    ##########
    ## This depends on the output of the initial part. If the splitting was successful, then load that instead of the Rs.
    ##########

    # Load all produced Rs (dirty) : IF splitting was not successful, aka. './Data/List_Splitted_Rs_dirty.bin' doesn't exist
#    R_sparse = sp.load_npz('./Data/AllRs_dirty.npz')

    # Load the list of R chunks    : IF splitting was successful, aka.  './Data/List_Splitted_Rs_dirty.bin' does exist
    sfile = bz2.BZ2File('./Data/List_Splitted_Rs_dirty.bin', 'r')
    Rs_list = pickle.load(sfile)


    writeLogs(f"\n * Everything loaded in: {time()-t0:.3f} s")

    t0 = time()
    writeLogs("\nFinding unique Rs...\n")

#    Rs_list = unique_mp(R_sparse,size,maxLenArray,time()) # If first case: splitting not successful
    Rs_list = unique_mp_case2(Rs_list,size,time(),last_chunk_processed) # If second case: splitting successful

    ###############
    #### Program execution ends here. Let's see what data we get before doing anything else.
    ###############

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
    with open('./Data/logs_sparse_getTabs.log','a') as f:
        f.write(string)
         

def findRs(cmpnds):
    indices = cmpnds.indices
    data = cmpnds.data
    indptr = cmpnds.indptr
    
    sz_cols = cmpnds.shape[1] # Number of elements
    
    Rs = []  
    cmpnd_num_Rs = []
    
    for c in range(indptr.shape[0]-1):
        indx = indices[indptr[c]:indptr[c+1]]
        sub_dat = data[indptr[c]:indptr[c+1]]
        
        for i in range(indx.shape[0]):
            c_data = sub_dat.copy()
            n = int(c_data[i])          
            
            for j in range(n):   #Loop through ith element's subindex
                c_data[i] -= 1   #Remove one
                n = j+1          #How many atoms of this element have been removed 

                #Append compound data with a reduced entry (R-X n-1)
                Rs.append(  np.append(c_data.copy(),n)  )

        cmpnd_num_Rs.append((sub_dat.sum(),indx))
   
    # Construct sparse matrix
    for_indics = list(chain(*[l[0]*[l[1]] for l in cmpnd_num_Rs]) )

    indptr = np.cumsum([0]+list(map(lambda x: len(x)+1 , for_indics)))
    indices = np.array(list(chain(*[list(l)+[sz_cols] for l in for_indics])))
    data = np.array(list(chain(*[l for l in Rs])))

    Rs = sp.csr_matrix((data, indices, indptr),
                        shape=(len(Rs), sz_cols+1),
                        dtype=np.short)

    return Rs

def validRs(Rs,ids,start_from_id=0):
    """Get all non-unique R-n vectors out of a chunk Rs."""
    if ids>=start_from_id:
    
        shap1 = Rs.shape[0]
        new_rs , c = np.unique(Rs.toarray(),axis=0, return_counts=True)
        
        if (c>1).sum()!=0:    
            new_rs = sp.csr_matrix(new_rs[c > 1],dtype=np.short) 
            writeLogs(f"\t{shap1}\t{new_rs.shape[0]}\n")

            ## Save as results are produced.
            sp.save_npz(f"./Data/CleanRs_{ids}.npz",new_rs)
            return new_rs

def unique_mp(Rs,size,N,t0):
    """Create data chunks for finding non-unique Rs in parallel.
    Test how this does on full DS. If pickling errors, implement recursive spliting.
    """
    dist_list = split_chunk(Rs,0,N=N)
    writeLogs(f"\nRs successfully splited in {len(dist_list)} chunks. Time: {time()-t0} s")
    writeLogs(f"\n Shape1\t Shape2\n")

    # Save this list so we can start with this from the start
    sfile = bz2.BZ2File('./Data/List_Splitted_Rs_dirty.bin', 'w')
    pickle.dump(dist_list, sfile)

    with mp.Pool(processes=size) as pool:
        R_results = [pool.apply_async(validRs,args=(r,ids,))
                     for ids,r in enumerate(dist_list)]        

        Rs_get = [r.get() for r in R_results]

    Rs_get = [r for r in Rs_get if r is not None]

    return Rs_get

def unique_mp_case2(dist_list,size,t0,last_chunk_processed):
    """Create data chunks for finding non-unique Rs in parallel.
    Test how this does on full DS. If pickling errors, implement recursive spliting.

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

def split_chunk(chunk,i,N=1e6):
    """Create data chunks so that resulting chunks occupy much less memory.
    Go through each of the created chunks and, if any is above N rows, further split it.
    Define recursive function to further split chunks"""

    maxSplit = 5  # Maximum number of sub-splits 

    # Split using ith index:
    split = []
    step = [0,1,2,4,6,10]

    for j in range(maxSplit):
        lower,upper = step[j],step[j+1] 
        ith_col = chunk[:,i].toarray()[:,-1]

        if j < maxSplit-1:         tmp_splt = chunk[(ith_col>=lower) & (ith_col<upper)]  # Entries that are either j or j+1
        else:                      tmp_splt = chunk[ith_col>=lower]   # Entries that are maxSplit-1 or larger
        split.append(tmp_splt)

    # Now recursively further split here
    newList = []
    for l in split:
        if l.shape[0] > N and i+1 < l.shape[1]:
            newList = newList + split_chunk(l,i+1,N)  # Further split l by next i
        elif l.shape[0]>1:      newList.append(l)  # Append only if chunk contains more than one entry

    return newList

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
    if 'logs_sparse_getTabs.log' in os.listdir('./Data/'):    os.remove('./Data/logs_sparse_getTabs.log')
    
    Matches,Rs = main()
    writeLogs(f"\n\nTotal runtime: {time()-t0:.3f} s")

    sfile = bz2.BZ2File('./Data/AllMatches.bin', 'w')
    pickle.dump(Matches, sfile)

    sp.save_npz('./Data/Rs_clean.npz',Rs)
