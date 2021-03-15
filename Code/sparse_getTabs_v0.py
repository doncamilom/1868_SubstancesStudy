#! /usr/bin/env python3

# Instructions: Run as
# ./MPIgetTabs_v2.py <datafile.tsv> NP MaxNRowsMP NC
# NP=Number of processors to run in parallel
# NC=Number of compunds to use from dataset
# MaxNRowsMP=Maximum number of rows of the data being transfered among subprocesses. It should be about 1000000

import multiprocessing as mp
import numpy as np
from scipy import sparse as sp
from time import time
from itertools import chain
from makeVecs_sparse import allVecs_sparse
import pickle
import bz2
import sys

def main():
    global NMax, size, maxLenArray

    DataFile = sys.argv[1]
    size = int(sys.argv[2])
    maxLenArray = int(sys.argv[3])
    if len(sys.argv)<5:     NMax = None   # If no max file length is passed, default to all
    else:                   NMax = int(sys.argv[4])

    # Preprocess compound data (make cmpnd vecs, years and ID) + produce element list.
    cmpnds,years,subsID, FullElemntList , NMax = allVecs_sparse(DataFile,NMax) 

    #Construct a dict to go from element symbol to an index
    elemDict = {}
    for i,elem in enumerate(FullElemntList):
        elemDict[i] = elem
    
    avgN = np.mean(np.sum(cmpnds,axis=1))
    print(f"\n * Average num. of atoms per compound: {round(avgN,4)}")
    print(f" * Max number of fragments R: {avgN*cmpnds.shape[0]:.0f}")
    print(f"\n * {cmpnds.shape[0]} unique compounds out of {NMax} provided...")
    print(f"\t Which means we dropped {NMax-cmpnds.shape[0]} compounds for being non-stoichiometric or being repeated (isomers).\n")

    t0=time()
    R_sparse = findRs(cmpnds)  # Returns a scipy.sparse.csr_matrix containing all possible (R,n)s.

    print(f" * All possible Rs were produced in: {time()-t0:.3f} s")

    t0 = time()
    print("Finding unique Rs...\n")
    Rs_list = unique_mp(R_sparse,size,maxLenArray)

    sz = sum([i.shape[0] for i in Rs_list])     # Calculate total amount of Rs.
    Rs_list = rejoin_chunks(Rs_list,sz)         # Rejoin some chunks

    t = time()-t0
    print(f"\n * Found a total of {sz} non-unique Rs in {t:.3f} s")
    print(f"\nStarting finding commonalities.\n")
    t = time()
    
    with mp.Pool(processes=size) as pool:
        R_results = [pool.apply_async(get_matches,args=(rs,cmpnds,years,subsID,elemDict,))
                     for rs in Rs_list ]
        R_get = [r.get() for r in R_results]
    
    Matches = list(chain(*R_get))
    Rs = sp.vstack(Rs_list)

    print("\nDone!")
    print(f"Time: {time()-t0:.3f} s\t NProc: {size}\t NMax: {NMax} ")


    ##### Calculate stats
    lensR = [len(x[0]) for x in Matches]
    meanLen = sum(lensR)/len(Matches)
    print(f"\tMean number of compounds per R: {meanLen}. Total R(n)s found: {len(Matches)}")  
    print(f"\tMax number of compounds per R (first 3 max) : {sorted(set(lensR))[-3:]}")
        
    return Matches,Rs


#######################
# Auxiliary functions #
#######################

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

def validRs(Rs):
    """Get all non-unique R-n vectors out of a chunk Rs."""
    new_rs , c = np.unique(Rs.toarray(),axis=0, return_counts=True)
    new_rs = new_rs[c > 1]
    if new_rs.shape[0] > 0:    return sp.csr_matrix(new_rs,dtype=np.short) 

def unique_mp(Rs,size,N):
    """Create data chunks for finding non-unique Rs in parallel.
    Test how this does on full DS. If pickling errors, implement recursive spliting.
    """
    dist_list = split_chunk(Rs,0,N=N)

    with mp.Pool(processes=size) as pool:
        R_results = [pool.apply_async(validRs,args=(r,))
                     for r in dist_list]        

        Rs_get = [r.get() for r in R_results]

    Rs_get = [r for r in Rs_get if r is not None]

    print(list(map(lambda x:x.shape[0], Rs_get)))

    return Rs_get

def split_chunk(chunk,i,N=1e6):
    """Create data chunks so that resulting chunks occupy much less memory.
    Go through each of the created chunks and, if any is above N rows, further split it.
    Define recursive function to further split chunks"""

    maxSplit = 3  # Maximum number of sub-splits 

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
        if i%1000==0:       print( f"\t{i}/{totalRs}\t R(n)s evaluated..." )
            
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
    Matches,Rs = main()
    print(f"\n\nTotal runtime: {time()-t0:.3f} s")

    sfile = bz2.BZ2File('./Data/AllMatches.bin', 'w')
    pickle.dump(Matches, sfile)

    sp.save_npz('./Data/AllRs.npz',Rs)
