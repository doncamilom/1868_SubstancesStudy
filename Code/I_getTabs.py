#! /usr/bin/env python3

# Instructions: Run as
# ./I_getTabs.py <datafile.tsv> MaxNRowsMP NC
# MaxNRowsMP=Maximum number of rows of the data being transfered among subprocesses. It should be about 1000000
# NC=Number of compunds to use from dataset

# This program gets all Rs found (dirty). 
# Hopefully returns as well a list of these splitted in chunks

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
    global NMax, maxLenArray

    DataFile = sys.argv[1]
    maxLenArray = int(sys.argv[2])
    if len(sys.argv)<4:     NMax = None   # If no max file length is passed, default to all
    else:                   NMax = int(sys.argv[3])

    # Preprocess compound data (make cmpnd vecs, years and ID) + produce element list.
    cmpnds,years,subsID, FullElemntList , NMax = allVecs_sparse(DataFile,NMax) 

    sp.save_npz('./Data/Cmpnds_sparse.npz',cmpnds)
    with bz2.BZ2File('./Data/year_ID_elems_nmax.bin', 'w') as f:
        pickle.dump([years,subsID,FullElemntList,NMax], f)

    avgN = np.mean(np.sum(cmpnds,axis=1))
    writeLogs(f"\n * Average num. of atoms per compound: {round(avgN,4)}")
    writeLogs(f"\n * Max number of fragments R: {avgN*cmpnds.shape[0]:.0f}")
    writeLogs(f"\n * {cmpnds.shape[0]} unique compounds out of {NMax} provided...")
    writeLogs(f"\n\t Which means we dropped {NMax-cmpnds.shape[0]} compounds for being non-stoichiometric or being repeated (isomers).\n")

    t0=time()
    R_sparse = findRs(cmpnds)  # Returns a scipy.sparse.csr_matrix containing all possible (R,n)s.
    sp.save_npz('./Data/AllRs_dirty.npz',R_sparse)

    writeLogs(f"\n * All possible Rs were produced in: {time()-t0:.3f} s")

    t0 = time()
    writeLogs("\nFinding unique Rs...\n")

    Rs_list = split_chunk(R_sparse,0,maxLenArray)

    with bz2.BZ2File('./Data/Chunks_Rs_dirty.bin', 'w') as f:
        pickle.dump(Rs_list, f)
    
    writeLogs(f"\nRs successfully splited in {len(Rs_list)} chunks. Time: {time()-t0} s")


#######################
# Auxiliary functions #
#######################
    
def writeLogs(string):
    with open('./Data/logs_I_getTabs.log','a') as f:
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


if __name__ == '__main__':
    t0 = time()

    # Remove log file if it existed before this execution
    if 'logs_I_getTabs.log' in os.listdir('./Data/'):    os.remove('./Data/logs_I_getTabs.log')
    
    main()
    writeLogs(f"\n\nTotal runtime: {time()-t0:.3f} s")
