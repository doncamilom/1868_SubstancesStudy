#! /usr/bin/env python3

# Instructions: Run as
# ./I_getTabs.py <datafile.tsv> NP NC
# NP=Number of procesors to use
# NC=Number of compunds to use from dataset

# This program gets all Rs found (dirty) and converts them into strings, producing a number 'NP' of files. 

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

def main(t0):
    global NMax,FullElemntList

    DataFile = sys.argv[1]
    sz = int(sys.argv[2])
    if len(sys.argv)<4:     NMax = None   # If no max file length is passed, default to all
    else:                   NMax = int(sys.argv[3])

    writeLogs('\nStarting run: Getting all compound data')
    # Preprocess compound data (make cmpnd vecs, years and ID) + produce element list.
    cmpnds,years,subsID, FullElemntList , NMax = allVecs_sparse(DataFile,NMax) 
    writeLogs(f'\nAll compounds loaded, saving... Total time: {time()-t0}')

    sp.save_npz('./Data/Cmpnds_sparse.npz',cmpnds)
    with bz2.BZ2File('./Data/year_ID_elems_nmax.bin', 'w') as f:
        pickle.dump([years,subsID,FullElemntList,NMax], f)

    avgN = np.mean(np.sum(cmpnds,axis=1))
    writeLogs(f"\n * Average num. of atoms per compound: {round(avgN,4)}")
    writeLogs(f"\n * Max number of fragments R: {avgN*cmpnds.shape[0]:.0f}")
    writeLogs(f"\n * {cmpnds.shape[0]} unique compounds out of {NMax} provided...")
    writeLogs(f"\n\t Which means we dropped {NMax-cmpnds.shape[0]} compounds for being non-stoichiometric or being repeated (isomers).\n")

    writeLogs(f"\n Starting finding all Rs. Total time: {time()-t0}")
    
    R_sparse = findRs(cmpnds)  # Returns a scipy.sparse.csr_matrix containing all possible (R,n)s.
    sp.save_npz('./Data/AllRs_dirty.npz',R_sparse)

    writeLogs(f"\n * All possible Rs were produced in: {time()-t0:.3f} s")

    ### Start converting newly found Rs into strs
    chunksz = R_sparse.shape[0]//sz 
    Rs_list = []
   
    for i in range(sz):
        init = chunksz*i 
        end  = chunksz*(i+1) 
        if i<sz-1:      Rs_list.append(R_sparse[init:end]) 
        else:      Rs_list.append(R_sparse[init:]) 

    writeLogs(f"\nWriting Rs as strings into {sz} files")

    with mp.Pool(processes=sz) as pool:
        R_results = [pool.apply_async(getFormulas,args=(r,i,))
                     for i,r in enumerate(Rs_list)]        

        Rs_get = [r.get() for r in R_results]
    writeLogs(f"\nAll Rs converted. Time = {time()-t0:.3f} s")


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

def getFormulas(Rs,i):
    """Convert R sparse vector into string composition:
    sparse([0,1,0,...,4,6]) --> Ti4X6 for instance """

    data = Rs.data
    indi = Rs.indices
    iptr = Rs.indptr
    

    # Remove file if exists already
    if f'strs_{i}_.txt' in os.listdir('./Data/'):       os.remove(f'./Data/strs_{i}_.txt')

    with open(f'./Data/strs_{i}_.txt','a') as f:
        for i in range(len(Rs.indptr)-1):
            init = iptr[i]
            end  = iptr[i+1]

            this_data = data[init:end]
            this_indi = indi[init:end]

            this_form = ''
            for ind,n in zip(this_indi[:-1],this_data[:-1]):
                if n > 0:               this_form += FullElemntList[ind] 
                if n > 1:               this_form += str(int(n))
    
            f.write(this_form + f'X{int(this_data[-1]) if this_data[-1]!=1 else ""}\n')

if __name__ == '__main__':
    t0 = time()

    # Remove log file if it existed before this execution
    if 'logs_I_getTabs.log' in os.listdir('./Data/'):    os.remove('./Data/logs_I_getTabs.log')
    
    writeLogs(f"Entering program\n\n")
    main(t0)
    writeLogs(f"\n\nTotal runtime: {time()-t0:.3f} s")
