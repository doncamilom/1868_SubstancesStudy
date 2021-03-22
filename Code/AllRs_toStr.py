#! /usr/bin/env python3

# Instructions: Run as
# ./I_getTabs.py <datafile.tsv> MaxNRowsMP NC
# MaxNRowsMP=Maximum number of rows of the data being transfered among subprocesses. It should be about 1000000
# NC=Number of compunds to use from dataset

# This program gets all Rs found (dirty). 
# Hopefully returns as well a list of these splitted in chunks

from scipy import sparse as sp
from time import time
import multiprocessing as mp
from itertools import chain
import sys
import os

def main(t0):
    global elemList
    sz = int(sys.argv[1])    

    print('\nLoading AllRs_dirty')

    Rs = sp.load_npz('./Data/AllRs_dirty.npz')
    print(f'\nAllRs_dirty loaded. Time: {time()-t0}')
    
    elemList = []
    with open(f"../Code/Data/ElementList.txt",'r') as f:
        for line in f:
            elemList.append(line.strip())

    chunksz = Rs.shape[0]//sz 
    Rs_list = []
   
    for i in range(sz):
        init = chunksz*i 
        end  = chunksz*(i+1) 
        if i<sz-1:      Rs_list.append(Rs[init:end]) 
        else:      Rs_list.append(Rs[init:]) 

    with mp.Pool(processes=sz) as pool:
        R_results = [pool.apply_async(getFormulas,args=(r,i,))
                     for i,r in enumerate(Rs_list)]        

        Rs_get = [r.get() for r in R_results]
    

   
#######################
# Auxiliary functions #
#######################

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
                if n > 0:               this_form += elemList[ind] 
                if n > 1:               this_form += str(int(n))
    
            f.write(this_form + f'X{int(this_data[-1]) if this_data[-1]!=1 else ""}\n')


if __name__ == '__main__':
    t0 = time()

    print(f"Entering program\n\n")
    main(t0)
    print(f"\n\nTotal runtime: {time()-t0:.3f} s")
