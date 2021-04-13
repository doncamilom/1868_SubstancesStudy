#! /usr/bin/env python3

import pickle
import bz2
from itertools import chain
import numpy as np
import scipy.sparse as sp 

def main():
    # Load element list
    elemList = []
    with open(f"./Data/ElementList.txt",'r') as f:
        for line in f:
            elemList.append(line.strip())
    elemList.append('X')

    # Load clean Rs
    with open("./Data/allRs_clean.tsv","r") as f:
        lines = f.readlines()

    matches, Rs = [],[]

    for l_ in lines:
        l = l_.split("\n")[0].split("\t")

        # Convert Rs (str) to Rs (sparse matrix)
        Rs.append(l[0])
        
        matches.append([l[1].split(":"),l[3].split(":"),l[2].split(":")])

    sparseRs = Rs_str_to_sparse(Rs,elemList)

    # Save results
    sfile = bz2.BZ2File('./Data/matches_pickle.bin', 'w')
    pickle.dump(matches,sfile)

    # Save sparse matrix of Rs
    sp.save_npz('./Data/Rs_sparse.npz',sparseRs)


def getVec(strRs,elemList):
    Li = strRs.split('-') 
    li = [fs.split(':') for fs in Li]
        
    # Construct two lists: input for sparse matrix construction
    col  = [elemList.index(i[0]) for i in li]  # Index of element i to put correspondent data
    data = [int(i[1]) for i in li]           # Num. atoms of element i

    return col,data

def Rs_str_to_sparse(strRs,elemList):
    # List of lists [col,data]
    colXdata = list(map(lambda x: getVec(x,elemList) , strRs))
    colXdata = [l for l in colXdata]
    
    indptr = np.cumsum([0]+list(map(lambda x: len(x[0]) , colXdata)))
    indices = np.array(list(chain(*[l[0] for l in colXdata])))
    data = np.array(list(chain(*[l[1] for l in colXdata])))
    
    sparseRs = sp.csr_matrix((data, indices, indptr),
                           shape=(len(colXdata), len(elemList)),
                           dtype=np.short)

    return sparseRs


if __name__=='__main__':

    main()
