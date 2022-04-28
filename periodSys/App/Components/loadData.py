#! /usr/bin/python3

"""
Load all necessary data for the visualizations.
This includes:
    - Similarity Matrices
    - Optimized Permutations of elements
    - Families of elements
etc.
"""

import numpy as np
import pickle

root='../'
dataPath = root+'Data/'
resultsPath = root+ 'Results/'


def costPerm(P, perm):
    i = perm.reshape(1,-1)
    diffs = np.abs(i - i.T) + np.eye(i.shape[1])
    F = -np.sum(P/diffs)
    return F

def getElemList(dataPath):
    elemList = []
    with open("{}ElementList.txt".format(dataPath),'r') as f:
        for line in f:
            elemList.append(line.strip())
    return elemList

# Load element list
elemList = getElemList(dataPath)

# Load similarity matrices
simMats = np.load(dataPath + 'simMat.npy')
min_yr, max_yr = 1800,2022

# Load optimized permutations of elements
fh = open(resultsPath + 'optim_permut_yearly.gen', 'rb') 
opt_permut = pickle.load(fh)

# Load families of elements.
fh = open(resultsPath + 'foundGroups.bin','rb')
FEs = dict(pickle.load(fh))
