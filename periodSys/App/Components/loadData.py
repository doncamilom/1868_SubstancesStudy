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

root='../'
dataPath = root+'Data/'
resultsPath = root+ 'Results/'


def getElemList(dataPath):
    elemList = []
    with open("{}/ElementList.txt".format(dataPath),'r') as f:
        for line in f:
            elemList.append(line.strip())
    return elemList

# Load element list
elemList = getElemList(dataPath)
simMats = np.load(dataPath + 'simMat.npy')
min_yr, max_yr = 1800,2022


