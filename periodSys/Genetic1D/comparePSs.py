#!/usr/bin/python3

import sys, os

# Add ../similarity to path, so we can import functions from simMat.py
path = os.path.join(os.path.dirname( __file__ ), '..', 'similarity')
sys.path.append(path)

from simMat import *
from genetic1D import *

global root


def main():
    global resultsPath, elemList, simMat_yr, min_yr, Indivs_yr

    resultsPath = root+'Results/'
    dataPath = root+'Data/'
    elemList = getElemList(dataPath)

    simMat_yr = np.load(dataPath+'simMat.npy',allow_pickle=True)

    fh = open(resultsPath + 'optim_permut_yearly.gen', 'rb')
    Indivs_yr = pickle.load(fh)

    min_yr = 1800

    historyMatrix(p=2, r=4)


def historyMatrix(p, r):
    """
    This is a highly vectorized version of `comparePerms`.
    In each iteration, compares a single permutation, simultaneously with one for each year (222 total).

    r is radius to consider elements to be close in a sequence.
    p is number of permutations to use. Use only the best `p` found each year.
    """

    # Each year, sort optimized permutations by cost
    Indivs_yr_sort = {}
    for yr in Indivs_yr.keys():
        Indivs_yr_sort[yr] = []

        S = simMat_yr[yr - min_yr].copy()
        P = symmetrize(S)

        a = np.array([cost(P,ind) for ind in Indivs_yr[yr]])

        # Sort the list of Individuals for each year, based on cost
        order = np.argsort(a)
        Indivs_yr_sort[yr] = np.array(Indivs_yr[yr])[order]


    # Calc if diag entry of element `i` is > 0 in simMat for every year
    exist_elems = np.einsum('ijj->ij',simMat_yr)>0

    # Convert to array of dim (222, 50, 103, 1)
    a = np.array([Indivs_yr_sort[y] for y in range(1800,2022)])[...,np.newaxis]

    # Get matrix of distances between each pair of elements, in each permutation
    a = np.abs(a - a.reshape(222, 50, 1, 103))

    # Get only upper diag of matrix, so combs. of elements are not repeated
    a = np.triu(a,k=1)

    # And convert these elements two high values, so don't count in sum
    a[a==0] = 200

    # Initialize matrix of historical comparisons
    matrix = np.zeros((222,222))

    # Iterate over year, and two indexes of permutation.
    for i,y in enumerate(range(1800,2022)):
        for j in range(p):
            for k in range(p):

                # Mask elements that don't exist this year
                mask = exist_elems[i]
                mask = ~(mask * mask.reshape(-1,1)) * 1
                mask[mask==1] = 200 # mask[i,j] > 200 if either ith or jth element don't exist, else 0

                # Sum mask. If elements exist, values are not modified. Thus if val < r, condition is satisfied.
                l1 = (a[i,j] + mask) <= r
                l2 = (a[:,k] + mask) <= r

                sim = np.sum(l1*l2, axis=(1,2))

                # Normalize using one of this year's permutation
                nor = np.sum(l1*l1) # Compare perm to itself

                matrix[i] += sim/nor  # Sum now, divide by number of combinations later


    # This division gets the average.
    matrix /= p**2
    np.save(root+"Results/mathistory.npy",matrix)
    return matrix


if __name__=='__main__':
    global root
    import argparse
    parser = argparse.ArgumentParser(description="Compute similarity matrices.")

    parser.add_argument("--inp_dir","-i", required=False, type=str, default="./",
                        help="Root directory for all computations")
    args = parser.parse_args()
    root=args.inp_dir

    from time import time
    t0=time()
    main()

    print("Total time: {:.3f} s".format(time()-t0))
