#!/usr/bin/python3

import sys, os

# Add ../similarity to path, so we can import functions from simMat.py
path = os.path.join(os.path.dirname( __file__ ), '..', 'similarity')
sys.path.append(path)

from simMat import *
from genetic1D import *

global NP,root


def main():
    global resultsPath, elemList, simMat_yr, min_yr
    resultsPath = root+'Results/'
    dataPath = root+'Data/'
    elemList = getElemList(dataPath)
    simMat_yr = np.load(dataPath+'simMat.npy',allow_pickle=True)
    min_yr = 1800

    #baselinesGA(2020)
    History(percentil=0.5)

def baselinesGA(year):
    """
    Calculate costs for some baseline orderings. 
    e.g. AN, Original Pettifor, Glawe GA.
    """

    # Calculate P for calculating cost function
    S = simMat_yr[year - min_yr].copy()
    P = symmetrize(S)

    # Set a reference ordering: atomic number
    orderAO = np.arange(P.shape[0])

    print("Cost using AN order = {:.5f}".format(cost(P, orderAO)))

    # 2. Check how Pettifor's scale does on these data.
    pettif = ['He','Ne','Ar','Kr','Xe','Rn','Fr','Cs','Rb','K','Na','Li','Ra','Ba','Sr','Ca','Yb','Eu','Y',
              'Sc','Lu','Tm','Er','Ho','Dy','Tb','Gd','Sm','Pm','Nd','Pr','Ce','La','Lr','No','Md','Fm','Es',
              'Cf','Bk','Cm','Am','Pu','Np','U','Pa','Th','Ac','Zr','Hf','Ti','Nb','Ta','V','Mo','W','Cr',
              'Tc','Re','Mn','Fe','Os','Ru','Co','Ir','Rh','Ni','Pt','Pd','Au','Ag','Cu','Mg','Hg','Cd','Zn',
              'Be','Tl','In','Al','Ga','Pb','Sn','Ge','Si','B','Bi','Sb','As','P','Po','Te','Se','S','C','At',
              'I','Br','Cl','N','O','F','H']

    # List the position of each element in Pettifor scale
    order_pett = np.array([pettif.index(e) for e in elemList])
    # e.g. H: 102, He: 0, Li: 11, etc
    print("Cost Pettifor = {:.5f}".format(cost(P, order_pett)))

    # See how their GA solution works

    GA_ref = ['He','Ne','Ar','At','Rn','Fr','Es','Fm','Md','No','Lr','Kr','Xe','Pm','Cs','Rb','K','Na',
              'Li','Ra','Ba','Sr','Ca','Eu','Yb','Lu','Tm','Y','Er','Ho','Dy','Tb','Gd','Sm','Nd','Pr',
              'Ce','La','Ac','Am','Cm','Bk','Cf','Pu','Np','U','Th','Pa','Sc','Zr','Hf','Ti','Nb','Ta',
              'V','Cr','Mo','W','Re','Tc','Os','Ru','Ir','Rh','Pt','Pd','Au','Ag','Cu','Ni','Co','Fe',
              'Mn','Mg','Zn','Cd','Hg','Be','Al','Ga','In','Tl','Pb','Sn','Ge','Si','B','C','N','P','As',
              'Sb','Bi','Po','Te','Se','S','O','I','Br','Cl','F','H']

    # List the position of each element in Pettifor scale
    order_GA = np.array([GA_ref.index(e) for e in elemList])
    # e.g. H: 102, He: 0, Li: 11, etc
    print("Cost GA Glawe = {:.5f}\n".format(cost(P, order_GA)))


    # Performance of random configurations
    rand_ord_cost = [cost(P, np.random.permutation(orderAO)) for i in range(5000)]
    print(pd.Series(rand_ord_cost).describe())

def genToElem(gen):
    order = ['_' for i in range(103)]
    for i,idx in enumerate(gen):
        order[idx] = elemList[i]
    return order

def getNGrams(seq,N):
    """Returns a list of NGrams (sets) in seq"""
    ngrams = []
    for i in range(len(seq)-N+1):
        ngrams.append(set(seq[i:i+N]))
    return ngrams

def compareNGrams(seq1,seq2,N):
    """Compare NGrams of seq1 and seq2
    Calculate number of NGrams of seq1 that are also in seq2, divided by number of NGrams in seq1
    """
    ng1 = getNGrams(seq1,N)
    ng2 = getNGrams(seq2,N)
    
    """An N-Gram from ng1 can be compared with more than one of ng2.
    Normalization constant is calculated by comparing a seq to itself"""
    Norm = {2:102,3:301,4:494,5:681}
    
    total = 0
    for i in ng1:
        for j in ng2:
            if len(i.intersection(j)) > 1: total += 1
    return total/Norm[N]

def History(percentil):
    global Indivs_yr, normConsts_sim

    fh = open(resultsPath + 'optim_permut_yearly.gen', 'rb') 
    Indivs_yr = pickle.load(fh)

    ## Here, we define a general similarity function between orderings
    # First, precalc. normalization constants
    normConsts_sim = {}

    # These depend on: lenSeq and N
    for lenSeq in range(10,105):
        for N in range(3,10):
            ng = getNGrams(genToElem(Indivs_yr[1820][0])[:lenSeq],N)    # Pick any example, only important thing is number of elements
            Norm = 0
            for i in ng:
                for j in ng:
                    if len(i.intersection(j)) >1:
                        Norm += 1 
            normConsts_sim["{}:{}".format(N,lenSeq)] = Norm


    mat = compareSimYearly_all(Indivs_yr,simMat_yr,percentil=percentil,ngSize=3)
    print(mat)
    np.save(root+"Results/mathistory.npy",mat)


# Define similarity function
def compareNGrams_asym(seq1,seq2,N,lenSeq1):
    """Compare NGrams of seq1 and seq2
    Calculate number of NGrams of seq1 that are also in seq2, divided by number of NGrams in seq1
        As some elements are lacking for some years, 
        compare only the patterns formed by existing elements (lenSeq1 first elems)
    """
    ng1 = getNGrams(seq1[:lenSeq1],N)
    ng2 = getNGrams(seq2,N)

    total = 0
    for i in ng1:
        for j in ng2:
            if len(i.intersection(j)) > 1:
                total += 1
    return total


def simYearRow(i,yr0,ngSize):
    """
    Compute similarity between the individuals of year `yr0` and all other years.
    This function computes the `i`th row of similarity matrix.
    """
    N = len(Indivs_yr.keys())
    yearMatrix = np.zeros((N,N))

    for j,yr1 in enumerate(Indivs_yr.keys()):
        if yr0<=yr1:  # Calc only half the matrix
            pairs = product(Indivs_yr_elemSeq[yr0],Indivs_yr_elemSeq[yr1])
            Nelemsyr0 = elem_count[yr0-min_yr]
            ocho = np.mean([compareNGrams_asym(s0,s1,ngSize,Nelemsyr0) for s0,s1 in pairs]) 
            yearMatrix[i,j] = ocho

    return yearMatrix

from itertools import product
def compareSimYearly_all(Indivs_yr,simMat,percentil=0.3,ngSize=3):
    """
    Compute similarity between the PSs of all pairs of years,
    as the average similarity between all combinations of individuals from such years.
    """
    global Indivs_yr_elemSeq, elem_count

    # Make a copy of Indivs_yr, that contains element sequences.
    Indivs_yr_elemSeq = {}
    for yr in Indivs_yr.keys():
        Indivs_yr_elemSeq[yr] = []
        curr_ind = Indivs_yr[yr]

        # Get this year's top `percentile`% individuals only.
        S = simMat[yr - min_yr].copy()
        P = symmetrize(S)

        a = np.array([cost(P,ind) for ind in Indivs_yr[yr]])
        q = np.quantile(a,percentil)

        for i,seq in enumerate(np.array(Indivs_yr[yr])[a<=q]):
            Indivs_yr_elemSeq[yr].append(genToElem(seq))


    elem_count = (simMat_yr.sum(axis=1)>0).sum(axis=1)

    # Iterate through every pair of years. This is potentially paralelizable
    ######################
    import multiprocessing as mp
    
    with mp.Pool(processes=NP) as pool:
        results = [pool.apply_async(simYearRow,args=(i,yr0,ngSize,)) for i,yr0 in enumerate(Indivs_yr.keys())]
        res = [r.get() for r in results] 

    ######################
                
    yearMatrix = np.array(res).sum(axis=0)
    yearMatrix += yearMatrix.T

    # Set diag to half it's value
    d = np.diag_indices_from(yearMatrix)
    yearMatrix[d] /= 2
    yearMatrix[d] -= 1/len(Indivs_yr[1800])  #remove self-comparisons (1/number indivs per year) 

    # Build normalization matrix
    NormMat = [normConsts_sim["{}:{}".format(ngSize,elem_count[i-min_yr])] for i in Indivs_yr.keys()]
    NormMat = np.repeat(np.array(NormMat).reshape(-1,1),len(NormMat),axis=1)
    # Normalize
    yearMatrix /= NormMat

    return yearMatrix

if __name__=='__main__':
    global root
    import argparse
    parser = argparse.ArgumentParser(description="Compute similarity matrices.")

    parser.add_argument("--NumProc","-N", required=True, type=int,
                        help="Number of processors to use for computations.")
    parser.add_argument("--inp_dir","-i", required=False, type=str, default="./",
                        help="Root directory for all computations")
    args = parser.parse_args()
    NP=args.NumProc
    root=args.inp_dir

    from time import time
    t0=time()
    main()

    print("Total time: {:.3f} s".format(time()-t0))
