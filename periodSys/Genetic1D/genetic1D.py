#! /usr/bin/env python3

import numpy as np
from time import time
from TPs import TP
import multiprocessing as mp
import sys
import pickle
import os




def main():
    global P, elemList, root, dataPath,size


    print(root)
    dataPath = root+'Data/'
    resultsPath = root+ 'Results/'

    # Load element list
    elemList = getElemList(dataPath)
    simMat_yr = np.load(dataPath + 'simMat.npy')
    min_yr, max_yr = 1800,2022
    popSize, NGens = 1500, 200

    t0 = time()

    results_total = {}
    for yr in range(min_yr,max_yr,1): # Do this for every year (starting 1800)
        i=yr-min_yr
        logs("Year: {}".format(yr))

        S = simMat_yr[i].copy()
        P = symmetrize(S)    

        # Get elements that exist on this year
        mask = P.sum(axis=0)!=0

        P = P[mask][:,mask]   # Matrix with only elements that exist
        elems_i = np.array(elemList)[mask]  # List of existing elements

        # Evolve all populations in parallel for the given year
        with mp.Pool(processes=size) as pool:
            results = [pool.apply_async(optPop,
                                        args=(popSize,NGens,i,elems_i.shape[0],P) )
                       for i in range(size)]
            results_get = [r.get() for r in results]   

        # Get and transform results
        bestIndivs_yr = []
        for pop in results_get:
            if elems_i.shape[0]<103:
                filled_genes = completeList(genToElem(pop.bestIndiv.gen,ref=elems_i))
            else: filled_genes = pop.bestIndiv.gen

            logs(pop.bestIndiv.cost)
            bestIndivs_yr.append(filled_genes)

        logs("\tTime:  {}".format(time()-t0))

        results_total[i+min_yr]=bestIndivs_yr
    
    # Save results

    filehandler = open(resultsPath + 'optim_permut_yearly.gen', 'wb') 
    pickle.dump(results_total, filehandler)
 

def optPop(popSize,NGens,i,lenElems,Pi):
    T = 0.7#np.random.random()*0.35 + 0.35 # Random T between 0.35 and 0.7
    mutRate = 0.3#np.random.random()*0.2 + 0.2 # Random mutRate btwn 0.2 and 0.4
    breed = Population(popSize,lenElems,seed=i,P=Pi).evolve(NGens,T=T,mutRate=mutRate,Round=0)
    for i in range(2):
        T = 0.7*T
        breed = breed.evolve(NGens,T=T,mutRate=mutRate,Round=i+1)
    return breed

def getElemList(dataPath):
    elemList = []
    with open("{}/ElementList.txt".format(dataPath),'r') as f:
        for line in f:
            elemList.append(line.strip())
    return elemList



def symmetrize(S):
    """Get symmetrized version (P) of matrix S"""
    # Set all nans to 0
    S = np.nan_to_num(S,0)    

    diag = np.diag(S)
    n = S.shape[0]
    Sum0 = S.sum(axis=0).reshape(-1,1).repeat(n,axis=1) - np.diag(diag)
    Sum1 = S.sum(axis=1).reshape(1,-1).repeat(n,axis=0) - np.diag(diag)

    P = np.sqrt(S**2/(Sum0*Sum1+1e-5))

    P = np.nan_to_num(P,0)  # If entry is nan, then element doesn't exist. Make weight == 0.

    # Set diagonal elements to 0
    inds = np.arange(0,n)
    P[inds,inds] = 0
    return P

def cost(P, order):
    i = order.reshape(1,-1)
    diffs = np.abs(i - i.T) + np.eye(i.shape[1])
    F = -np.sum(P/diffs)
    return F

class Individual:
    def __init__(self,N,gen=False,P=False):
        self.N = N
        self.P = P
        if type(gen)!=bool:            self.gen = gen
        else:              self.gen = np.random.permutation(self.N)
        self.cost = cost(self.P,self.gen)

    def PMX(self,parent):
        """Partially-mapped crossover (PMX)
        Assumes NParents=2
        parent is an `Individual` object containing the other parent"""
        NParents = 2
        
        offsp = -np.ones((NParents,self.N),dtype=np.short)
        cutpts = np.sort(np.random.randint(1,self.N-1,2))
        
        m0,m1 = self.gen[cutpts[0]:cutpts[1]] , parent.gen[cutpts[0]:cutpts[1]]
        # Define initial mapping
        offsp[0,cutpts[0]:cutpts[1]] = m1
        offsp[1,cutpts[0]:cutpts[1]] = m0
        
        ### Start filling
        offsp[0,:cutpts[0]] = self.gen[:cutpts[0]]
        offsp[1,:cutpts[0]] = parent.gen[:cutpts[0]]
        offsp[0,cutpts[1]:] = self.gen[cutpts[1]:]
        offsp[1,cutpts[1]:] = parent.gen[cutpts[1]:]
        
        map0 = dict(zip(m1,m0))
        map1 = dict(zip(m0,m1))
        
        for i in range(cutpts[0]):
            while offsp[0,i] in m1:
                offsp[0,i] = map0[offsp[0,i]]
            while offsp[1,i] in m0:
                offsp[1,i] = map1[offsp[1,i]]
              
        for i in range(cutpts[1],self.N):
            while offsp[0,i] in m1:
                offsp[0,i] = map0[offsp[0,i]]
            while offsp[1,i] in m0:
                offsp[1,i] = map1[offsp[1,i]]
                
        return [Individual(self.N,gen=offsp[i],P=self.P) for i in range(2)]
    
    def mutate(self):
        """Move a random slice of the genes to a random place"""
        
        indiv = self.gen
        cutpts = np.sort(np.random.randint(0,self.N,2))
        Slice = indiv[cutpts[0]:cutpts[1]].copy()
        
        Left = np.concatenate([indiv[:cutpts[0]],indiv[cutpts[1]:]])              
        move_to = np.random.randint(0,Left.shape[0])
        
        self.gen = np.concatenate([Left[:move_to],Slice,Left[move_to:]]) 
        self.cost = cost(self.P,self.gen) # Update cost
    
class Population:
    def __init__(self,size,N,NParents=2,gens = False,Indivs=False,bestIndiv=False,seed=0,P=False):
        self.size = size  # Number of individuals
        self.N = N        # Number of elements in dataset
        self.NParents = NParents # Number of parents for a single breed
        self.seed = seed
        self.P = P
        
        np.random.seed(int(time()+seed))
        if not Indivs:
            if gens:        self.indivs = [Individual(N,gen=gens[i],P=self.P) for i in range(size)]
            else:           self.indivs = [Individual(N,P=self.P) for i in range(size)]
        else:               self.indivs = Indivs
            
        self.costs = [i.cost for i in self.indivs]
        
        ### Store best global configuration
        if bestIndiv:     
            self.bestIndiv = bestIndiv
            self.bestCost = bestIndiv.cost
        else:
            self.bestIndiv = False
            self.bestCost = 0
        
    def minCost(self):
        return np.min(self.costs)
    def meanCost(self):
        return np.mean(self.costs)
    def getBestIndiv(self):
        return self.indivs[np.argmin(self.costs)]
    
    def partnerUp(self):
        """Create a group of individuals for crossover"""
        parents = []
        idx = np.random.choice(np.arange(self.size),size=self.NParents,replace=False,p=self.probs)
        for i in idx:
            parents.append(self.indivs[i])
        return parents
    
    def newBreed(self,crossover,mutRate):
        newIndivs = []
        for i in range(self.size//2):
            parents = self.partnerUp()  # Select couples
            
            # Do crossover on these parents
            if crossover == 'PMX':        indivs_i = parents[0].PMX(parents[1])
    
            # Mutations
            for i in range(self.NParents):
                # Mutate each child with probability mutRate
                if np.random.random() < mutRate:    indivs_i[i].mutate()
                
            newIndivs = newIndivs + indivs_i
        return newIndivs
            
    def evolve(self,NGeners,T=0.7,crossover='PMX',mutRate=0.1,Round=0):
        for i in range(NGeners):
            # Probabilites of breeding for each individual
            # Use Boltzmann's distrib. with T. A lower T means an more uneven distribution (~ 0.x)
            self.probs = np.array([np.exp(-cost/T) for cost in self.costs])
            Z = np.sum(self.probs)
            self.probs = self.probs/Z

            self = Population(self.size,self.N,bestIndiv=self.bestIndiv,seed=self.seed,
                              Indivs=self.newBreed(crossover,mutRate),P=self.P)
            
            minc,meanc = self.minCost(), self.meanCost()
            if minc < self.bestCost:
                self.bestIndiv = self.getBestIndiv()
                self.bestCost = self.bestIndiv.cost
                
            if i%50==0:
                logs("** Round {}, Iter No. {}, Mean = {:.3f}, Min = {:.3f}, P{}".format(Round,i,meanc,minc,self.seed))
        return self

# Utility to convert gen into element list
def genToElem(gen,ref=False):
    """ref: list of elements of reference. When all elements present, ref = elements by AN"""
    if type(ref)==bool:
        ref = elemList
    order = ['_' for i in range(len(ref))]
    for i,idx in enumerate(gen):
        order[idx] = ref[i]
    return order

def completeList(incomp):
    """incomp: the incomplete, optimized list of elements when some element does not exist
    Returns a `genes` list, 
    containing the information of the inputed list and filled with the other elems"""
    
    sample_perm = np.random.permutation(np.arange(len(incomp),103))
    new_list = ['_' for i in range(103)]
    
    j=0
    for i in range(103):
        if elemList[i] in incomp:  # If element exists
            new_list[i] = incomp.index(elemList[i])  # Use same position as in incomplete list
        else:
            new_list[i] = sample_perm[j] # Else, fill with random index
            j+=1
    return np.array(new_list)

def logs(text):
    with open(root+'log/genetic1D.log','a') as f:
        f.write(str(text)+'\n')


if __name__ == '__main__':
    global root
    import argparse
    parser = argparse.ArgumentParser(description="Compute similarity matrices.")

    parser.add_argument("--NumProc","-N", required=True, type=int,
                        help="Number of processors to use for computations.")
    parser.add_argument("--inp_dir","-i", required=False, type=str, default="./",
                        help="Root directory for all computations")
    args = parser.parse_args()
    size=args.NumProc
    root=args.inp_dir

    main()
