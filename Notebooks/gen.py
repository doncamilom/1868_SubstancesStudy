#! /usr/bin/env python3

import numpy as np

class Individual:
    def __init__(self,N,gen=False):
        self.N = N
        if type(gen)!=bool:            self.gen = gen
        else:              self.gen = np.random.permutation(self.N)
        self.cost = cost(P,self.gen)

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
                
        return [Individual(self.N,gen=offsp[i]) for i in range(2)]
    
    def mutate(self):
        """Move a random slice of the genes to a random place"""
        
        indiv = self.gen
        cutpts = np.sort(np.random.randint(0,self.N,2))
        Slice = indiv[cutpts[0]:cutpts[1]].copy()
        
        Left = np.concatenate([indiv[:cutpts[0]],indiv[cutpts[1]:]])              
        move_to = np.random.randint(0,Left.shape[0])
        
        self.gen = np.concatenate([Left[:move_to],Slice,Left[move_to:]]) 
        self.cost = cost(P,self.gen) # Update cost
    
class Population:
    def __init__(self,size,N,NParents=2,T=0.7,gens = False,Indivs=False):
        self.size = size  # Number of individuals
        self.N = N        # Number of elements in dataset
        self.NParents = NParents # Number of parents for a single breed
        self.T = T        # Temperature for probabilities.
        
        if not Indivs:
            if gens:        self.indivs = [Individual(N,gen=gens[i]) for i in range(size)]
            else:           self.indivs = [Individual(N) for i in range(size)]
        else:               self.indivs = Indivs
            
        self.costs = [i.cost for i in self.indivs]
        
        # Probabilites of breeding for each individual
        # Use Boltzmann's distrib. with T. A lower T means an more uneven distribution (~ 0.x)
        self.probs = np.array([np.exp(-cost/T) for cost in self.costs])
        Z = np.sum(self.probs)
        self.probs = self.probs/Z
        
        ### Store best global configuration
        self.bestCost = 0
        self.bestIndiv = False
        
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
            
    def evolve(self,NGeners,crossover='PMX',mutRate=0.1):
        Mins,Means = [], []
        for i in range(NGeners):
            self = Population(self.size,self.N,T=self.T,
                              Indivs=self.newBreed(crossover,mutRate))
            
            minc,meanc = self.minCost(), self.meanCost()
            if minc < self.bestCost:
                self.bestCost = minc
                self.bestIndiv = self.getBestIndiv()
                
            if i%10==0:
                Mins.append(minc)
                Means.append(meanc)
                print(f"Iter No. {i}, Mean = {meanc:.3f}, Min = {minc:.3f}",end='\r')
        plt.plot(Mins)
        plt.plot(Means)
        return self
