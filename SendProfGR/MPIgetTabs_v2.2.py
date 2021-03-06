#! /usr/bin/env python3

# Instructions: Run as
# ./MPIgetTabs_v2.py <datafile.tsv> NP NC
# With NP=Number of processors to run in parallel
# NC=Number of compunds to use from dataset

import multiprocessing as mp
from time import time
import numpy as np
import makeVecs
import TPs   #import periodic tables
import sys

def main():
    global TP, TPshape, templateTable, NMax, size

    DataFile = sys.argv[1]
    size = int(sys.argv[2])
    if len(sys.argv)<4:     NMax = 0   # If no max file length is passed, default to all
    else:                   NMax = int(sys.argv[3])

    # Preprocess compound data (make cmpnd vecs, years and ID) + produce element list.
    cmpnds,years,subsID, FullElemntList  = makeVecs.allVecs(DataFile,NMax) 


    #Construct a dict to go from element symbol to an index
    elemDict = {}
    for i,elem in enumerate(FullElemntList):
        elemDict[i] = elem
    
    #### Load periodic table schemes for building representation
    TP = TPs.TP
    TPshape = getTPShape()
    
    #### Make a template table so this doesn't have to be computed over and over
    templateTable = -np.ones(TPshape,dtype=np.intc)
        # Map -1 for every element that is in DS 
    indexTempTab = np.array(list(map(lambda x: TP[x], FullElemntList)))
    templateTable[indexTempTab[:,0],indexTempTab[:,1]] = -1

    avgN = np.mean(np.sum(cmpnds,axis=1))
    print(f"\n * Average num. of atoms per compound: {round(avgN,4)}")
    print(f" * Max number of fragments R: {avgN*cmpnds.shape[0]:.0f}")
    print()

    print("Finding unique Rs...")

    t0=time()
    Rlist = findRs(cmpnds)

    print(f"\n * {cmpnds.shape[0]} unique compounds out of {NMax} requested.\n")
    print(" Total unique (usable) Rs: ", Rlist.shape[0])

    print(f"\nStarting finding commonalities. Time so far: {time()-t0}")


    R_distrib = distributeRs(Rlist,size)
    CommList, Rlist = run_getCommonal_distrib([R_distrib,cmpnds,years,subsID,elemDict,True])

    print("Done!")
    print(f"Time: {time()-t0:.3f} s\t NProc: {size}\t NMax: {NMax} ")




    #### Let's calculate some useful values
    CommListTotal = []
    for i in range(size):
        CommListTotal = CommListTotal+CommList[i]

    # Mean number of element per list?
    LenR = [ len(set(x)) for x in CommListTotal ] 
    meanLen = sum(LenR)/len(LenR)
    print("\tMean number of compounds per R :",meanLen)
    print("\tMax number of compounds per R (first 3 max) :",sorted(set(LenR))[-3:])

    # End of main program


def run_getCommonal_distrib(args):
    # Run getCommonal in parallel for the distributed array
    with mp.Pool(processes=size) as pool:
        # starts the sub-processes without blocking
        # pass the chunk to each worker process
        R_results = [pool.apply_async(getCommonal,
                                      args=(Rlist_i,rank,*args[1:],))
                     for rank,Rlist_i in enumerate(args[0])]

        # blocks until all results are fetched
        R_results_get = [r.get() for r in R_results]   # A list [  [CommList0, Rlist0] , [CommList1, Rlist1]  ]

    # Recover all results into a couple lists:   CommList,Rlist
    CommList, Rlist = [],[]
    for proc in range(size):
        CommList.append(R_results_get[proc][0])
        Rlist.append(R_results_get[proc][1])

    return CommList, Rlist

def distributeRs(Rlist,size):
    # Split R list for parallel processing
    R_distrib = []
    chunk_sz = len(Rlist)//size
    for proc in range(size):
        start = proc * chunk_sz
        end = (proc + 1) * chunk_sz if proc < size - 1 else None
        
        R_distrib.append(Rlist[start:end])
    return R_distrib

def getTPShape():
    """Return the shape of the global PT"""
    XY = np.array(list(TP.values()))  #Positions of elements in PT
    return tuple(np.max(XY,axis=0) + 1)
    
# Function to construct a TPR given a TP layout and a list of elements.
def getTable(listElem,years,subsID,useID=False):
    """This gets as input: a periodic system (table) and a list of elements. 
        Returns a 2D array coloured at the positions of these elements"""
    table = templateTable.copy()
    
    # What values to put in array. Either substance ID or year
    if useID: dataArray = subsID
    else:     dataArray = years
 
    for i,elem in enumerate(listElem):
        x,y = TP[elem]
        table[x,y] = dataArray[i]
    return table

def getRepeated(Rs):
    #global new_rs   # Making it global so it avoids pickling crashes
    new_rs , c = np.unique(Rs,axis=0, return_counts=True)
    new_rs = new_rs[c > 1]
    return new_rs
    
def findRs(cmpnds):
    # Find all unique Rs
    n = len(cmpnds[0])
    Rs = []
    
    cShape = cmpnds.shape[1]
    for c in cmpnds:
        indx = np.nonzero(c)  #Get index of non-zero entries
        for i in indx[0]:     #Loop through these
            c_ = np.zeros(cShape+1,dtype=np.short)
            c_[:-1] = c
            n = int(c_[i])
            for j in range(n):  #Loop through ith element's subindex
                c_[i] -= 1      #Remove one
                c_[-1] = j+1    #How many atoms of this element have been removed 

                Rs.append(c_.copy())  #Append the compound with a reduced entry (R-Xn-1)

    # Get only the Rs that were produced more than once (meaning it's guaranteed at least 2 elems per R)

    # Split Rs so np.unique runs in parallel + faster as new subprocesses are much lighter
    # Split by n (R[-1])

    Rs = np.array(Rs)
    max_n = 40               # Choose wisely, this may bias results a little (the bigger the better). Read next comment
    Rs = Rs[Rs[:,-1]<max_n]  # Cut Rs by the n. If n>max_n it's very (very) likely no two compounds share it

    global Rs_distrib_list   # Making it global so it avoids pickling crashes
    # Create data chunks
    Rs_distrib_list = [Rs[Rs[:,-1]==i] for i in range(1,max_n)]
    
    for i in Rs_distrib_list:
        print(i.shape)

    with mp.Pool(processes=size) as pool:
        # starts the sub-processes without blocking
        # pass the chunk to each worker process

        R_results = [pool.apply_async(getRepeated,
                                      args=(Rs_i,))
                     for Rs_i in Rs_distrib_list]

        # blocks until all results are fetched
        R_results_get = [r.get() for r in R_results]   # A list of arrays

    # Merge filtered Rs by concatenating the resulting list
    Rs = np.concatenate(R_results_get,axis=0)

    return Rs
 
def getCommonal(Rs,rank,cmpnds,years,subsID,elemDict,useID=False):
    """Get Commonalities. For each R(n) find all elements X such that compound R-Xn exists in dataset. 
    Build a list of these for each R(n)
    rank is the number of processors in which the operation is to be run"""
    
    if rank==0:   print("\n Finding commonalities (finding sets for each (R,n) pair)...\n")
    Comm_list = []
    R_list = []
    Table_list = []
    # Now generate the R representations on TP (TPR).

    j=0 #Counter for Rs with more than one appearence
    sumCmpnds = cmpnds.sum(axis=1)

    for i,R_ in enumerate(Rs):
        if i%1000==0 and rank==0:       print( f"\t{i}th R evaluated..." )

        n = R_[-1]  #Take subindex 
        R = R_[:-1] #The actual R
        curr_list = []
        
        # Encode a condition to search only within a subset of compounds fulfulling certain conditions based on R
        # 1. R is contained in compound
        cond1 = ((cmpnds - R) >= 0).all(axis=1)        
        # 2. sum of atoms in cmpnd == sum of atoms in R_ (sum(R) + n)
        cond2 = (cmpnds.sum(axis=1) == R_.sum())

        subsetCmpnds = cmpnds[cond1 & cond2]  # Select subset of cmpnds
        curr_years = years[cond1 & cond2]
        curr_subsID = subsID[cond1 & cond2]

        cmpnds_no_R = (subsetCmpnds - R)
        # Now select only those cmpnds where residual is due to one element only (X_n)
        subsetCmpnds = subsetCmpnds[(cmpnds_no_R!=0).sum(axis=1)==1]
        curr_years = curr_years[(cmpnds_no_R!=0).sum(axis=1)==1]
        curr_subsID = curr_subsID[(cmpnds_no_R!=0).sum(axis=1)==1]

        if subsetCmpnds.shape[0] > 1:

            curr_list = list(map(lambda x: elemDict[x]  , (subsetCmpnds - R).nonzero()[1] ))
            Table_list.append(getTable(curr_list,curr_years,curr_subsID,useID=useID)) 
            Comm_list.append(curr_list)
            R_list.append(R_)
            j+=1
    
    if rank==0:    print("Saving...")
    table_list_arr = np.array(Table_list,dtype=np.intc)
    np.save(f'./Data/TablesID_NMax{NMax}_P{rank}.npy',table_list_arr)
    np.save(f'./Data/RVector_NMax{NMax}_P{rank}.npy',np.array(R_list,dtype=np.short))


    ### Make a function here that produces a new array exchanging ID for year, just a mapping
    if useID:
        if rank==0:  print("\n\tProducing array filled with years from ID array...")
        
        # Mapping going from ID to year
        mapping = dict(zip(list(subsID),list(years)))
        mapping[-1] = -1
        # Apply mapping 
        yearArray = np.vectorize(mapping.__getitem__)(table_list_arr)
        np.save(f'./Data/TablesYears_NMax{NMax}_P{rank}.npy',np.array(yearArray,dtype=np.short))
            
    return Comm_list,R_list  #This list contains lists (one for each R) of elements X such that R-X exist in dataset.



def plotExample(listElem,TPs=[],IncognElem='H'):
    """ This code generates an example image of the representation we're looking for
    Expects as inputs:
        TP: A list of periodic tables
        IncognElem: An element to ask for prediction. e.g. does R-(this element) exist?"""

    import matplotlib.pyplot as plt
    
    if len(TPs)>1:
        fig,ax = plt.subplots(len(TPs),1)
        for i,TPi in enumerate(TPs):
            arr = getTable(TPi,listElem,IncognElem)
            ax[i].imshow(arr)
            ax[i].axis('off')
    else:
        fig,ax = plt.subplots()
        arr = getTable(TPs[0],listElem,IncognElem)
        ax.imshow(arr)
        ax.axis('off')
    plt.show()


if __name__ == '__main__':
    main()


