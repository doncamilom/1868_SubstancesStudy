#! /usr/bin/env python3

# Instructions: Run as
# ./MPIgetTabs_v2.py <datafile.tsv> NP MaxMemMP NC
# NP=Number of processors to run in parallel
# NC=Number of compunds to use from dataset
# MaxMemMP=Maximum size of the data being transfered among subprocesses, in MB. It should be about 700

import multiprocessing as mp
from scipy import sparse as sp
from time import time
from itertools import chain
import numpy as np
from makeVecs_sparse import allVecs_sparse
import TPs   #import periodic tables
import sys

def main():
    global TP, TPshape, templateTable, NMax, size, maxWeightArray

    DataFile = sys.argv[1]
    size = int(sys.argv[2])
    maxWeightArray = float(sys.argv[3])
    if len(sys.argv)<5:     NMax = None   # If no max file length is passed, default to all
    else:                   NMax = int(sys.argv[4])

    # Preprocess compound data (make cmpnd vecs, years and ID) + produce element list.
    cmpnds,years,subsID, FullElemntList , NMax = allVecs_sparse(DataFile,NMax) 



    #########
    # DEV ###  In real DS, isomers are cleaned already so there's no need to do np.unique on cmpnds

    cmpnds = sp.csr_matrix(np.unique(cmpnds.toarray(),axis=0))
    years , subsID = years[:cmpnds.shape[0]] , subsID[:cmpnds.shape[0]] 

    #########



    #Construct a dict to go from element symbol to an index
    elemDict = {}
    for i,elem in enumerate(FullElemntList):
        elemDict[i] = elem
    
    #### Load periodic table schemes for building representation
    TP = TPs.TP
    TPshape = getTPShape()

    avgN = np.mean(np.sum(cmpnds,axis=1))
    print(f"\n * Average num. of atoms per compound: {round(avgN,4)}")
    print(f" * Max number of fragments R: {avgN*cmpnds.shape[0]:.0f}")
    print(f"\n * {cmpnds.shape[0]} unique compounds out of {NMax} requested...")
    print(f"\t Which means we dropped {NMax-cmpnds.shape[0]} compounds for being non-stoichiometric.\n")

    t0=time()
    R_sparse = findRs(cmpnds)  # Returns a scipy.sparse.csr_matrix containing all possible (R,n)s.

    print(f" * All possible Rs were produced in: {time()-t0:.3f} s")

    t0 = time()
    print("Finding unique Rs...")
    Rs_uniq = unique_mp(R_sparse,size)
    t = time()-t0

    sz = sum([i.shape[0] for i in Rs_uniq])
    print(f" * Found a total of {sz} non-unique Rs in {t:.3f} s")

    print(f"\nStarting finding commonalities.")

    t0 = time()
    Matches = list(chain(*[get_matches(rs,cmpnds,years,subsID,elemDict) for rs in Rs_uniq]))

    print("Done!")
    print(f"Time: {time()-t0:.3f} s\t NProc: {size}\t NMax: {NMax} ")



    ##### Calculate stats
    print(len(Matches))
    lensR = [len(x[0]) for x in Matches]
    print(sum(lensR)/len(Matches))
        
    ########### DEV
    return Matches



    R_distrib = distributeRs(Rlist,size)
    CommList, Rlist = run_getCommonal_distrib([R_distrib,cmpnds,years,subsID,elemDict,True])



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

#######################
# Auxiliary functions #
#######################

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

def validRs(Rs):
    """Get all non-unique R-n vectors out of a chunk Rs."""
    new_rs , c = np.unique(Rs.toarray(),axis=0, return_counts=True)
    new_rs = new_rs[c > 1]
    return sp.csr_matrix(new_rs,dtype=np.short)

def unique_mp(Rs,size):
    """Create data chunks for finding non-unique Rs in parallel.
    Test how this does on full DS. If pickling errors, implement recursive spliting.
    """
    max_n = 4
    # First split by n, the most natural choice.
    dist_list = [Rs[(Rs[:,-1]==i).toarray()[:,-1] ]
                  for i in range(1,max_n)]   + [Rs[(Rs[:,-1]>=max_n).toarray()[:,-1]]]

    with mp.Pool(processes=size) as pool:
        R_results = [pool.apply_async(validRs,args=(r,))
                     for r in dist_list]        

        Rs_get = [r.get() for r in R_results]

    return Rs_get

def get_matches(Rs,cmpnds,years,subsID,elemDict):
    """Get matches. For each R(n) find all elements X such that compound R-Xn exists in dataset. 
    Build an element set for each R(n).
    """
    ns = Rs[:,-1].data
    R  = Rs[:,:-1]
    
    sumCmpnds = cmpnds.sum(axis=1)
    sumRaxis1 = np.array(    R.sum(axis=1).flatten() + ns    ).flatten()
    
    Matches = []
    for i,n in enumerate(ns):
        if i%1000==0:       print( f"\t{i}th R evaluated..." )
            
        r = R[i] #The actual R
        
        """Encode a condition to search only within a subset of compounds
        fulfulling certain conditions based on R"""
        # 1. R is contained in compound        
        cond1 = ((cmpnds - r.toarray())>=0).all(axis=1)
        # 2. sum of atoms in cmpnd == sum of atoms in R_ (sum(R) + n)
        cond2 = (sumCmpnds == sumRaxis1[i])
        
        cond = np.array(cond1 & cond2).flatten()  # Combine conditions
        subsetCmpnds = cmpnds[cond]  # Select subset of cmpnds
        curr_years = years[cond]
        curr_subsID = subsID[cond]
        
        cmpnds_no_R = (subsetCmpnds - r.toarray())
        
        # Now select only those cmpnds where residual is due to one element only (X_n)
        # Only useful for n!=1
        if n!=1:
            cond = np.array((cmpnds_no_R!=0).sum(axis=1)==1).flatten()
            subsetCmpnds = subsetCmpnds[cond]
            curr_years = curr_years[cond]
            curr_subsID = curr_subsID[cond]
            
        # At this point, subsetCmpnds contains all compounds that match with R(n).
        elemIndex = (subsetCmpnds - r.toarray()).nonzero()[1]
        curr_list = set(map(lambda x: elemDict[x]  , elemIndex))  # Map dict to above list of elems

        Matches.append( [curr_list, curr_years, curr_subsID] )

        ###########
        ## Deal with this after you get all data created above
        #Table_list.append(getTable(curr_list,curr_years,curr_subsID,useID=useID))      
        ###########
    return Matches      # [list_elems,list_years, list_ids] for each R(n)


###########################








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




def distribRs_forUnique(Rs,max_n):
    # Create data chunks so that resulting chunks weigh at most 700 MB so that pickling isn't a problem when using Pool.
 
    # Go through each of the created chunks and, if any is above 700 MB, further split it.

    # Define recursive function to further split chunks
    def split_chunk(chunk,i):
        maxSplit = 5  # Maximum number of subsplits you want 

        # Split using ith index:
        split = []
        step = 2
        for j in range(maxSplit):
            lower,upper = j*step,(j+1)*step 
            if j < maxSplit-1:         tmp_splt = chunk[(chunk[:,i]>=lower) & (chunk[:,i]<upper)]  # Entries that are either j or j+1
            else:                      tmp_splt = chunk[chunk[:,i]>=lower]   # Entries that are maxSplit-1 or larger
            split.append(tmp_splt)

        # Now recursively further split here
        newList = []
        for l in split:
            if l.nbytes/1e6 > maxWeightArray:
                print("\t** Found a very large chunk! (Inside recursive function)")
                newList = newList + split_chunk(l,i+1)  # Further split l by next i
            elif l.shape[0]>1:      newList.append(l)  # Append only if chunk contains more than one entry
        #    else:   newList.append(l)  # Append only if chunk contains more than one entry

        return newList


    # First split by n, the most natural choice.
    dis_list = [Rs[Rs[:,-1]==i] for i in range(1,max_n)]

    print("\nStarting recursive splitting of Rs")

    new_chunks = []
    for l in dis_list:
        if l.nbytes/1e6 > maxWeightArray:
            print("\t** Found a very large chunk!")
            new_chunks = new_chunks + split_chunk(l,0)
        else:        new_chunks.append(l)
    
    print("Ended recursion")

    return new_chunks

if __name__ == '__main__':
    t0 = time()
    main()
    print(f"\n\nTotal runtime: {time()-t0:.3f} s")
    

