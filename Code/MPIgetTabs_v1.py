#! /usr/bin/env python3

from mpi4py import MPI
from time import time
import numpy as np
import makeVecs
import TPs   #import periodic tables
import sys

def main():

    global FullElemntList, TP, TPshape
    global size,rank

    comm = MPI.COMM_WORLD
    rank=comm.Get_rank() #process id
    size=comm.Get_size() #number of procesors

    DataFile = sys.argv[1]
    NMax = int(sys.argv[2])

    cmpnds = makeVecs.allVecs(DataFile,NMax,rank)  #List of vectors representing all compounds in dataset

    FullElemntList = makeVecs.getElemList(DataFile,NMax)

    #Construct a dict to go from element symbol to an index
    elemDict = {}
    for i,elem in enumerate(FullElemntList):
        elemDict[i] = elem
    
    #### Load periodic table schemes for building representation (pick one)
    TP = TPs.TP
    TPshape = getTPShape()

    if rank==0:
        avgN = np.mean(np.sum(cmpnds,axis=1))
        print("\n * Average num. of atoms per compound: ", round(avgN,4))
        print(" * Max num. ligands R: ", avgN*cmpnds.shape[0])
        print()

        print("Finding unique Rs...")

    t0=time()
    Rlist = findRs(cmpnds)


    if rank==0: 
        print(f"\n * {cmpnds.shape[0]} unique compounds out of {NMax} requested.\n")
        print(" Unique Rs: ", Rlist.shape[0])

    #################
    ## Following code taken from https://gist.github.com/krischer/2c7b95beed642248487a
    def split(container, count):
        return [container[_i::count] for _i in range(count)]

    if rank == 0:    Rlist = split(Rlist, size)
    else:            Rlist = None

    Rlist = comm.scatter(Rlist, root=0)
    ## Up to here
    ################


    print(f"  -- P{rank} Finding commonalities...\n")
    CommList,Rlist = getCommonal(Rlist,cmpnds,elemDict)
    if rank==0:
        print(f"Time: {time()-t0:.3f} s\t NProc: {size}\t NMax: {NMax} ")

    
    # Collect data and merge
    CommList_ = MPI.COMM_WORLD.gather(CommList,root=0)    

    if rank==0:
        CommListTotal = []
        for i in range(size):
            CommListTotal = CommListTotal+CommList_[i]
        
        #### Let's calculate some useful values
        # Mean number of element per list?
        LenR = [ len(set(x)) for x in CommListTotal ] 
        meanLen = sum(LenR)/len(LenR)
        print("Mean number of compounds per R :",meanLen)
        print("Max number of compounds per R :",max(LenR))


#        plotExample(CommList[np.argmax(LenR)],[TP],'H')



def getTPShape():
    """Return the shape of the global PT"""
    XY = np.array(list(TP.values()))  #Positions of elements in PT
    return tuple(np.max(XY,axis=0) + 1)
    
# Function to construct a TPR given a TP layout and a list of elements.
def getTable(TP,listElem,incogn=False):
    """This gets as input: a periodic system (table) and a list of elements. 
        incogn. is the element I'll be asking for (does compound R-incogn exist?)
        Returns a 2D array coloured at the positions of these elements
        Shape is (#cols, #rows) #cols = 6, #rows = 32 if Lanthanides included, 18 if not"""
    table = -np.ones(TPshape)
    for elem in FullElemntList:   # ONLY EXISTENT ELEMENTS ARE SHOWN IN COLOR 
        x,y = TP[elem]
        if elem in listElem:   table[x,y] = 1.  #If in element list: one color
        elif elem==incogn:     table[x,y] = -2
        else:                  table[x,y] = 0.  #If not, other.
    return table

def findRs(cmpnds):
    # Find all unique Rs
    n = len(cmpnds[0])
    Rs = []
    ns = []
    
    for c in cmpnds:
        indx = np.nonzero(c)  #Get index of non-zero entries
        for i in indx[0]:     #Loop through these
            c_ = c.copy()
            n = int(c_[i].copy())
            for j in range(n):  #Loop through ith element's subindex
                c_[i] -= 1      #Remove one
                Rs.append(c_.copy())   #Append the compound with a reduced entry (R-Xn-1)
                ns.append(j+1)  #How many atoms of this element have been removed 
    
    Rs = np.array(Rs)
    ns = np.array(ns).reshape(-1,1)
    
    Rs = np.concatenate([Rs,ns],axis=1)   # First columns represent an R ligand. last column represents subind. n
    Rs = np.unique(Rs,axis=0) #This array contains all unique Rs.

    return Rs
 
    
def getCommonal(Rs,cmpnds,elemDict):
    """Get Commonalities. For each R(n) find all elements X such that compound R-Xn exists in dataset. 
    Build a list of these for each R(n)
    rank is the number of processors in which the operation is to be run"""
    
    Comm_list = []
    R_list = []
    Table_list = []
    # Now generate the R representations on TP (TPR).

    j=0 #Counter for Rs with more than one appearence
    sumCmpnds = cmpnds.sum(axis=1)

    for i,R_ in enumerate(Rs):
        if i%1000==0:       print( f"\tPast {i}th R" )

        n = R_[-1]  #Take subindex 
        R = R_[:-1] #The actual R
        curr_list = []
        
        # Encode a condition to search only within a subset of compounds fulfulling certain conditions based on R
        # 1. R is contained in compound
        cond1 = ((cmpnds - R) >= 0).all(axis=1)        
        # 2. sum of atoms in cmpnd == sum of atoms in R_ (sum(R) + n)
        cond2 = (cmpnds.sum(axis=1) == R_.sum())

        subsetCmpnds = cmpnds[cond1 & cond2]  # Select subset of cmpnds

        cmpnds_no_R = (subsetCmpnds - R)
        # Now select only those cmpnds where residual is due to one element only (X_n)
        subsetCmpnds = subsetCmpnds[(cmpnds_no_R!=0).sum(axis=1)==1]

        if subsetCmpnds.shape[0] > 1:

            curr_list = list(map(lambda x: elemDict[x]  , (subsetCmpnds - R).nonzero()[1] ))
            #np.save(f'TPR1/R{j}.npy', getTable(TP,curr_list)) 
            Table_list.append(getTable(TP,curr_list)) 
            Comm_list.append(curr_list)
            R_list.append(R_)
            j+=1


    print(len(R_list))
    print("Saving...")
    np.save(f'../Data/R_Tables_P{rank}_new.npy',np.array(Table_list))
    np.save(f'../Data/R_Vector_P{rank}_new.npy',np.array(R_list))
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


