#! /usr/bin/env python3

import numpy as np
import makeVecs
import TPs   #import periodic tables
import sys

def main():
    global FullElemntList
    global TP
    global TPshape

    DataFile = sys.argv[1]
    cmpnds = makeVecs.allVecs(DataFile)  #List of vectors representing all compounds in dataset
    FullElemntList = makeVecs.getElemList(DataFile)

    #Construct a dict to go from element symbol to an index
    elemDict = {}
    for i,elem in enumerate(FullElemntList):
        elemDict[i] = elem
    
    avgN = np.mean(np.sum(cmpnds,axis=1))
    print("Average num. of atoms per compound: ", round(avgN,4))
    print("Max num. ligands R: ", avgN*cmpnds.shape[0])
    print()

    #### Load periodic table schemes for building representation (pick one)
    #TP = TPs.TP
    TP = TPs.TP
    TPshape = getTPShape()

    from time import time
    t0=time()

    print("Finding unique Rs...")
    Rlist = findRs(cmpnds)
    print(f"Done in {time()-t0:.3f} seconds!")

    print("\nFinding commmonalities...")
    CommList,Rlist = getCommonal(Rlist,cmpnds,elemDict)
    print(len(CommList))
#    plotExample(CommList[0],[TP],'H')
    
    print(f"This operation took {time()-t0:.3f} seconds!")

## Auxiliary functions:

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
    print("Unique Rs: ", Rs.shape[0])

    return Rs 
    
def getCommonal(Rs,cmpnds,elemDict):
    """Get Commonalities. For each R(n) find all elements X such that compound R-Xn exists in dataset. 
    Build a list of these for each R(n)"""
    
    Comm_list = []
    R_list = []
    Table_list = []
    # Now generate the R representations on TP (TPR).

    j=0 #Counter for Rs with more than one appearence
    for i,R_ in enumerate(Rs):
        n = R_[-1]  #Take subindex 
        R = R_[:-1] #The actual R
        curr_list = []  #List to hold elements X, where R-X exist, for this R
        
        for comp in cmpnds:  #Loop through all compounds
            elemInCmpnd = np.nonzero(comp)[0]
            for i in elemInCmpnd:  #loop through elements present in this compound
                tmp = np.zeros(comp.shape[0])        
                tmp[i] = n  
                if np.all((comp - tmp)==R): #If removing element i (in proportion n) from this compound gives the same R, then R-in exists 
                    curr_list.append(elemDict[i])  #Add this element to the list
        
        if len(curr_list)>1:  #Only consider Rs that exist in more than 1 compound
            #np.save(f'TPR1/R{j}.npy', getTable(TP,curr_list)) 
            Table_list.append(getTable(TP,curr_list)) 
            Comm_list.append(curr_list)
            R_list.append(R_)
            j+=1

    np.save('Data/R_Tables.npy',np.array(Table_list))
    np.save('Data/R_Vector.npy',np.array(R_list))
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


