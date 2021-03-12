#! /usr/bin/env python3

import numpy as np
import pandas as pd
from itertools import chain
from scipy import sparse as sp
import re

#This function extracts the elements present on dataset data.csv
def getElems(DataFile,NMax=None):
    """NMax specifies number of rows of dataset to be used"""
    col_names = ['ID','formula','year']
    sep = '\t'

    df = pd.read_csv(DataFile,header=None,sep=sep,nrows=NMax,names=col_names)  #Load data
    df['formula'] = df['formula'].str.strip()   #Remove white spaces at begginning and end of string 

    elems = set()
    for cmpnd in df['formula']:
        txt = ''.join(re.findall(r'[A-z]',cmpnd))   #Remove all subindices (there must be a regex to this but who knows)
        elems = elems.union(  set(re.split(r"(?<!^)(?=[A-Z])",txt))  )  # Add elements of this set to the set of known elements
    elems = sorted(list(elems)) # Convert to list and sort

    # Save this list of elements so it doesn't have to be calculated every time
    with open("./Data/ElementList.txt", "w") as f:
        for A in elems:
            f.write(str(A) +"\n")
    return elems  # This returns a list with all sorted elements in dataset

def getVec_sparse(tx,elemList): #Tx = molecular formula, e.g. C4H2ClBr
    #### This regex handles non-integer subindices: C6H16Na3O12.5 (which happens in DS) 
    Li = re.split(r"(?<!^)(?=[A-Z])",tx)  #Split as ['H2','O']
    
    # Adds 1 if no subindex. Result is ['H2','O1']. 
    # Right after, split chem symbol from subindex as [['H',2],['O',1]]
    
    li = [re.split(r"([A-z]+)(([0-9]*[.])?[0-9]+)",i)
          if bool(re.match(r'[A-z]*([0-9]*[.])?[0-9]+',i))
          else re.split(r"([A-z]+)(([0-9]*[.])?[0-9]+)",i+'1') for i in Li]  
    
    # Construct two lists: input for sparse matrix construction
    col  = [elemList.index(i[1]) for i in li]  # Index of element i to put correspondent data
    data = [float(i[2]) for i in li]           # Num. atoms of element i
    
    for i in data:
        if float(i)!=int(i):
            return None  # Return empty lists or better None?
    return col,data


def allVecs_sparse(DataFile,NMax=None):
    col_names = ['ID','formula','year']
    sep = '\t'

    df = pd.read_csv(DataFile,header=None,sep=sep,nrows=NMax,names=col_names)  #Load data
    NMax = df.shape[0]

    df['formula'] = df['formula'].str.strip()   #Remove white spaces at begginning and end of string 

    elemList = getElems(DataFile,NMax)
    
    # List of lists [col,data]
    colXdata = list(map(lambda x: getVec_sparse(x,elemList) , df['formula'].values))
    index = [i for i, l in enumerate(colXdata) if l is not None]
    colXdata = [l for l in colXdata if l is not None]
    # DB is assumed to be already cleaned for repeated entries. That is, we already got rid of isomers.
    
    # See docs for scipy.sparse.csr_matrix to understand the syntaxis
    indptr = np.cumsum([0]+list(map(lambda x: len(x[0]) , colXdata)))
    indices = np.array(list(chain(*[l[0] for l in colXdata])))
    data = np.array(list(chain(*[l[1] for l in colXdata])))

    cmpnds = sp.csr_matrix((data, indices, indptr), 
                           shape=(len(colXdata), len(elemList)),
                           dtype=np.short)
       
    years = df['year'].values[index]
    subsID = df['ID'].values[index]

    return cmpnds,years,subsID, elemList, NMax
