#! /usr/bin/env python3

import numpy as np
import pandas as pd
import re

#This function extracts the elements present on dataset data.csv
def getElems(DataFile,NMax):
    """NMax specifies number of rows of dataset to be used"""
    
    col_names = ['ID','formula','year']
    sep = '\t'

    if NMax == 0:  df = pd.read_csv(DataFile,header=None,sep=sep,names=col_names)  #Load all data
    else:          df = pd.read_csv(DataFile,header=None,sep=sep,nrows=NMax,names=col_names)  #Load data

    df['formula'] = df['formula'].str.strip()   #Remove white spaces at begginning and end of string 

    elems = set([])

    for cmpnd in df['formula']:
        txt = ''.join(re.findall(r'[A-z]',cmpnd))   #Remove all subindices (there must be a regex to this but who knows)
        elems = elems.union(  set(re.split(r"(?<!^)(?=[A-Z])",txt))  )  # Add elements of this set to the set of known elements

    elems = sorted(list(elems)) # Convert to list and sort

    # Save this list of elements so it doesn't have to be calculated every time
    with open("./Data/ElementList.txt", "w") as f:
        for A in elems:
            f.write(str(A) +"\n")

    return elems  # This returns a list with all sorted elements in dataset


def getVec(tx,elemList): #Tx = molecular formula, e.g. C4H2ClBr

    #### This regex handles non-integer subindices: C6H16Na3O12.5 (which happens in DS) 
    Li = re.split(r"(?<!^)(?=[A-Z])",tx)  #Split as ['H2','O']
    li = [i if bool(re.match(r'[A-z]*([0-9]*[.])?[0-9]+',i)) else i+'1' for i in Li]  #Adds 1 if no subindex. Result is ['H2','O1']
    
    # Construct vector where each entry corresponds to the num. of atoms of that element in the compound.
    vec = np.zeros(len(elemList))
    for i in li: #Loop through elements in compound
        #i is a string: Xn, where X is an element and n its subindex on the compound
        a = re.split(r"([A-z]+)(([0-9]*[.])?[0-9]+)",i) #Split these two components (X, n)
        elem = a[1] #Get element
        mul = float(a[2])  #Get subindex
        vec[elemList.index(elem)] += mul #Assign this index to the position of the element in elemList

    if np.all(vec == vec.astype(int)):     return vec
#    elif np.all(2*vec == (2*vec).astype(int)): 
#        #print(f" * Duplicating indices of *some substance*")#{tx}")
#        return 2*vec  #If some n==0.5 then an integer composition is *2.
    else:
        return 0*vec   #Problem can't be easily solved

def allVecs(DataFile,NMax):
    col_names = ['ID','formula','year']
    sep = '\t'

    if NMax == 0:  df = pd.read_csv(DataFile,header=None,sep=sep,names=col_names)  #Load all data
    else:          df = pd.read_csv(DataFile,header=None,sep=sep,nrows=NMax,names=col_names)  #Load data

    df['formula'] = df['formula'].str.strip()   #Remove white spaces at begginning and end of string 
    
    elemList = getElems(DataFile,NMax)     
    totalVecs = np.array(list(map(lambda x: getVec(x,elemList),df['formula'].values)),dtype=np.short)
    totalVecs,index = np.unique(totalVecs,axis=0,return_index=True)
    years = df['year'].values[index]    
    subsID = df['ID'].values[index]    

    return totalVecs,years,subsID, elemList, df.shape[0]

