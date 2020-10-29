#! /usr/bin/env python3

import numpy as np
import pandas as pd
import re
import getElems


#This program extracts the elements present on dataset data.csv
def getElemList(DataFile,NMax): 
    """NMax specifies how many rows of the dataset are to be considered"""

    elemList = getElems.getElems(DataFile,NMax) #Run getElems module. This returns a list with all sorted elements in dataset
    return elemList  #Run getElems module. This returns a list with all sorted elements in dataset


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
    elif np.all(2*vec == (2*vec).astype(int)): return 2*vec  #If some n==0.5 then an integer composition is *2.
    else: 
        print(  f" *** Trouble {tx}" )
        return 0*vec   #Problem can't be easily solved

def allVecs(DataFile,NMax):
    df = pd.read_csv(DataFile,header=None,nrows=NMax,names=['year','form'])  #Load data
    df['form'] = df['form'].str.strip()   #Remove white spaces at begginning and end of string 
    
    elemList = getElemList(DataFile,NMax)     
    totalVecs = np.array(list(map(lambda x: getVec(x,elemList),df['form'].values)))
    return np.unique(totalVecs,axis=0)


