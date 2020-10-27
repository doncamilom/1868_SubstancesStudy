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

    Li = re.split(r"(?<!^)(?=[A-Z])",tx)  #Split as ['H2','O']
    li = [i if bool(re.match(r'[A-z]*\d+$',i)) else i+'1' for i in Li]  #Adds 1 if no subindex, else nothing. Result is ['H2','O1']
    
    # Construct vector where each entry corresponds to the num. of atoms of that element in the compound.
    vec = np.zeros(len(elemList))
    for i in li: #Loop through elements in compound
        #i is a string: Xn, where X is an element and n its subindex on the compound
        a = re.split(r"([A-z]+)(\d+)",i) #Split these two components (X, n)
        elem = a[1] #Get element
        mul = int(a[2])  #Get subindex
        vec[elemList.index(elem)] += mul #Assign this index to the position of the element in elemList
    return vec
    
def allVecs(DataFile,NMax):
    df = pd.read_csv(DataFile,header=None,nrows=NMax,names=['year','form'])  #Load data
    df['form'] = df['form'].str.strip()   #Remove white spaces at begginning and end of string 
    
    elemList = getElemList(DataFile,NMax)     
    return np.array(list(map(lambda x: getVec(x,elemList),df['form'].values)))


