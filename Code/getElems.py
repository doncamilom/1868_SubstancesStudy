#! /usr/bin/env python3

import pandas as pd
import re

#This program extracts the elements present on dataset data.csv

def getElems(DataFile,NMax):
    """NMax specifies number of rows of dataset to be used"""
    
    col_names = ['ID','year','formula']
    sep = '\t'

    if NMax == 0:  df = pd.read_csv(DataFile,header=None,sep=sep,names=col_names)  #Load all data
    else:          df = pd.read_csv(DataFile,header=None,sep=sep,nrows=NMax,names=col_names)  #Load data

    df['formula'] = df['formula'].str.strip()   #Remove white spaces at begginning and end of string 

    elems = set([])

    for cmpnd in df['formula']:
        txt = ''.join(re.findall(r'[A-z]',cmpnd))   #Remove all subindices (there must be a regex to this but who knows)
        elems = elems.union(  set(re.split(r"(?<!^)(?=[A-Z])",txt))  )  # Add elements of this set to the set of known elements

    elems = sorted(list(elems)) # Convert to list and sort

    return elems
