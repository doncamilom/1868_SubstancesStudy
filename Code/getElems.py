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
    
    text = ''.join(df['formula'])  #Join all compounds into one single string
    
    text2 = ''.join(re.findall(r'[A-z]',text))  #Remove all subindices (there must be a regex to this but who knows)
    
    elems = sorted(list(set(re.split(r"(?<!^)(?=[A-Z])",text2))))  #All unique elements found in dataset
    return elems
