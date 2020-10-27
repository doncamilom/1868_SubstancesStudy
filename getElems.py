#! /usr/bin/env python3

import pandas as pd
import re

#This program extracts the elements present on dataset data.csv

def getElems(DataFile,NMax):
    """NMax specifies number of rows of dataset to be used"""

    df = pd.read_csv(DataFile,header=None,nrows=NMax,names=['year','form'])  #Load data
    df['form'] = df['form'].str.strip()   #Remove white spaces at begginning and end of string 
    
    text = ''.join(df['form'])  #Join all compounds into one single string
    
    text2 = ''.join(re.findall(r'[A-z]',text))  #Remove all subindices (there must be a regex to this but who knows)
    
    elems = sorted(list(set(re.split(r"(?<!^)(?=[A-Z])",text2))))  #All unique elements found in dataset
    return elems
