#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set() ; sns.set_style('white')
import pandas as pd
from itertools import chain
from TPs import TP

import re
import bz2
import pickle
from scipy import sparse as sp
from time import time

def getElemList(dataPath):
    elemList = []
    with open(f"{dataPath}/ElementList.txt",'r') as f:
        for line in f:
            elemList.append(line.strip())
    return elemList

dataPath = '../Data/'
elemList = getElemList(dataPath)

elemList_AO = [] # Atomic weight ordered element list
for e in TP.keys():
    if e in elemList:
        elemList_AO.append(e)
        
def main():
    global elemList
    
#    dataPath = "./Data/"
    dataPath = "../Preprocess/Data/"

    # Load element list
    elemList = getElemList(dataPath)
    print(f"Number of elements: {len(elemList)}")

    t0 = time()
    sim_mat = calc_simMats_yearly(dataPath)
    print(time()-t0)
    np.save(dataPath+'history_simMat.npy',sim_mat)


def getTable(TP, elems, fill_data,emptyTab):
    """This gets as input: a PT and a list of elements. 
        Returns a 2D array coloured at the positions of these elements"""
    for elem,fill in zip(elems,fill_data):
        x,y = TP[elem]
        emptyTab[x,y] = fill
    return emptyTab

def getSimilarities_yr(Tyr,element,mask,year):
    """Get array of the similarities between elements and the given element
    based on 'replaceability' in chemical formulas, for the CS existent at a given year."""
    # Filter out by year
    T = (Tyr <= year)&(Tyr>0)
    
    # When filtering by year, we need only sum over those tables where at least 1 element is present!!!
    T = T[T.sum(axis=(1,2))>1]

    # Now start actually calculating simMats
    X,Y = TP[element] # Get coords of elem in PT
    # similarity: number of times other elements appear in same table
    simil = T[T[:,X,Y]].sum(axis=0)*1.
    
    # Apply mask: convert any "off table" entry into nan
    with np.errstate(invalid='ignore',divide='ignore'):
        simil /= mask
    # Extract similarity data for every element
    data_arr = {e:simil[TP[e][0],TP[e][1]] for e in TP.keys() if e in elemList}

    data_arr = pd.Series(data_arr).rename(element)
    return data_arr

def simMat(Tyr,mask,year):
    mat = pd.DataFrame(columns=[e for e in TP.keys() if e in elemList])
    for e in TP.keys():
        if e in elemList: # Calculate in order, but subjected to what's in elemList
            mat[e] = getSimilarities_yr(Tyr,e,mask,year)
    return mat  # Better return array, not dataframe

def calc_simMats_yearly(dataPath = "../OnlyStrsCode/Data/"):
    """Calculate all simMats yearly.
    First load all the data that was produced during preprocessing.
    Takes as input: dataPath, and number of elements in dataset"""

    # Load matches_pickle.bin
    Match_file = bz2.BZ2File(dataPath+'matches_pickle.bin', 'r')
    Matches = pickle.load(Match_file)
    miny,maxy = getMinMax(Matches)
    print(miny,maxy)

    Tyr = np.zeros((len(Matches),7,32))      # This gets the arrays filled with years
    for i,match in enumerate(Matches):
        getTable(TP, match[0], match[1], Tyr[i] )  # Match[1] to fill with years


    # mask: 7x32 bool: If at x,y an element exists, mask[x,y]=1. Else, mask[x,y]=0
    mask = (Tyr>0).sum(axis=0)>0

    # Generate data yearly
    year_span = maxy-miny+1
    numElems = len(elemList)
    simMats_yr = np.zeros((year_span,numElems,numElems))
    
    for i,yr in enumerate(range(miny,maxy+1)):
        simMats_yr[i] = simMat(Tyr,mask,yr)

    return simMats_yr


def plot_SimPTBar(simMat_yr,year,element,min_yr,save=False):
    # Select simMat for this year
    arr_yr = simMat_yr[year-min_yr].copy()
    # Select a particular element
    X,Y = TP[element]

    # Generate a list of elements present at the given year
    c,elems_yr = 0,[]  # counter, element list
    for e in TP.keys():
        if e in elemList:
            # If all entries at this place are nan: they don't exist
            if (~np.isnan(arr_yr[:,c].all())): 
                c+=1
                elems_yr.append(e)

    # Select the array for the given element, for the given year
    arr_thisElem = arr_yr[:,elems_yr.index(element)]
    arr_thisElem[elems_yr.index(element)] = 0  # Remove this element's value, so it's white as well

    img = np.zeros((7,32))
    mask = img.copy()
    # Create a mask to wipe out nan entries, so they appear as white
    c = 0
    for e in elems_yr:
        x,y = TP[e]
        mask[x,y] = 1
        if ~np.isnan(arr_thisElem[c]):
            img[x,y] = arr_thisElem[c]
            c+=1

    mask[X,Y] = 0
    with np.errstate(invalid='ignore',divide='ignore'):
        img /= mask

    # Plot similarity PT for elem
    fig = plt.figure(figsize=(18,7))
    gs = fig.add_gridspec(2,2,  width_ratios=(99, 1), height_ratios=(6, 4),
                  #    left=0, right=0.9, bottom=0.1, top=0.9,
                      wspace=0, hspace=0.2)

    ax = fig.add_subplot(gs[0, :])
    ax1 = fig.add_subplot(gs[1, :])
    cbar = fig.add_subplot(gs[0,1])

    cmap = sns.color_palette("magma", as_cmap=True)
    sns.heatmap(img,ax=ax,cbar_ax=cbar,
                vmin=0,vmax=np.nanmax(img),
                cmap=cmap)

    ax.set_title(f'Replaceability of {element} in chemical formulas, Year = {year}', fontsize=20)

    ax.axis('off')

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    #ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    # Plot barplot
    df = (pd.Series(arr_thisElem,index=elems_yr)
          .reset_index()
          .rename(columns={'index':'Element',
                           0:'Occurences'}))
    df = df[df.Element != element]

    df['color'] = (df['Occurences']/df['Occurences'].max()).apply(cmap)
    ax1.bar(x=range(df.shape[0]) ,
            height=df['Occurences'],
            color=df['color'],edgecolor = "k")

    ax1.set_xticks(range(df.shape[0]))
    ax1.set_xticklabels(df['Element'],fontsize=8)
    ax1.set_xlim(-1,df.shape[0])

    # Put the element's symbol at it's position
    if len(element)==1:  tab = 0.2
    else:                tab = 0.02
    ax.text(Y+tab,X+0.7,element,fontsize=17)
    
    if save: plt.savefig(save,dpi=400,bbox_inches='tight')
    
def plot_simMat_yr(simMat_yr,year,min_yr,save=False,raw=True,palette='magma_r',ordering=False,scale=15,show=True):
    """Plot similarity matrix for a given year
    year: which year to plot
    raw:  plot the raw normalized matrix (non-symmetric)
        if raw=False: plot symmetrized version
    ordering: lists position of elements ordered by atomic number. 
        e.g. [40,21,10] means H is in position 40, He in 21 and Li is 10th.
    """
    S = simMat_yr[year - min_yr].copy()
    
    # First change order, then clean empty rows + cols
    if type(ordering)!=bool: # Use new order
        indices = ['_' for i in range(103)]
        labels = indices.copy()
        for i,idx in enumerate(ordering):
            indices[idx] = i
            labels[idx] = elemList_AO[i]
        S = S[indices][:,indices]
        
    else: labels = elemList_AO
    
    # Remove non-existent elements (diag==0)
    diag = np.diag(S)
    isn = diag!=0
    S = S[isn][:,isn]
    n = isn.sum()
    diag = diag[isn]
    
    if raw:    P = (S/diag).T        
    else:
        Sum0 = S.sum(axis=0).reshape(-1,1).repeat(n,axis=1)
        Sum1 = S.sum(axis=1).reshape(1,-1).repeat(n,axis=0)
        P = np.sqrt(S**2/(Sum0*Sum1))
        
    ## Replace diagonal with 0, so that important features are evident
    inds = np.arange(0,n)
    P[inds,inds] = 0
    
    if show: 
        fig,ax = plt.subplots(figsize=(scale,scale))
        ax.set_title(f"Similarity matrix between elements ordered by atomic number, Year = {year}",
                     fontsize=20)
        sns.heatmap(P,ax=ax,cbar=False,
                    cmap=sns.color_palette(palette, as_cmap=True))
        
        labl = np.array(labels)[isn]      
        tick = [i+0.5 for i in range(len(labl))]
        ax.set_xticks(tick)
        ax.set_yticks(tick)
        ax.set_xticklabels(labl,fontsize=8)
        ax.set_yticklabels(labl,fontsize=8)
        
        if save: plt.savefig(save,dpi=400,bbox_inches='tight')

    return P



if __name__ == '__main__':
    main()
