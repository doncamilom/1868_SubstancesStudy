#! /usr/bin/env python3

import numpy as np
import multiprocessing as mp
import pandas as pd
from itertools import chain
from TPs import TP
import re
import bz2
import pickle
from scipy import sparse as sp
from time import time

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import LogNorm
    from matplotlib import cm
except: pass

def getElemList(dataPath):
    with open(dataPath + 'ElementList.txt','r') as f:
        elemList=f.read().split('\n')
    elemList.remove('')
    return elemList


def main():
    global elemList,Matches,NP

    dataPath = root+'Data/'   

    # Load element list
    elemList = getElemList(dataPath)
    print("Number of elements: {}".format(len(elemList)))

    # Load matches_pickle.bin
    Match_file = bz2.BZ2File(dataPath+'matches_pickle.bin', 'r')
    Matches = pickle.load(Match_file)

    t0 = time()
    sim_mat = calc_simMats_yearly()
    print("Time: {:.3f} s".format(time()-t0))

    np.save(dataPath+'simMat.npy',sim_mat)


########
# Functions necessary for computation of similarity matrices
########

def getTable(TP, elems, fill_data, emptyTab):
    """
    This func maps from an element list into an array filled with years at the elements positions.
    e.g. input [Al, Kr, Cl] and [1860, 2010, 1952]. This returns [1860, 0, 0, ... , 1952, 0, ... , 2010, .., 0]
    """
    indx=list(map(elemList.index,elems))
    for i,fill in zip(indx,fill_data):
        emptyTab[i] = fill
    return emptyTab


def getSimilarities_yr(Tyr,element,year):
    """Get array of the similarities between elements and the given element
    based on 'replaceability' in chemical formulas, for the CS existent at a given year."""
    # Filter out by year
    T = (Tyr <= year)&(Tyr>0)
    
    # When filtering by year, we need only sum over those tables where at least 2 elements are present!!!
    T = T[T.sum(axis=1)>1]

    # Now start actually calculating simMats
    X = elemList.index(element)   # Get element's position in list
    # similarity: number of times other elements appear in same table
    simil = T[T[:,X]].sum(axis=0)*1.
 
    data_arr = pd.Series(simil,index=elemList).rename(element)
    return data_arr


def getTyr(i):
    """ Compute Tyr for entries in Matches from i to i+1 (in units of len(Matches)/NP"""
    ch_sz = len(Matches)//NP
    lim1,lim2 = ch_sz*i, ch_sz*(i+1)
    if i<NP-1:  sub_matches=Matches[lim1:lim2]
    else:       sub_matches=Matches[lim1:]
    
    Tyr = np.zeros((len(sub_matches),len(elemList)))      # This gets the arrays filled with years
    for j,match in enumerate(sub_matches):
        getTable(TP, match[0], match[1], Tyr[j] )  # Match[1] to fill with years
    return Tyr


def simMat(year):
    mat = pd.DataFrame(columns=[e for e in TP.keys() if e in elemList])
    for e in elemList: # Calculate in order, but subjected to what's in elemList
        mat[e] = getSimilarities_yr(Tyr,e,year)
    return year,mat  # Better return array, not dataframe


def calc_simMats_yearly():
    """Calculate all simMats yearly.
    First load all the data that was produced during preprocessing.
    Takes as input: dataPath, and number of elements in dataset"""
    global Tyr

    miny,maxy=1800,2022
    print(miny,maxy)

    with mp.Pool(processes=NP) as pool:
        res = pool.map(getTyr, range(NP))
    print("Done computing Tyr")

    # Concat results
    Tyr = np.concatenate(res,axis=0)

    # Generate data yearly
    with mp.Pool(processes=NP) as pool:
        results = [pool.apply_async(simMat,args=(yr,)) for yr in range(miny,maxy+1)]
        res = [r.get() for r in results] 
    
    # Collect all results into a dict so we can sort them by year in a following step
    simMats={yr:mat for yr,mat in res}

    year_span = maxy-miny+1
    numElems = len(elemList)
    simMats_yr = np.zeros((year_span,numElems,numElems))

    for i,yr in enumerate(range(miny,maxy+1)):
        simMats_yr[i] = simMats[yr][elemList] # Resort by order of elemList

    return simMats_yr


##########
# Other functions used in the notebook
##########

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

    ax.set_title('Replaceability of {} in chemical formulas, Year = {}'.format(element,year), fontsize=20)

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
    
def plot_simMat_yr(simMat_yr,year,min_yr,save=False,raw=True,cmap=False,ordering=False,scale=15,show=True,EL=False):
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
            labels[idx] = EL[i]
        S = S[indices][:,indices]
        
    else: labels = EL
    
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
        
    if show: 
        fig,ax = plt.subplots(1,2,figsize=(scale,scale),
                              gridspec_kw={"width_ratios":[100,1],"wspace":0.05})
        ax = ax.ravel()
        ax[0].set_title("Similarity matrix between elements ordered by atomic number, Year = {}".format(year),
                     fontsize=20)
        sns.heatmap(P,norm=LogNorm(),ax=ax[0],cbar_ax=ax[1],cmap=cmap)
        
        labl = np.array(labels)[isn]      
        tick = [i+0.5 for i in range(len(labl))]
        ax[0].set_xticks(tick)
        ax[0].set_yticks(tick)
        ax[0].set_xticklabels(labl,fontsize=8)
        ax[0].set_yticklabels(labl,fontsize=8)
        
        if save: plt.savefig(save,dpi=400,bbox_inches='tight')

    return P

if __name__ == '__main__':
    global root
    import argparse
    parser = argparse.ArgumentParser(description="Compute similarity matrices.")
    parser.add_argument("--NumProc","-n", required=True, type=int,
                        help="Number of processors to use for computations.")
    parser.add_argument("--inp_dir","-i", required=False, type=str, default="./",
                        help="Root directory for all computations")
    args = parser.parse_args()

    NP=args.NumProc
    root=args.inp_dir

    main()
