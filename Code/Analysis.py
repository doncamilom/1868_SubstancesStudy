#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from TPs import TP
import seaborn as sns; sns.set() ; sns.set_style('white')
import pandas as pd
from itertools import chain

import re
import bz2
import pickle

from matplotlib.animation import FuncAnimation, PillowWriter
from scipy import sparse as sp
import matplotlib.patheffects as PathEffects

def main():
    global Matches

    ####### First load all the data that was produced during preprocessing.
    #######
    dataPath = "./Data/"

    # Load AllMatches.bin
    Match_file = bz2.BZ2File(dataPath+'AllMatches.bin', 'r')
    Matches = pickle.load(Match_file)

    # Load AllRs.npz
    Rs = sp.load_npz(dataPath+'AllRs_clean_sparse.npz')

    print(Rs.shape)
    # Load element list
    elemList = []
    with open(f"{dataPath}/ElementList.txt",'r') as f:
        for line in f:
            elemList.append(line.strip())
            
    print(elemList[:10])

    #######
    ####### Done loading data

    # Generate data for animation
    year_span = 30
    dataGIF = np.zeros((year_span,7,32))
    for i in range(year_span):
        dataGIF[i] = getmat_heat(Matches,i*2+1700)


    

    grid_kws = {'width_ratios': (0.99, 0.01), 'wspace': 0.0}
    fig, (ax, cbar_ax) = plt.subplots(1, 2, gridspec_kw = grid_kws, figsize = (18, 4.5))
    anim = FuncAnimation(fig = fig, func = animate_func, fargs = (dataGIF,ax,cbar_ax,) ,frames = year_span, interval = 50, blit = False)

    writer = PillowWriter(fps=2)  
    anim.save("./PT_evol.gif", writer=writer) 





#########################
##### Auxiliary functions
#########################

# Function to construct a TPR given a TP layout and a list of elements.
def getTable(TP, elems, fill_data,emptyTab):
    """This gets as input: a PT and a list of elements. 
        Returns a 2D array coloured at the positions of these elements"""
    for elem,fill in zip(elems,fill_data):
        x,y = TP[elem]
        emptyTab[x,y] = fill
    return emptyTab

def getFormula(R,elemList):
    """Convert R sparse vector into string composition:
    sparse([0,1,0,...,4,6]) --> Ti4X6 for instance """

    form = ''
    for ind,n in zip(R.indices[:-1],R.data[:-1]):
        form += elemList[ind] 
        if n != 1:
            form += str(int(n))
    return form + f'X{int(R.data[-1]) if R.data[-1]!=1 else ""}'

# Do stats

def makeBoxPlot(dt,TP):
    """ Use redefined PT, so that it matches with dt dimensions"""
    names = []  #Will be the names of df columns. Will appear in the box's labels
    inv_PT = {v: k for k, v in TP.items()}  #Inverted map (positions -> element name)
    for y in range(dt.shape[1]):
        for x in range(dt.shape[2]):         
            key = (y,x)
            if key in inv_PT.keys():     names.append(inv_PT[key])  #Return element at this position
            else:                        names.append('None')       #If no element, return None
    
    dt = pd.DataFrame(dt.reshape(dt.shape[0],-1),columns=names)  #Reshape into [N_Tables, TPdims] DF
    dt = dt.dropna(axis=1,how='all')    #Select only entries for which elements exist (non-nan)
    fig,ax = plt.subplots(figsize=(18,5))
    
    sns.boxplot(data=dt,ax=ax)#,saturation=1,fliersize=4)
    plt.show()
    
def makeTPPlot(dt,norm,title="Mean Euclidean Distance",save=False,elementLabels=False,fs=22,PT=TP):
    with np.errstate(invalid='ignore'):
        dt = np.nan_to_num(dt*norm).sum(axis=0)/norm.sum(axis=0)
        
    grid_kws = {'width_ratios': (0.99, 0.01), 'wspace': 0.0}
    fig, (ax, cbar_ax) = plt.subplots(1, 2, gridspec_kw = grid_kws, figsize = (18, 4.5))
    
    if title:     ax.set_title(title,fontsize=20)
        
    sns.heatmap(dt,
        ax = ax,
        cbar = True,
        cbar_ax = cbar_ax,
        vmin = np.nanmin(dt),  #max and min ignoring nans
        vmax = np.nanmax(dt))

    putLabelsInTables(~np.isnan(dt),ax,fs=18,pt=TP,pathEffects=True)
    ax.set_xticks([])
    ax.set_yticks([])

    plt.tight_layout()
    ax.axis('off')
    if save:    plt.savefig(save,dpi=300,bbox_inches='tight')
    plt.show()

def putLabelsInTables(ta,ax,fs=11,pt=TP,pathEffects=False):
    """
    Set the elements' chemical symbol on top of the table for visualization
    ta is the table in question, ax is the ax where the table is to be drawn.
    """
    for element in pt.keys():
        x,y = pt[element]
        if ta[x,y] >= 1:
            # Put the element's symbol at it's position
            if len(element)==1:  dy = -0.25
            else:                dy = -0.15
            txt = ax.text(y-dy,x+0.64,element,fontsize=fs,color='k') 
            if pathEffects:
                txt = ax.text(y-dy,x+0.64,element,fontsize=fs,color='w') 
                txt.set_path_effects([PathEffects.withStroke(linewidth=2.5, foreground='k')])

def meanDist(Matches,max_year):
    """For more efficiency, compute only one of the distances instead of both.
    For that, function has to be changed."""
    # First, filter out data. We only care about cmpnds discovered before max_year
    filt_mat = list(map(lambda x: [x[0][i] for i in range(len(x[0])) if x[1][i] <= max_year] 
                       if (x[1]<=max_year).sum()>1 else None,   Matches))
    filt_mat = [i for i in filt_mat if i is not None]
    
    def dists(match):
        """Get a mapping [elems] --> [distances]"""
        # /np.array([6,32]) is an optional normalization
        xy = np.array(list(map(lambda x: TP[x], match)))  /np.array([6,1])
        N = len(match)
        A = np.repeat(xy,N).reshape(N,2,N)
        vert,hori = np.abs(A - A.T).sum(axis=0)

        return vert,hori, N
    
    # Translate above mappings into arrays
    arr = np.zeros((2,len(filt_mat),7,32))
    Ns  = np.zeros(len(filt_mat))
    for i,match in enumerate(filt_mat):
        vert,hori,N = dists(match)
        getTable(TP, match, vert, arr[0,i])  # Vert dist
        getTable(TP, match, hori, arr[1,i])  # Hori dist
        Ns[i] = N
        
    with np.errstate(invalid='ignore'):
        """Divide by boolean array: is element in a particular position?
        arr.sum in axis 0 == horiz + vert dists. 
        If elem exists, either can be zero, but not both at the same time 
        otherwise element'd be alone in table, which is not possible."""
        norm = (arr.sum(axis=0)!=0)*Ns.reshape(1,-1,1,1)
        arr /= norm
    return arr,norm[0]   # arr = [vert,hori]

def getmat_heat(Matches,i):
    """Obtain distances for data up to a given year"""
    # First calculate distances
    dis,norm = meanDist(Matches,i)
    vert,hori = dis
    with np.errstate(invalid='ignore'):
        dt = np.nan_to_num(hori*norm).sum(axis=0)/norm.sum(axis=0)
    return dt

def animate_func(i,data,ax,cbar_ax):
    ax.cla()
    ax.set_title(f"Mean Horizontal Distance. Year = {i*2+1700}",fontsize=20)
    
    sns.heatmap(data[i, ...],
            ax = ax,
            cbar = True,
            cbar_ax = cbar_ax,
            vmin = np.nanmin(data),  #max and min ignoring nans
            vmax = np.nanmax(data))

    putLabelsInTables(~np.isnan(data[i]),ax,fs=18,pt=TP,pathEffects=True)
    ax.set_xticks([])
    ax.set_yticks([])

def makeBoxPlot(dt,TP):
    dataList = []
    elems = []
    # Extract distance data for every element
    for el in TP.keys():
        x,y=TP[el]
        dt_i = dt[:,x,y]
        # Only append if there is at least a non-nan entry
        if np.isnan(dt_i).sum()<dt_i.shape[0]: 
            dataList.append(dt_i) 
            elems.append(el)
            
    dt = pd.DataFrame(np.array(dataList).T,columns=elems)
    fig,ax = plt.subplots(figsize=(18,5))
    sns.boxplot(data=dt,ax=ax)#,saturation=1,fliersize=4)
    plt.show()

if __name__ == '__main__':
    main()
