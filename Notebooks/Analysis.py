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

#from matplotlib.animation import FuncAnimation, PillowWriter
from scipy import sparse as sp
#import matplotlib.patheffects as PathEffects

def main():
    global Matches

    ####### First load all the data that was produced during preprocessing.
    #######
    dataPath = "../Code/OnlyStrs/Data/"
#    dataPath = "./Data/"

    # Load matches_pickle.bin
    Match_file = bz2.BZ2File(dataPath+'matches_pickle.bin', 'r')
    Matches = pickle.load(Match_file)

    # Load Rs_sparse.npz
    Rs = sp.load_npz(dataPath+'Rs_sparse.npz')

    print(Rs.shape)

    # Load element list
    elemList = []
    with open(f"{dataPath}/ElementList.txt",'r') as f:
        for line in f:
            elemList.append(line.strip())

    #######
    ####### Done loading data

    miny,maxy = getMinMax(Matches)
    print(miny,maxy)

    # Generate data for animation
    year_span = maxy-miny+1
    dataGIF = np.zeros((year_span,7,32))

    for i,yr in enumerate(range(miny,maxy+1)):
        dataGIF[i] = getmat_heat(Matches,yr)

    np.save(dataPath+'history_HND.npy',dataGIF)

    
#    grid_kws = {'width_ratios': (0.99, 0.01), 
#                'height_ratios':(0.3,0.7),
#                'wspace': 0.05,'hspace':0.23}
#    fig, ax = plt.subplots(2, 2, gridspec_kw = grid_kws, figsize = (13, 4))
#
#    anim = FuncAnimation(fig = fig, func = animate_func, fargs = (dataGIF,ax,np.nanmin(dataGIF),np.nanmax(dataGIF),) ,
#                            frames = 100, interval = 50, blit = False)
#
#    print("Saving animation...")
#    writer = PillowWriter(fps=2)  
#    anim.save("./PT_evol.gif", writer=writer) 





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
    

def makeTPPlot(dt,title="Mean Euclidean Distance",save=False,elementLabels=False,fs=22,PT=TP,min_scale=0,max_scale=10,ax=0):

    cmap = sns.color_palette("magma", as_cmap=True)
    
    # Make boxplot
    epc = np.sum(~np.isnan(dt),axis=0) # epc: number of elements per column
    ys = np.nan_to_num(dt).sum(axis=0)/(epc+(epc==0))

    # Plot barplot
    df = (pd.Series(ys)
          .reset_index()
          .rename(columns={'index':'Element',
                           0:'Mean'}))
    
    df['color'] = (df['Mean']/max_scale).apply(cmap)
    ax[0,0].clear()
    ax[0,0].bar(x=range(df.shape[0]) ,
            height=df['Mean'],
            color=df['color'],edgecolor = "k")
    ax[0,0].set_ylim(0,max_scale)
    
    # Invisibilize ax[0,1]
    ax[0,1].set_xticks([])
    ax[0,1].set_yticks([])
    ax[0,1].axis('off')
    
    # Labels for boxplot
    labels = ['IA','IIA','   l- - -']+ ['- - -']*6 + ['-IIIB-'] +['- - -']*6 + \
                ['- - -l  ','IVB','VB','VIB','VIIB','  l- -','-VIIIB-', '- -l  ','IB','IIB','IIIA',\
                'IVA','VA','VIA','VIIA','VIIIA']
    
    ax[0,0].set_xticks(range(32))
    ax[0,0].set_xticklabels(labels,fontsize=9,color='k')
    ax[0,0].set_xlim(-0.5,31.5)
    
    if title:     ax[0,0].set_title(title,fontsize=20)
    
    sns.heatmap(dt,
        ax = ax[1,0],
        cbar = True,
        cbar_ax = ax[1,1],
        vmin = min_scale,  #max and min ignoring nans
        vmax = max_scale,
        cmap = cmap)

    putLabelsInTables(~np.isnan(dt),ax[1,0],fs=12,pt=TP,pathEffects=True)
    ax[1,0].set_xticks([])
    ax[1,0].set_yticks([])

    ax[0,0].spines['top'].set_visible(False)
    ax[0,0].spines['right'].set_visible(False)
    ax[0,0].spines['left'].set_visible(False)
    
    
    if save:    plt.savefig(save,dpi=300,bbox_inches='tight')
    #plt.show()

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
    if max_year>0:
        filt_mat = list(map(lambda x: [x[0][i] for i in range(len(x[0])) if int(x[1][i]) <= max_year] 
                           if (np.array(x[1],dtype=np.int)<=max_year).sum()>1 else None,   Matches))
        filt_mat = [i for i in filt_mat if i is not None]
    else:
        filt_mat = [x[0] for x in Matches]

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

def animate_func(i,data,ax,min_scale,max_scale):
    
    makeTPPlot(data[i],title=f"Mean Horizontal Distance {i*2+1700}",PT=TP,fs=18,max_scale=max_scale,min_scale=min_scale,ax=ax)

def getMinMax(Matches):
    years = list(chain(*[x[1] for x in Matches] ))
    years = list(map(int,years))

    return min(years),max(years)


if __name__ == '__main__':
    main()
