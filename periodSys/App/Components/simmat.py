#!/usr/bin/python3

"""
This is code to plot a similarity matrix using plotly
"""

import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np
from Components import loadData

year = 2020
ordering=False
elemList=loadData.elemList

S = loadData.simMats[year-loadData.min_yr].copy()

# First change order, then clean empty rows + cols
if type(ordering)!=bool: # Use new order
    indices = ['_' for i in range(103)]
    labels = indices.copy()
    for i,idx in enumerate(ordering):
        indices[idx] = i
        labels[idx] = elemList[i]
    S = S[indices][:,indices]
    
else: labels = elemList

# Remove non-existent elements (diag==0)
diag = np.diag(S)
isn = diag!=0
S = S[isn][:,isn]
n = isn.sum()
diag = diag[isn]

Sum0 = S.sum(axis=0).reshape(-1,1).repeat(n,axis=1)
Sum1 = S.sum(axis=1).reshape(1,-1).repeat(n,axis=0)
P = np.sqrt(S**2/(Sum0*Sum1))
    
## Replace diagonal with 0, so that important features are evident
inds = np.arange(0,n)
#P[inds,inds] = 0




def plot_simMat_yr(simMat_yr,year,min_yr,save=False,raw=True,cmap=False,ordering=False,scale=15,show=True,elemList=False):
    
    if show: 
        fig,ax = plt.subplots(1,2,figsize=(scale,scale),
                              gridspec_kw={"width_ratios":[100,1],"wspace":0.05})
        ax[0].set_title("Similarity matrix between elements ordered by atomic number, Year = {}".format(year),
                     fontsize=20)
        sns.heatmap(P,norm=LogNorm(),ax=ax[0],cbar_ax=ax[1],cmap=cmap)
        
        labl = np.array(labels)[isn]      
        tick = [i+0.5 for i in range(len(labl))]
        ax.set_xticks(tick)
        ax.set_yticks(tick)
        ax.set_xticklabels(labl,fontsize=8)
        ax.set_yticklabels(labl,fontsize=8)
        
        if save: plt.savefig(save,dpi=400,bbox_inches='tight')

    return P


def colorbar(zmin):
    lowlabel = np.ceil(np.log10(zmin))
    tickvals = np.arange(lowlabel,0,1)
    return dict(
        tickmode = "array",
        tickvals = list(tickvals)+[0],
        ticktext = ['10^-{}'.format(i) for i in -tickvals.astype(int)] + ['10^0'],
        thickness = 8,
        len = 0.7,
        x = 0.99
    )

P /= P.max()
zmin = P[P>0.].min()

fig = go.Figure()
fig.add_trace(go.Heatmap(x=elemList,y=elemList,z=np.log10(P), text = P,
                         colorscale='Jet',
                         colorbar = colorbar(zmin),
                         hovertemplate =
                            "<b>%{x}~%{y}</b><br>" +
                            "<i>%{text:.4f}</i>" +
                            "<extra></extra>",
                        hoverlabel = dict(
                            bgcolor="white",
                            font_size=20,
                            font_family="Rockwell"
                        ),

    ))

fig.update_layout(yaxis = dict(scaleanchor = 'x'),
                  margin=go.layout.Margin(
                        l=0, #left margin
                        r=0, #right margin
                        b=0, #bottom margin
                        t=0  #top margin
                ))

fig.update_layout(plot_bgcolor='rgba(0,0,0,0)')
fig['layout']['yaxis']['autorange'] = "reversed"
