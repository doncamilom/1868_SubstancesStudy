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
                  margin=dict(l=0,r=0,b=0,t=0),
                  plot_bgcolor='rgba(0,0,0,0)'
                  )

fig['layout']['yaxis']['autorange'] = "reversed"
