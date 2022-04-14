#!/usr/bin/python3

"""
This is code to plot a similarity matrix using plotly
"""

import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np
from Components import loadData

# Compute all (symmetric) similarity matrices at once
P = loadData.simMats.copy()
Sum0 = P.sum(axis=1).reshape(P.shape[0], P.shape[1], 1).repeat(103,axis=2)
Sum1 = P.sum(axis=2).reshape(P.shape[0], 1, P.shape[1]).repeat(103,axis=1)

with np.errstate(divide='ignore', invalid='ignore'):
    P = np.sqrt(P**2/(Sum0*Sum1))
    P /= np.nanmax(P)

zmin = P[P>0.].min()
min_yr = 1800


def getElemList(dataPath):
    elemList = []
    with open("{}/ElementList.txt".format(dataPath),'r') as f:
        for line in f:
            elemList.append(line.strip())
    return elemList

def getSimMat(year):
    return P[year-min_yr]

def colorbar(zmin):
    lowlabel = np.ceil(np.log10(zmin))
    tickvals = np.arange(lowlabel,0,1)
    return dict(
        tickmode = "array",
        tickvals = list(tickvals)+[0],
        ticktext = ['10^-{}'.format(i) for i in -tickvals.astype(int)] + ['10^0'],
        thickness = 8,
        len = 0.7,
        x = 0.99,
    )

elemList = getElemList('../Data')

fig = go.Figure()
def plotSimMat(year, update=True):
    np.seterr(divide='ignore')
    trace = go.Heatmap(x=elemList,y=elemList,
                         z=np.log10(getSimMat(year)), text = getSimMat(year),
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
        )

    if update:      fig.update_traces(trace)
    else:           fig.add_traces(trace)

    fig.update_layout(yaxis = dict(scaleanchor = 'x',autorange='reversed'),
                      xaxis = dict(side='top',tickangle=0),
                      margin=dict(l=0,r=0,b=0,t=0),
                      plot_bgcolor='rgba(0,0,0,0)'
                      )
    return fig

fig = plotSimMat(2022,update=False)