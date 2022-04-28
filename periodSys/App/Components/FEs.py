#!/usr/bin/python3

"""
Put here everything related to visualizations of FEs. 
Same or different window?
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

