#!/usr/bin/python3

"""
Put here everything related to visualizations of families (FEs)
"""

import plotly.graph_objects as go
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
from Components import loadData

from matplotlib import cm

cmap=cm.get_cmap("jet")

def Tanimoto(A1,A2):        
    A1,A2 = set(A1),set(A2)
    return len(A1.intersection(A2))/len(A1.union(A2))

FEs = loadData.FEs
elemList = loadData.elemList
simMats = loadData.simMats
min_yr=1800


# Compute DF with information of all FEs found over time.
FEs_df = pd.DataFrame()
for y in range(1800,2022):
    for fe in FEs[y]:
        FEs_df.loc[str(sorted(fe)),y] = 1

FEs_df = (FEs_df
          .replace(np.nan,0)
          .cumsum(axis=1)
          .sort_values(2021))

# Filter: Consider only FEs that have been found more than only twice in history.
FEs_df = FEs_df[FEs_df[2021]>2]

# Clean names (remove ], [, and ')
FEs_df.index = FEs_df.index.str.replace("]|\[|'","",regex=True) 

# Get year of newest element for each FE
elem_yr = pd.Series(3000,index=elemList)
FEs_df["yr_find_FE"] = 3000

for y in range(1800,2022):
    # List of CEs existing in year `yr`
    subs_elems = np.array(elemList)[simMats[y-min_yr].sum(axis=0)>0]
        
    # Update CEs that exist in year `y`, but current reported year is > `y`
    elem_yr[(elem_yr.index.isin(subs_elems)) & (elem_yr>y)] = y
        
    # Add year of first appearence for each FE.
    FE_yr = list(map(lambda x: re.sub("'|\[|]","", str(sorted(x))) , FEs[y]))
    FEs_df.loc[FEs_df.index.isin(FE_yr) &
               (FEs_df["yr_find_FE"]>y),"yr_find_FE"] = y  
    
# Map year information to each set, and select max year
FEs_df['yr_newest_elem'] = (FEs_df.index.to_series()
                            .str.split(', ')
                            .apply(lambda x: elem_yr[x].values.max()))

# Compute historical relevance, relative to each year.
# Difference between each year, and year of discovery of latest CE in each FE.
denom = (np.arange(1800,2022)[np.newaxis] - FEs_df.yr_newest_elem.values[:,np.newaxis])
denom*(denom>=0)  # Only count positive values.

# This division yields historical relevance for each FE, relative to each year.
hist_relev = FEs_df.loc[:,1800:2021]/(denom*(denom>=0)+1)

fig_closure = go.Figure()



def compMat(
        year_relat,
        FEs_df, 
        hist_relev,
        incl_ce=False, 
        show_nr = False,
        relev = "percent_max_elem",
        h=0.4,
        update=False,
):
    def closure(x,C):
        """
        Compute closure of set x given collection C.
            Closure of x given C, is smallest set y in C, so that x is contained in y. 
            if x not contained in C, return inf.
            Else, return |y|-|x|
        """
        if type(x)==str:    x = set(x.split(", "))
        tani = 0
        for y in C:
            t = Tanimoto(y,x)
            if t>tani:
                tani=t
        clos = np.inf

        for y in C:
            if set(x).issubset(y):
                c=len(y)-len(x)
                if c<clos:
                    clos=c

        if clos>0 and not np.isinf(clos):
            clos=1

        return clos

        for y in C:
            if set(x).issubset(y):
                if len(y-x)<clos:
                    clos = len(y-x)
        return clos

    # Get precomputed `relev`
    FEs_df = (FEs_df
              .merge(hist_relev[year_relat]
                     .rename("relev"), 
                     left_index=True, right_index=True))



    # Filter results to only those containing elems in `incl_ce`
    if incl_ce:
        FEs_df = FEs_df[FEs_df
                        .index.to_series()
                        .apply(lambda x: (set(x.split(', '))
                                          .issuperset(incl_ce)))]

    # For this subset of FEs, show only first `show_nr`
    if show_nr:
        FEs_df = FEs_df.iloc[-show_nr:]

    # For each FE in list, compute closure as defined above, given each years' collection.
    FEs_df = (pd.Series(range(1800,2022))    # For every year
              .apply(lambda y: (FEs_df
                                .index.to_series()        # To every FE
                                .apply(lambda x: closure(x, FEs[y]))  # Apply closure
                               )
                    )
              .T    # Make FE name be rows, year columns
              .rename(columns={i:y for i,y in enumerate(range(1800,2022))})     # Rename cols using year
              .merge(FEs_df["yr_find_FE"],
                     left_index=True, right_index=True)    # Merge with FEs_df to get relevances
             )

    hmap = (FEs_df
            .sort_values("yr_find_FE")
            .loc[:,1800:2021]
            .replace(np.inf,np.nan)
           )

    gg = np.zeros([hmap.shape[0]])
    for i in range(hmap.shape[0]):
        g = hmap.iloc[i].loc[~np.isnan(hmap.iloc[i])]
        gg[i] = np.sum(g==0)/len(g)

    def spl(x, i):
        """
        Reformat strings if they're too long. Add a "\n"
        """
        n=20
        if len(x)>n:
            x = x.split(", ")
            l = len(x)
            x = ", ".join(x[:l//2+1]) + "<br>" + ", ".join(x[l//2+1:])
        x = "{}".format(x)
        return x

    # Clip values not to show very high values. Anything above 4 is same color (>4)
    hmap[hmap>4] = 4

    labls = hmap.index.to_series()
#    labls = [i
#             #spl(lab,i)
#             for i,lab in enumerate(labls)]


    trace = go.Heatmap(y=labls,
                       x=hmap.columns,
                       z=hmap.values,
                       colorscale=[
                           [0,'#000000'],
                           [1,'#9b2850']
                       ],
                       showscale=False
                       )

    fig_closure.update_layout(
        yaxis = dict(scaleanchor = 'x',
                     autorange='reversed',
                     tickvals=list(range(len(labls))),
                     ticktext=["" for i in labls]
                     ),
        xaxis = dict(side='top',
                     tickangle=-90,
                     tickvals=np.arange(1800,2022,10)),
        margin=dict(l=0,r=0,b=0,t=0),
        plot_bgcolor='rgba(0,0,0,0)'
    )
    fig_closure.update_xaxes(showgrid=True, gridwidth=1, gridcolor='black')

    if update:
        fig_closure.update_traces(trace)
    else:
        fig_closure.add_traces(trace)

    return fig_closure


closure_fig = compMat(
    2021,
    FEs_df,
    hist_relev,
    incl_ce = {'K'},
    show_nr=10,
)
