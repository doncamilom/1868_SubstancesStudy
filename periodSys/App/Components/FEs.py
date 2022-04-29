#!/usr/bin/python3

"""
Put here everything related to visualizations of FEs. 
Same or different window?
"""

import plotly.graph_objects as go
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
from Components import loadData

from matplotlib import cm
cmap=cm.get_cmap("jet")

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

# Filter: Consider only FEs that have been found more than only once in history.
FEs_df = FEs_df[FEs_df[2021]>1]

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
def HeatmapClosure(year_relat,
                   FEs_df, hist_relev,
                   incl_ce=False,
                   notincl_ce=False,
                   thresh_relev = 0.1,
                   update = True,):
    """
    Compute relevancy descriptors for each FE
    based on year of discovery of FE/discovery of newest element in FE.

    Filter results using:
        incl_ce: set. Analize only FEs containing the CEs in this set. e.g. {'H','Na'}
        notincl_ce: set. Analize only FEs NOT containing the CEs in this set. e.g. {'K'}
        thresh_relev: Limit FEs analyzed to only those with relevance > thresh_relev

    Returns pd.DataFrame with all data of relevance and yearly closure for filtered FEs.
    """
    def closure(x,C):
        """
        Compute closure of set x given collection C.
            Closure of x given C, is smallest set y in C, so that x is contained in y.
            if x not contained in C, return inf.
            Else, return |y|-|x|
        """
        if type(x)==str:    x = set(x.split(", "))
        clos = np.inf

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

    # First filter FEs, so only relevants are considered
    FEs_df = FEs_df[FEs_df["relev"]>thresh_relev]

    # Filter results to only those containing elems in `incl_ce`
    if incl_ce:
        FEs_df = FEs_df[FEs_df
                        .index.to_series()
                        .apply(lambda x: (set(x.split(', '))
                                          .issuperset(incl_ce)))]

    # Filter results to include only FEs NOT containing elems in `notincl_ce`
    if notincl_ce:
        FEs_df = FEs_df[FEs_df
                        .index.to_series()
                        .apply(lambda x: (len(set(x.split(', '))
                                          .intersection(notincl_ce)) == 0))]


    # For each FE in list, compute closure as defined above, given each years' collection.
    FEs_df = (pd.Series(range(1800,2022))    # For every year
              .apply(lambda y: (FEs_df
                                .index.to_series()        # To every FE
                                .apply(lambda x: closure(x, FEs[y]))  # Apply closure
                               )
                    )
              .T    # Make FE name be rows, year columns
              .rename(columns={i:y for i,y in enumerate(range(1800,2022))})     # Rename cols using year
              .merge(FEs_df["relev"],
                     left_index=True, right_index=True)    # Merge with FEs_df to get relevances
             )

    # Parse title
    if not incl_ce and not notincl_ce:
        title = f"Relative to {year_relat}. Threshold: {thresh_relev}.\n"
    elif incl_ce and not notincl_ce:
        title = f"Relative to {year_relat}. Threshold: {thresh_relev}.\nRestriction: Containing {incl_ce}\n"
    elif not incl_ce and notincl_ce:
        title = f"Relative to {year_relat}. Threshold: {thresh_relev}.\nRestriction: Not containing {notincl_ce}\n"
    else:
        title = f"Relative to {year_relat}. Threshold: {thresh_relev}.\n\
Restriction: Containing {incl_ce} but not {notincl_ce}\n"

    title = re.sub("'|\{|\}","",title)


    hmap = (FEs_df
            .sort_values("relev",
                         ascending=False)
            .loc[:,1800:2021]
            .replace(np.inf,np.nan)
           )


    trace = go.Heatmap(y=hmap.index.values,x=hmap.columns,
                         z=hmap.values,
                         colorscale='Jet',
                         colorbar = dict(
                                tickmode = "array",
                                thickness = 8,
                                len = 0.7,
                                x = 1,
                            ),
                #         text = p,
                #         hovertemplate =
                #            "<b>%{x}~%{y}</b><br>" +
                #            "<i>%{text:.4f}</i>" +
                #            "<extra></extra>",
                #        hoverlabel = dict(
                #            bgcolor="white",
                #            font_size=20,
                #            font_family="Rockwell"
                        #),
        )

    fig_closure.update_layout(yaxis = dict(scaleanchor = 'x',
                              autorange='reversed'),
                      xaxis = dict(side='top',tickangle=0),
                      margin=dict(l=0,r=0,b=0,t=0),
                      plot_bgcolor='rgba(0,0,0,0)'
                      )

    if update:  fig_closure.update_traces(trace)
    else:       fig_closure.add_traces(trace)

    return fig_closure



def PlotHistoricalRelev(incl_ce=False,
                        notincl_ce=False,thresh=0.1):
    """
    Plot historical relevance for each FE, as a function of time.
        Consider only FEs containing all CEs in `incl_ce`.
                 only FEs containing none of the CEs in `notincl_ce`.
                 only FEs whose Historical Relevance is > `thresh` at some point in history.
    """
    relev = hist_relev.copy()

    # Filter results to only those containing elems in `incl_ce`
    if incl_ce:
        relev = relev[(relev
                       .index.to_series()
                       .apply(lambda x: (set(x.split(', '))
                                         .issuperset(incl_ce))
                             )
                      )]

    # Filter results to include only FEs NOT containing elems in `notincl_ce`
    if notincl_ce:
        relev = relev[(relev
                       .index.to_series()
                       .apply(lambda x: (not set(x.split(', '))
                                         .issuperset(notincl_ce))
                             )
                      )]

    relev = (relev
             .sort_values(2021,ascending=False))


    N=relev.shape[0]
    fig,ax=plt.subplots(N,1,figsize=(17,N*1.2),sharex=True,gridspec_kw={"hspace":0.3})

    for num,i in enumerate(relev.index):
        relev.columns = range(1800,2022)
        ax[num].plot(relev.loc[i],color='b')

        min_yr = relev.loc[i,~relev.loc[i].isna()].index.min()
        y_text = (relev.loc[i,2021])


        ax[num].text(2022,y_text,relev.loc[i].name,fontsize=20,rotation=0)

        for _ in ['top','left','right']:    ax[num].spines[_].set_visible(False);
        ax[num].set_yticks([])
        ax[num].grid("x")

        ax[num].set_ylabel(num+1)
        #ax[num].set_ylim(0,1.05)
        ax[num].set_ylabel('')

    ax[num].set_xlim(1800,2022)
    return fig


closure_fig = HeatmapClosure(2021,
                           FEs_df,hist_relev,
                           incl_ce = {'K'},
                           notincl_ce = {'H'},
                           thresh_relev=0.,update=False)


relev_fig = PlotHistoricalRelev(incl_ce = {'Kr'},
                                thresh=0.6)#, update=False)
