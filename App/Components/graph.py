from jax import grad, jit

import pickle
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import plotly.graph_objects as go

def M2(A1,A2):         
    A1,A2 = set(A1),set(A2)
    return len(A1.intersection(A2))/np.sqrt(len(A1)*len(A2))

class HistGraph():
    """The whole graph, not limited to MIN/MAX_YR. 
    Methods plot and optimize are restricted tho. So opt can be carried out in windows.
    
    nodes: dict(yr: [list of groups])
    """
    def __init__(self,nodes):
        # First df for nodes (groups) and initialize x var for each
        node_idd = [(f"{yr}_{id}",nodes[yr][id],yr) for yr in nodes.keys() # List of (id,grp) for all groups
                                    for id,_ in enumerate(nodes[yr])] 
        node_id,node_unwr,node_yr = list(zip(*(node_idd))) # Unwrap in full list of ids and list of groups
        self.node_df = pd.DataFrame({"name":node_id,
                                     "x":np.random.random(len(node_id))*2,
                                     "elems":node_unwr,
                                     "yr":node_yr})
        self.node_df["z"] = np.random.random(len(node_id))*2 # Add z var for flexibility
        map_name_ind = dict(zip(node_id,self.node_df.index))
        
        # Second, df for edges
        edges = []
        for yr in nodes.keys():
            for id1,gr1 in enumerate(nodes[yr]):  
                # Include connections between different years
                if yr<max(nodes.keys()): # So yr+2 won't be a problem
                    for id2,gr2 in enumerate(nodes[yr+2]): 
                        m2 = M2(gr1,gr2)
                        dist = 1-m2 + 1.2 # 1.2 for scaling. if 1-dist==1 -> k~2.2~sqrt(5): diagonal 
                        v = [f"{yr}_{id1}",f"{yr+2}_{id2}",yr,yr+2,dist,m2]  
                        edges.append(v)

                # Include connections between groups of same year 
                id2=id1
                for gr2 in nodes[yr][id1+1:]:
                    id2+=1
                    # Add 0.05 so groups are still distinguishable.
                    m2=M2(gr1,gr2)
                    dist = 1-m2 + 0.05
                    v = [f"{yr}_{id1}",f"{yr}_{id2}",yr,yr,dist,m2] 
                    edges.append(v)   
        
        self.edges = pd.DataFrame(edges,columns=['yr_id1','yr_id2','yr1','yr2','k','m2'])
        # Map the indices for then mapping the x coords
        self.edges["grp1_ind"] = self.edges["yr_id1"].apply(lambda x: map_name_ind[x])
        self.edges["grp2_ind"] = self.edges["yr_id2"].apply(lambda x: map_name_ind[x])  

    def Optimize(self,MIN_YR,MAX_YR,EPOCHS=1000,LR=1e-3,v=500):
        """Optimize locally points in the range [MIN_YR,MAX_YR] (including both!)
        v controls verbosity: Frequency of output, cost is printed every `v` iterations.
        v=0: no output"""
        
        # Query only edges involving the year range of interest
        quer1 = (self.edges.yr1>=MIN_YR)&(self.edges.yr2>=MIN_YR)
        quer2 = (self.edges.yr1<=MAX_YR)&(self.edges.yr2<=MAX_YR)
        opt_edges = self.edges[quer1&quer2]
        
        # Get inputs for energy
        grp1_ind,grp2_ind = opt_edges.loc[:,["grp1_ind","grp2_ind"]].values.T
        yr1,yr2 = opt_edges.loc[:,['yr1','yr2']].values.T  #ys
        z1,z2 = self.node_df["z"].values[grp1_ind],self.node_df["z"].values[grp2_ind]
        ks = opt_edges.loc[:,'k'].values

        # Opt loop (just basic grad descent) ## modify to use more powerful methods
        costs = []
        for e in range(EPOCHS):
            df = dEdx(self.node_df["x"].values,grp1_ind,grp2_ind,yr1,yr2,z1,z2,ks)  # Calc grad
            self.node_df["x"] -= LR*df  # Update xs_all using grad

            E = energy(self.node_df["x"].values,grp1_ind,grp2_ind,yr1,yr2,z1,z2,ks)
            costs.append(E)
            if v!=0: 
                if e%v==0:  print(e,E)

    def plotGraph(self,MIN_YR,MAX_YR,THRESH=0.0,seed=False,
                  figsize=(20,5),title="",save=False):
        """Plot nodes and edges in range [MIN_YR, MAX_YR]. 
        Edges are plotted above THRESH [0,1], and color intensity is based on k.
        
        seed: list of [grp_ids] used as a seed to propagate the graph. 
            Only significant connections (M2>THRESH) to such groups will be plotted.
        """
        ## Query on seed: get only relevant nodes+edges
        if seed:      yr_seed, edges, nodes = self.__queryOnSeed(seed,THRESH)
        else:         yr_seed, edges, nodes = 0,self.edges,self.node_df
        
        fig = go.Figure(layout=dict(
            showlegend=False,
            plot_bgcolor='white',
            yaxis=dict(showticklabels=False)),
            layout_xaxis_range=[1840,1860])
            
        # Draw edges.
        quer1 = (edges.yr1>=MIN_YR)&(edges.yr2>=MIN_YR)
        quer2 = (edges.yr1<=MAX_YR)&(edges.yr2<=MAX_YR)
        plot_edges = edges[quer1&quer2].copy() # Query only relevant ones
        plot_edges["x1"] = nodes.loc[plot_edges.grp1_ind.values,'x'].values  # Map xs for each node
        plot_edges["x2"] = nodes.loc[plot_edges.grp2_ind.values,'x'].values

        for e in range(plot_edges.shape[0]):
            yr1,yr2,x1,x2,k = plot_edges.iloc[e].loc[["yr1","yr2","x1","x2","k"]].values
            
            if yr1!=yr2 and 1-(k-1.2)>THRESH: # plot edges only between different years. and sim>THRESH
                alpha = round((1-(k-1.2)-THRESH)/(1-THRESH),5)
                fig.add_trace(go.Scatter(x=[yr1,yr2],y=[x1,x2],
                                         mode='lines',
                                         marker_color='rgba(0, 0, 0, 0)', # Invisible markers (alpha=0)
                                         line_color=f'rgba(0, 0, 0, {abs(alpha)})')) 


        # Draw nodes
        for yr in range(MIN_YR,MAX_YR+2,2):
        #    ax.axvline(yr,color='k',linewidth=0.5,zorder=-1)      # A line for each year
            fig.add_vline(x=yr,line_width=0.5,line_color='rgba(34, 204, 0, 0.5)')
            
            # Get coords of each group for this year.
            xs = nodes.loc[nodes.yr==yr,"x"].values
            ys = np.repeat(yr,xs.shape[0])
            
            if yr==yr_seed: kwargs=dict(marker_color='rgb(0, 42, 255)',marker_symbol='star',marker_size=10)
            else:           kwargs=dict(marker_color='rgb(220, 0, 0)',)

            #Plot groups for this year, represented each as a point
            fig.add_trace(go.Scatter(x=ys,y=xs,mode='markers',**kwargs,
                                     hovertemplate='%{text}<extra></extra>',
                                     text=nodes.loc[nodes.yr==yr,'elems'].astype(str).str.replace("'|{|}",""),
                                     ids=nodes.loc[nodes.yr==yr,'name'].values))  
    
        fig.layout.clickmode = 'event+select'

        return fig
        
    def __queryOnSeed(self,seed,THRESH): 
        """    Return new graph specification (edges,nodes) 
        containing only edges+nodes connected significantly to nodes specified in `seed`.
        `Significantly` meaning there's a chain of edges with M2>THRESH, 
            connecting an object to any of said nodes.
        All nodes in `seed` should be in the same year (for now at least)
        """
        yr_seed = int(seed[0].split("_")[0]) # Extract year of selection
        keep_nodes = seed   # List to store nodes that "pass" selection

        # First take all edges with relevant M2
        edges = self.edges[self.edges['m2']>THRESH].copy()
        # Take only edges associated with nodes in seed
        relev_ed = edges[(edges.yr_id1.isin(seed) | edges.yr_id2.isin(seed))]

        #Now propagate to the past. 
        #From these nodes, take only those in 2 yrs before, and search only in the past.
        #By construction, yr1<=yr2, so to look in the past you need only look at yr1. (yr2 bzw. for futur)
        p_yr = yr_seed
        past_seed = relev_ed.loc[relev_ed.yr1<p_yr,'yr_id1'].values
        while True:
            p_yr-=2
            keep_nodes += list(past_seed)
            # Find edges participating in relevant edges with nodes in past_seed, in year before
            past_seed = edges[edges.yr_id2.isin(past_seed) & (edges.yr1<p_yr)].yr_id1.unique()  
            if len(past_seed)==0: break

        # Propagate into future, same principle
        f_yr = yr_seed
        futr_seed = relev_ed.loc[relev_ed.yr2>f_yr,'yr_id2'].values
        while True:
            f_yr+=2
            keep_nodes += list(futr_seed)
            futr_seed = edges[edges.yr_id1.isin(futr_seed) & (edges.yr2>f_yr)].yr_id2.unique()  
            if len(futr_seed)==0: break

        # Now query edges, such that both nodes are in keep_nodes
        edges = self.edges[(self.edges.yr_id1.isin(keep_nodes) & 
                            self.edges.yr_id2.isin(keep_nodes) & 
                            self.edges.m2>THRESH)].copy()

        keep_nodes = np.unique(list(edges.yr_id1.values)+ list(edges.yr_id2.values))
        nodes = self.node_df[self.node_df.name.isin(keep_nodes)].copy()
        print(nodes.shape)
        
        return yr_seed,edges, nodes
        


##### END class definition
        
@jit  # Energy function definition
def energy(xs,grp1_ind,grp2_ind,yr1,yr2,z1,z2,ks): 
    """
    Calculate potential energy between nodes. 
    xs: The only variables. Shape is total number of groups. 
    grp{i}_ind: int array for indexing xs. 
        Local Energy calcs (thus opts) are controlled through this and yr{i}
    yr{i}: indexed year."""
    grp1_x = xs[grp1_ind]
    grp2_x = xs[grp2_ind]
    dists = (grp1_x - grp2_x)**2 + (yr1 - yr2)**2 + (z1-z2)**2
    energ = np.sum(2*dists/ks + ks**2/(dists)**2)  
    return energ

# Gradient of energy w.r.t. xs
dEdx = jit(grad(energy))
