#!/usr/bin/env python3

import flask
import dash
from dash import dcc,html
#import dash_bootstrap_components as dbc 
from dash.exceptions import PreventUpdate
from dash.dependencies import Input, Output, State, ClientsideFunction

import plotly.express as px
import plotly.graph_objects as go 
import base64
import os
import io

from Components import graph,loadData
import numpy as np


grps = loadData.grps

graph = graph.HistGraph(grps)
OriginalFig = graph.plotGraph(1840,1860,THRESH=0.)

#Create the app
server = flask.Flask(__name__)
app = dash.Dash(__name__,title="Evolution of element families")

app.layout = html.Div([html.H1("Evolution of element families",
                               style={"textAlign": "center",
                                      "background": "yellow"}),

                       html.Div(dcc.RangeSlider(id="year-slider",min=1800,max=2016,step=2,value=[1840,1860],
                                        tooltip={"placement": "bottom", "always_visible": True})),
                        dcc.Graph(figure=OriginalFig,id="mainGraph",style={"height":"400px","margin-top":"10px"}),

                        html.Div([html.Button("Apply selection",id="apply-button"),
                                  html.Button("Reset plot",id="reset-button"),
                                  html.Button("Optimize layout",id="opt-button"),
                                 ],style={"text-align":"center"}
                                ),

                        html.Div([html.Div("Select Threshold value",
                                            style={"padding": "15px 0px 0px 0px", # Add a little space above
                                                   "text-align":"center"}),

                                  html.Div(dcc.Slider(id="Thresh-slider",min=0.,max=1.,step=0.05,value=0.,
                                        tooltip={"placement": "bottom", "always_visible": True}),

                                            style={"padding": "0px 400px 0px 400px"} # Center
                                            )
                                ])
                      ])

###################################################### Callbacks ###################################################

seed = False      # If seed is not defined, use everything

@app.callback(
    Output("mainGraph","figure"),
    [Input("year-slider","value"),Input("apply-button","n_clicks"),Input("reset-button","n_clicks"),
     Input("opt-button","n_clicks"),Input("Thresh-slider","value")],
    [State("mainGraph","selectedData"),State("mainGraph","figure")],
)
def fig_updater(year_range,applyB,resetB,optB,THRESH,sel_data,fig):
    global OriginalFig,seed

    ctx_trig = dash.callback_context.triggered[0]

    if ctx_trig['prop_id']!='.':
        minyr,maxyr=year_range

        if ctx_trig["prop_id"]=="apply-button.n_clicks":
            seed = [grp['id'] for grp in sel_data['points']] # Extract information about selected points from `sel_data`
            fig = graph.plotGraph(minyr,maxyr,THRESH=THRESH,seed = seed)

        elif ctx_trig["prop_id"]=="reset-button.n_clicks":
            seed = False
            fig = OriginalFig

        elif ctx_trig["prop_id"]=="Thresh-slider.value":
            fig = graph.plotGraph(minyr,maxyr,THRESH=THRESH,seed = seed)

        elif ctx_trig["prop_id"]=="opt-button.n_clicks":
            yr_seed = int(seed[0].split("_")[0])
            graph.Optimize(yr_seed-4,yr_seed+4,seed=seed,THRESH=THRESH,EPOCHS=5000,lr=1e-1,v=0)
            fig = graph.plotGraph(minyr,maxyr,THRESH=THRESH,seed=seed)

        elif ctx_trig["prop_id"]=="year-slider.value":
            fig = graph.plotGraph(minyr,maxyr,THRESH=THRESH,seed = seed)

    return fig


@app.callback(
    Output("Thresh-slider","value"),
    [Input("reset-button","n_clicks")],
)
def reset_slider(butt):
    return 0.
    

# TO-DO: Add a slider for changing THRESH.


###################################################### Run app server ###################################################    
    
if __name__ == "__main__":
    app.run_server(debug=True,host='localhost')
