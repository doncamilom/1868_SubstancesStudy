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
OriginalFig.layout.clickmode = 'event+select'

#Create the app
server = flask.Flask(__name__)
app = dash.Dash(__name__)#, external_stylesheets = [dbc.themes.BOOTSTRAP],server=server) #USING BOOTSTRAP'S CSS LIBRARY

app.layout = html.Div([html.H1("Evolution of element families",
                               style={
                                      "textAlign": "center",
                                      "background": "yellow"}),
                       dcc.Graph(figure=OriginalFig,id="mainGraph"),
                        html.Div([html.Button("Apply selection",id="apply-button")]),
                        html.Div([html.Button("Reset plot",id="reset-button")])

                      ], #style={
                          #      "background": "#000080"}
                         )


###################################################### Callbacks ###################################################


@app.callback(
    Output("mainGraph","figure"),
    [Input("apply-button","n_clicks"),Input("reset-button","n_clicks")],
    [State("mainGraph","selectedData"),State("mainGraph","figure")],
)
def fig_updater(applyB,resetB,sel_data,fig):
    global OriginalFig
    ctx_trig = dash.callback_context.triggered[0]
    if ctx_trig['prop_id']!='.':

        print(ctx_trig)

        if ctx_trig["prop_id"]=="apply-button.n_clicks":
            # Extract information about selected points from `sel_data`
            seed = []
            for i in sel_data['points']:
                seed.append(f"{i['x']}_{i['pointIndex']}")
        
            fig = graph.plotGraph(1840,1860,THRESH=0.,seed = seed)

        elif ctx_trig["prop_id"]=="reset-button.n_clicks":
            fig = OriginalFig

    return fig

# TO-DO: Add a slider for changing THRESH.


###################################################### Run app server ###################################################    
    
if __name__ == "__main__":
    app.run_server(debug=True,host='localhost')
