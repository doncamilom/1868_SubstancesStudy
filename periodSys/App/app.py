#!/usr/bin/env python3

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
import numpy as np

#Create the app
app = dash.Dash(__name__,title="Evolution of element families")

app.layout = html.Div([html.H1("Visualization of an empirical Periodic System through history.",
                               style={"textAlign": "center",
                                      "background": "yellow"}),
    
                    ])

###################################################### Callbacks ###################################################



###################################################### Run app server ###################################################    
    
if __name__ == "__main__":
    app.run_server(debug=True,host='localhost')
