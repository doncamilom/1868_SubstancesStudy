#!/usr/bin/env python3

import dash
from dash import dcc,html
#import dash_bootstrap_components as dbc 
from dash.exceptions import PreventUpdate
from dash.dependencies import Input, Output, State, ClientsideFunction

import plotly.express as px
import plotly.graph_objects as go 
import numpy as np

from Components import loadData,simmat

#Create the app
app = dash.Dash(__name__,title="Evolution of element families")


matbox = html.Div([])

matrix_height = '650px'
matrix_width = '720px'
app.layout = html.Div([ html.H1("Visualization of an empirical Periodic System, over history.",
                               style={"textAlign": "center",
                                      "background": "#a6bddb"}),
                        html.Div([dcc.Graph(figure=simmat.fig, id='simmat-plot', 
                                style={'height':matrix_height,'margin-top':'0px'})],
                        style={'height':matrix_height,'width':matrix_width})
                         ])#,
                        #style={'height':'14000px','width':'15000px'})
    

###################################################### Callbacks ###################################################



###################################################### Run app server ###################################################    
    
if __name__ == "__main__":
    app.run_server(debug=True,host='localhost',port='8080')
