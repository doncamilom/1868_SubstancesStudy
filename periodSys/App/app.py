#!/usr/bin/env python3

import dash
from dash import dcc,html
import dash_bootstrap_components as dbc 
from dash.exceptions import PreventUpdate
from dash.dependencies import Input, Output, State, ClientsideFunction

import plotly.express as px
import plotly.graph_objects as go 
import numpy as np

from Components import loadData,simmat, periodTable

#Create the app
app = dash.Dash(__name__,title="Evolution of element families", external_stylesheets = [dbc.themes.BOOTSTRAP])


matbox = html.Div([])

matrix_height = '640px'
matrix_width = '690px'

title = html.Div([html.H1("Visualization of an empirical Periodic System, over history.",
                               style={"textAlign": "center",
                                      "background": "#a6bddb",
                                      "margin-top":"10px"})
                ])

matplot = dbc.Card(
            dbc.CardBody(
                    [dcc.Graph(figure=simmat.fig, id='simmat-plot', 
                                style={'height':matrix_height,'margin-top':'0px'})]
                )
        )

otherplot = dbc.Card(
            dbc.CardBody(
                    [dcc.Graph(figure=periodTable.fig, id="PT-plot",
                        style=dict(height='700px'))]
                )
        )

row = dbc.Row([matplot, otherplot])
tabs = dbc.Col([title,  row])



app.layout = html.Div([tabs])
    

###################################################### Callbacks ###################################################



###################################################### Run app server ###################################################    
    
if __name__ == "__main__":
    app.run_server(debug=True,host='localhost',port='8080')
