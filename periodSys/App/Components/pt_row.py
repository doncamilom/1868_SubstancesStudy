#!/usr/bin/env python3

import dash
from dash import dcc,html
import dash_bootstrap_components as dbc 
from dash.exceptions import PreventUpdate
from dash.dependencies import Input, Output, State, ClientsideFunction

import plotly.express as px
import plotly.graph_objects as go 
import numpy as np
import flask

from Components import loadData,simmat, periodTable, FEs


# Periodic Table
PTplot = dbc.Col(
    [
        dcc.Graph(figure=periodTable.fig,
                  id="PT-plot",
                  style={"width":'100%',
                         "margin-top":"30px"})
    ],
    align='center',
    style={'width': '100%',
           'margin-top': '30px'},
    width={"size": 7,
           "order": 2,
           "offset": 0}
)

# PT description
pt_descr = dbc.Col(
    [
        html.Div(
            [
                html.H1("The Periodic System encodes similarity."),
            ],
            style={'height': '100px'}
        ),
        html.P("Elements in the same column tend to be similar, \
        however this is not always true!",
               className="p-html-text",
               style={"margin-left": "30px",
                      "margin-top": "30px"}),
        html.P("Select an element in the Periodic Table,\
        to see what are the most similar elements to it.",
               className="p-html-instruction"),
        html.P("Note how this is time dependent.",
               className="p-html-highlight")
    ],
    width={"size": 4,
           "offset": 1,
           "order": 1},
    align="center"
)

# Build row
pt_row = html.Div(
    dbc.Row(
        [
            PTplot,
            pt_descr
        ],
        justify='end',
        style={'width': '100%',
               'margin-top': '20px',}
    ),
    style={"margin-left": '0px',
           "margin-bottom": '190px'}
)
