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



# Families evolution plot
family_plot = dbc.Col(
    [
        dcc.Graph(figure=FEs.closure_fig,
                  id="closure-plot",
                  style={"width":"100%",
                         "height": "300px",
                         "margin-top":"0px",
                         }
                  )
    ],
    width={"size": 8,
           "offset": 0,
           "order": 1
           },
    align="center"

)

# Descripting text
families_p = dbc.Col(
    [
        html.Div(
            [
                html.H1("Evolution of families"),
                html.P("Similar elements are clustered into families",
                       className='p-html-text',
                       style={"margin-left": "30px"}),
                html.P("Enter a (list of) elements to visualize its evolution in families",
                    className='p-html-instruction'),
                html.Div(dcc.Input(id="contain-clos",
                                   type="text",
                                   placeholder="Input elements:"
                                   ),
                         style={#'height': '30px',
                            'width': '250px',
                        }
                         ),
                html.Div("Each row shows the evolution of a single family.",
                         className='p-html-text'),
                html.Div("Here, black means the family was found exactly \
                in that year, while red, that the family was found in that \
                year as a subset of another family.",
                         className='p-html-highlight',
                         style={"margin-left": "30px"}),
                html.Div("Hover over it to see more information.",
                         className='p-html-instruction'),
            ],
            style={'height': '400px'}
        )
    ],
    width={"size": 4,
           "offset": 0,
           "order": "last"}
)

# Build row
families_row = html.Div(
    dbc.Row(
        [
            families_p,
            family_plot
        ],
        justify='end',
        style={'width': '100%',
               'margin-top': '50px',}
    ),
    style={"margin-left": '60px',
           "margin-bottom": '190px'}
)
