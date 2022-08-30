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


# Similarity matrix
matrix_width = '530px'

matrix_fig = dbc.Col(
    [
        dcc.Graph(figure=simmat.fig,
                  id='simmat-plot',
                  style={'height':'100%',
                         'margin-left':'10%',
                         'margin-right':'10%'})
    ],
    align="center"
)


tab_width = "550px"

optimize_col_div = html.Div(
    [
        html.H3("Optimize the sequence of elements"),
        html.Div(
            [
                html.P("Similarity information can be encoded in a sequence \
                of elements, so that similar elements are closer together.",
                       className="p-html-text"),
                html.P("Using genetic algorithms,\
                we find such optimal sequences.",
                       className="p-html-text"),
            ],
            style={'margin-left': '40px',
                   'margin-right': '70px',
                   "margin-top": "20px",
                   "margin-bottom": "10px"}
        ),
        html.Div(
            [
                html.Div("Press the button to optimize the sequence.",
                         className="p-html-instruction"),
                html.Div("This will bring high values of the matrix \
                closer to the diagonal.",
                        className="p-html-text"),
            ],
                    style={"width": tab_width,
                           "margin-top": "10px",
                           "margin-bottom": "20px"
                           }
        ),
        html.Div(
            [
                html.Button("Optimize sequence",
                            id='opt-button',
                            style={'text-align': 'center',
                                   'margin-bottom': '10px'
                                   }
                            ),
                html.Button("Atomic Number",
                            id='an-button',
                            style={'text-align': 'center',
                                   'margin-bottom': '10px'
                                   }
                            ),
            ],
            style={'text-align': 'center'}
        ),
        html.Div([],
                 id="show-cost",
                 style={'height':'90px',
                        'text-align':'center',
                        'font-size':20,
                        'width':'100%',
                        'position': 'relative'},
                 className="p-html-highlight"
                 )
    ],
    style={"width": tab_width}
)

simmat_descript = html.Div(
    [
        html.H3("Similarity between the chemical elements"),
        html.Div(
            [
                html.P("Chemical elements show resemblances to others in \
                the compounds they form.\
                Two elements are alike if they both form compounds \
                with similar compositions.",
                       className="p-html-text"),
                html.P("This similarity matrix encodes how similar each element \
                is to any other one.",
                       className="p-html-text"),
            ],
            style={'margin-left': '40px',
                   'margin-right': '70px',
                   "margin-top": "20px",
                   "margin-bottom": "30px"}
        ),

        html.P("Hover over any pixel to visualize the similarity \
        between a pair of elements",
                className="p-html-instruction"),
        html.Div("",
       #     display_hover,
            style={'margin-left': '40px',
                   'margin-right': '70px',
                   'height': '100px',
                   'text-align': 'center'}
        ),
    ],
    style={"width": tab_width,}
)


# Description and guide to similarity matrix
simmat_p = html.Div(
    dbc.Col(
        dcc.Tabs(
            [
                dcc.Tab(simmat_descript,
                        label="Similarity",
                        className='custom-tab',
                        selected_className='custom-tab--selected'
                        ),
                dcc.Tab(optimize_col_div,
                        label="Optimization",
                        className='custom-tab',
                        selected_className='custom-tab--selected'
                        ),
            ],
            parent_className='custom-tabs',
            className='custom-tabs-container',
        ),
        align='center',
        style={'margin-left':'0px'}
    ),
)


# Build similarity row
matplot_row = dbc.Row(
    [
        matrix_fig,
        simmat_p
    ],
    justify='between',
    style={'width': '100%',
           'margin-top': '60px',
           }
)
