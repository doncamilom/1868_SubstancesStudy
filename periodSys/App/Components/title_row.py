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


title_descrip = dbc.Col(
    [
        html.Div(
            [
                html.H1("The Evolving Periodic System",
                        style={
                            "margin-top":"20px",
                            'width':'100%'}
                        ),
                html.P("Chemistry evolves rapidly over the years, but...\
                how does the Periodic System react to these changes?",
                       style={"width":"70%",
                              "font-size":"110%",
                            "margin-left":'0'},
                       )
            ],
        )
    ],
    width=7,
    align='start'
)


mpi_title_img = dbc.Col(
    [
        html.Img(src="./assets/ENG_MPI_MiS_Bildwortmarke_farbig.png",
                 style={
                     'width':'100%',
                     'float':'right',
                     'margin-top':'0',
                     'padding-top':'0',
                     'padding-right':'0'
                 },
                 ),
    ],
    width=5,
    align='start'
)

title = html.Div(
    [
        dbc.Row(
            [
                title_descrip,
                mpi_title_img
            ],
            className='g-0',
            justify='between',
            style={'height':'140px'}
        )
    ],
)

year_slider = html.Div(
    [
        html.H1("Change the year to see how similarities are modified",
                style={'text-align':'center',
                       'margin-top':'80px',
                       'font-size':'200%',
                       }
        ),
        dcc.Slider(id="year-slider",
                   min=1800,
                   max=2021,
                   step=1,
                   value=2021,
                   marks={yr:str(yr) for yr in range(1800,2022,10)},
                   tooltip={"placement": "top",
                            "always_visible": False}),
    ],
    style={'margin-bottom':'20px',
           'margin-top':'30px',
           'margin-right':'140px',
           'margin-left':'140px',
           }
)


footer = html.P(
    [
        html.A("Source code", href="https://github.com/doncamilom/Reaxys_PS/tree/master/periodSys/App"),
        ". Read the ",
        html.A("paper", href="https://arxiv.org/"),
        ". App built by ",
        html.A("Andres M Bran ", href="https://twitter.com/drecmb"),
        ":)"
    ],
)
