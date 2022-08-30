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
        html.A(href = "https://www.mis.mpg.de/",
               target = "_blank",
               children = [
                   html.Img(src="./assets/ENG_MPI_MiS_Bildwortmarke_farbig.png",
                            style={
                                'width':'100%',
                                'float':'right',
                                'margin-top':'0',
                                'padding-top':'0',
                                'padding-right':'0'
                            },
                            ),
               ]
               )
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


footer_left = dbc.Col(
    [
        html.A("Source code",
               href="https://github.com/doncamilom/Reaxys_PS/tree/master/periodSys/App",
               target="_blank"),
        ". Read the ",
        html.A("paper",
               href="https://arxiv.org/",
               target="_blank"),
        ". App built by ",
        html.A("Andres M Bran ",
               href="https://twitter.com/drecmb",
               target = "_blank"),
        ":)"
    ],
)

footer_right = dbc.Row(
    [
        html.Div("Data analysis with data from ",
                 style={"width": "230px"}),
        html.A(href = "https://www.reaxys.com",
               target = "_blank",
               children = [
                   html.Img(src="./assets/reaxys.png",
                            style={
                                'width':'100%',
                                'float':'right',
                                'margin-bottom':'10px',
                                'margin-left': '60px',
                            },
                            ),
               ],
            style={"width": "20%"}
        )
    ],
    style={"width": "650px",
           'margin-right': '60px',
           },
    #width=5,
    justify="end"
    #align='end'
)

footer = dbc.Row(
    [
        footer_left,
        footer_right
    ]
)
