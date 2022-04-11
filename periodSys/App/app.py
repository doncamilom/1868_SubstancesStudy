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


title = html.Div([html.H1("Visualization of an empirical Periodic System, over history.",
                               style={"textAlign": "center",
                                      "background": "#a6bddb",
                                      "margin-top":"10px"})
                ])


matrix_height = '600px'
matrix_width = '670px'

year_slider = html.Div([
                    #html.H1("Year",
                    #        style={'text-align':'center',
                    #                'font-size':'1'}),
                    dcc.Slider(id="year-slider",min=1800,max=2022,step=1,value=2022,
                            marks={yr:str(yr) for yr in range(1800,2022,10)},
                            tooltip={"placement": "top", "always_visible": True}),
                    ],
                    style={'margin-bottom':'40px'})

matplot = dcc.Graph(figure=simmat.fig, id='simmat-plot', 
        style={'height':matrix_height,'width':matrix_width,'margin-top':'0px'})


otherplot = dcc.Graph(figure=periodTable.fig, id="PT-plot",
                style=dict(height='400px'))

col2 = dbc.Col([
            dbc.Row([
                    html.Div("Some information box :)",id = 'test-box',
                            style=dict(height='200px')),
                    html.Div("Some information box 2 :)",id = 'test-box2',
                            style={'height':'200px','text-align':'center'}),
                    ]),
            otherplot
            ])

row = dbc.Row([matplot, col2])
tabs = dbc.Col([title, year_slider,  row],
        style={'margin-bottom':'100px'})



app.layout = html.Div([tabs])
    

###################################################### Callbacks ###################################################

@app.callback(
        Output('test-box','children'),
        [Input('PT-plot','clickData')]
        )
def test_pt(inp):
    ctx_trig = dash.callback_context.triggered
    print(ctx_trig )

    return "Update"


@app.callback(
        Output('simmat-plot','figure'),
        [Input('simmat-plot','clickData'), Input('year-slider','value')],
        [State('simmat-plot','figure')]
        )
def test_simmat(inp,year,fig):
    ctx_trig = dash.callback_context.triggered
    print(ctx_trig)

    #if ctx_trig[0]['prop_id'] == 'year-slider.value':
    fig = simmat.plotSimMat(year)

    return fig

###################################################### Run app server ###################################################    
    
if __name__ == "__main__":
    app.run_server(debug=True,host='localhost',port='8080')
