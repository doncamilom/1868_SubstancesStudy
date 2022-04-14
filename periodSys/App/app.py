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
                    dcc.Slider(id="year-slider",min=1800,max=2021,step=1,value=2021,
                            marks={yr:str(yr) for yr in range(1800,2022,10)},
                            tooltip={"placement": "top", "always_visible": True}),
                    ],
                    style={'margin-bottom':'40px'})

matplot = dcc.Graph(figure=simmat.fig, id='simmat-plot', 
        style={'height':matrix_height,'width':matrix_width,'margin-top':'0px'})


PTplot = dcc.Graph(figure=periodTable.fig, id="PT-plot",
                style=dict(height='300px'))

col2 = dbc.Col([
            dbc.Row([
                    html.Button("Optimize", id='opt-button',),
                    html.Div("Cost: -3.6",id = 'test-box2',
                            style={'height':'200px','text-align':'center',
                                    'font-size':25,'width':'200px'}),
                    ]),
            PTplot
            ])

row = dbc.Row([matplot, col2])
tabs = dbc.Col([title, year_slider,  row],
        style={'margin-bottom':'100px'})



app.layout = html.Div([tabs])
    

###################################################### Callbacks ###################################################

@app.callback(
        Output('PT-plot','figure'),
        [Input('PT-plot','clickData'), Input('year-slider','value')],
        [State('PT-plot','figure')]
        )
def test_pt(inp, yr, fig):
    ctx_trig = dash.callback_context.triggered
    if ctx_trig:
        try:
            s_elem = inp['points'][0]
            x,y = s_elem['x'], s_elem['y']
            elem = periodTable.symbol[::-1][y,x]
        except:      elem = None

        fig = go.Figure(fig)
        s_data = periodTable.plotSimPT(yr, elem)    # Get necessary data for update

        # Update color pattern
        fig['data'][0]['z'] = s_data[0][::-1]

        # Update labels, to account for elements that don't exist at some given year
        new_annots = [e for e in periodTable.annotations if e.text.split('<')[1][2:] in s_data[1]]
        fig['layout']['annotations'] = new_annots

    return fig


@app.callback(
        Output('simmat-plot','figure'),
        [Input('year-slider','value'),Input('opt-button','n_clicks')],
        [State('simmat-plot','figure')]
        )
def update_simmat(year, _, fig):
    global perm
    ctx_trig = dash.callback_context.triggered

    if ctx_trig[0]['prop_id'] == 'year-slider.value':
        fig = simmat.plotSimMat(year, perm)

    if ctx_trig[0]['prop_id'] == 'opt-button.n_clicks':
        # Select a random pre-optimized permutation for this year. 
        indivs_yr = loadData.opt_permut[year]
        perm = indivs_yr[ np.random.choice(len(indivs_yr)) ]

        print(perm)


        fig = simmat.plotSimMat(year, perm)

    return fig

#ctx_trig = dash.callback_context.triggered
#print(ctx_trig)
#if ctx_trig[0]['prop_id'] == 'year-slider.value':
#Input('simmat-plot','clickData'), 
#[State('simmat-plot','figure')]

###################################################### Run app server ###################################################    
    
if __name__ == "__main__":
    app.run_server(debug=True,host='localhost',port='8080')
