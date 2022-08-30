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
from Components.similarity_row import matplot_row
from Components.families_row import families_row
from Components.pt_row import pt_row
from Components.title_row import title, year_slider, footer


#Create the app
server = flask.Flask(__name__)
app = dash.Dash(__name__,
                title="The Evolving Periodic System",
                external_stylesheets = [dbc.themes.BOOTSTRAP],
        server=server)


main_col_plots = dbc.Col(
    [
        matplot_row,
        year_slider,
        pt_row,
        families_row,
    ],
    style={'margin-bottom': '100px'}
)

tabs = dbc.Col(
    [
        title,
        html.Hr(),
        main_col_plots,
        html.Hr(),
        footer
    ],
)

app.layout = html.Div(
    [
        tabs,
        dcc.Store(id='current-perm')
    ]
)


################################## Callbacks ###############################

@app.callback(
        Output('PT-plot','figure'),
        [Input('PT-plot','clickData'), Input('year-slider','value')],
        [State('PT-plot','figure')]
        )
def update_PT(click, yr, fig):
    ctx_trig = dash.callback_context.triggered

    if ctx_trig[0]['prop_id'] in [ 'PT-plot.clickData','year-slider.value' ]:
        try:    # If element exists in clicked position
            s_elem = click['points'][0]
            x, y = s_elem['x'], s_elem['y']
            elem = periodTable.symbol[::-1][y,x]
        except:      elem = None

        fig = go.Figure(fig)
        # Get necessary data for update
        s_data = periodTable.plotSimPT(yr, elem)    

        # Update color pattern
        fig['data'][0]['z'] = s_data[0][::-1]

        # Update labels, to account for elements that don't exist at some given year
        new_annots = [e for e in periodTable.annotations if e.text.split('<')[1][2:] in s_data[1]]
        fig['layout']['annotations'] = new_annots

    return fig


@app.callback(
        Output('current-perm','data'),
        [Input('year-slider','value'),
         Input('opt-button','n_clicks'),
         Input('an-button','n_clicks')],
        [State('current-perm','data')]
        )
def select_perm(year, _, __, current):
    """
    Upon pressing button, select a random index of permutation, and store in dcc.Store
    """
    ctx_trig = dash.callback_context.triggered

    if ctx_trig[0]['prop_id'] == 'opt-button.n_clicks':
        indivs_yr = loadData.opt_permut[2020]
        n = np.random.choice(len(indivs_yr))
        return year, n
    elif ctx_trig[0]['prop_id'] == 'an-button.n_clicks':
        return year, 999
    else:
        return current

@app.callback(
        Output('simmat-plot','figure'),
        [Input('year-slider','value'),
         Input('opt-button','n_clicks'),
         Input('an-button','n_clicks'),
         Input('current-perm','data')],
        [State('simmat-plot','figure')]
        )
def update_simmat(year, _, __, store, fig):
    """
    Update what similarity matrix is being shown.
    User has the option to select a year, and to optimize permutation.
    Depending on the order, various outcomes are possible, e.g.
        Select year, optimize:
            Show simMat for year, with ordering opt for that year.
        Select year1, optimize, select year2:
            Show simMat of year2, with ordering optimized for year1.
    """

    ctx_trig = dash.callback_context.triggered

    # If only year is updated: don't change permutation
    if ctx_trig[0]['prop_id'] == 'year-slider.value':   
        if store:
            prev_yr, perm_n = store
            if perm_n != 999:
                indivs_yr = loadData.opt_permut[prev_yr]
                perm = indivs_yr[perm_n]
            else:
                perm = np.arange(103)
        else:
            perm = np.arange(103)

        fig = simmat.plotSimMat(year, perm)

    # If opt button was activated: change permutation, to this year's
    if ctx_trig[0]['prop_id'] in ['opt-button.n_clicks',
                                  'an-button.n_clicks',]:
        # Select a random pre-optimized permutation for this year. 
        _, perm_n = store

        print(store)
        indivs_yr = loadData.opt_permut[year]
        if type(perm_n)==int:
            if perm_n != 999:
                perm = indivs_yr[perm_n]
            else:
                perm = np.arange(103)
        else:
            perm = np.arange(103)

        fig = simmat.plotSimMat(year, perm)

    return fig



@app.callback(
        Output('show-cost','children'),
        [Input('year-slider','value'),
         Input('opt-button','n_clicks'),
         Input('an-button','n_clicks'),
         Input('current-perm','data')]
        )
def update_cost(year, _, __, store):
    """
    Update box showing the cost of permutation, on the current simmat.
    """

    if store:
        perm_year, perm_n = store
        indivs_yr = loadData.opt_permut[perm_year]

        if type(perm_n)==int:
            if perm_n != 999:
                perm = indivs_yr[perm_n]
            else:
                perm = np.arange(103)
        else:
            perm = np.arange(103)

        S = loadData.simMats[year - 1800].copy()
        P = simmat.symmetrize(S)
        cost = loadData.costPerm(P, perm)

        ctx_trig = dash.callback_context.triggered
        if ctx_trig[0]['prop_id'] in ['year-slider.value',
                                      'opt-button.n_clicks',
                                      'an-button.n_clicks',
                                      ]:
            return [
                html.Div("Sequence optimized in {}".format(perm_year)),
                html.Div("Cost in {} is {:.3f}".format(year, cost))
            ]

    return [""]

# Select elems to modify closure plot
@app.callback(
        Output('closure-plot','figure'),
        [Input('contain-clos','n_submit')],
        [State('contain-clos','value')],
        )
def update_closure(_,incl_ce,):

    ctx_trig = dash.callback_context.triggered[0]["prop_id"]

    if ctx_trig == 'contain-clos.n_submit':
        if incl_ce is not None:
            incl_ce = set(incl_ce.split(","))
    else:
        incl_ce = {'Na'}

    fig = FEs.compMat(
        2021,
        FEs.FEs_df,
        FEs.hist_relev,
        incl_ce = incl_ce,
        update=True,
        show_nr=10,
    )

    return fig


################################## Run app server ###############################
if __name__ == "__main__":
    app.run_server(debug=True,host='0.0.0.0',port=8050)
