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

#Create the app
server = flask.Flask(__name__)
app = dash.Dash(__name__,
                title="The Evolving Periodic System",
                external_stylesheets = [dbc.themes.BOOTSTRAP],
        server=server)


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
        html.H1("Select the year you want to explore.",
                style={'text-align':'center',
                       'margin-top':'80px',
                       'font-size':'200%'}
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
    style={'margin-bottom':'60px',
           'margin-top':'30px',
           'margin-right':'140px',
           'margin-left':'140px',
           }
)


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
    align="center")


display_hover = [
    dbc.CardHeader("Display on hover"),
    dbc.CardBody(
        [
            html.H5("Display stuff", className="card-title"),
            html.P(
                "This is some card content that we'll reuse",
                className="card-text",
            ),
        ]
    ),
]



tab_width = "550px"
optimize_col_div = html.Div(
    [
        html.H3("Optimize the sequence of elements"),
        html.Div(
            [
                html.Div("Press the button to optimize the sequence.",
                         className="p-html-instruction"),
                html.Div("This will bring high values of the matrix \
                closer to the diagonal.",
                        className="p-html-text"),
            ],
                    style={"width": tab_width,
                           "margin-top": "50px",
                           "margin-bottom": "30px"
                           }
        ),
        html.Div(
            [
                html.Button("Optimize sequence",
                            id='opt-button',
                            style={'text-align': 'center',
                                   'margin-bottom': '50px'
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
    ]
)

simmat_descript = html.Div(
    [
        html.H3("Similarity between the chemical elements."),
        html.Div(
            [
                html.P("Chemical elements show resemblances to others\
                in the way they react, and in the compounds they form.",
                        className="p-html-text"),
                html.P("This similarity matrix encodes how similar each element \
                is to any other one.",
                    className="p-html-text"),
            ],
            style={'margin-left': '40px',
                'margin-right': '70px'}
        ),

        html.P("Hover over any pixel to visualize the similarity \
        between a pair of elements",
                className="p-html-instruction"),
        html.Div(
            display_hover,
            style={'margin-left': '40px',
                'margin-right': '70px',
                   'height': '150px',
                   'text-align': 'center'}
        ),
    ],
    style={"width": tab_width}
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


# Families evolution plot
family_plot = dbc.Col(
    [
        dcc.Graph(figure=FEs.closure_fig,
                  id="closure-plot",
                  style={"width":"100%",
                         "height": "250px",
                         "margin-top":"0px",
                         }
                  )
    ],
    width={"size": 8,
           "offset": 0,
           "order": "last"
    }
)

# Descripting text
families_p = dbc.Col(
    [
        html.Div(
            [
                html.H1("Evolution of families"),
                html.P("Similar elements are clustered into families.",
                    className='p-html-highlight'),
                html.P("Select an element to visualize its evolution in families.",
                    className='p-html-instruction'),
                html.Div(dcc.Input(id="contain-clos",
                                type="text",
                                placeholder="Input elements:"
                                ),
                        style={#'height': '30px',
                            'width': '250px',
                        }
                        )
            ],
            style={'height': '300px'}
        )
    ],
    width={"size": 4,
           "offset": 0,
           "order": 1}
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
               'margin-top': '120px',}
    ),
    style={"margin-left": '60px'}
)



# Periodic Table
PTplot = dbc.Row(
    [
        dcc.Graph(figure=periodTable.fig,
                  id="PT-plot",
                  style={"width":'750px',
                         "margin-top":"30px"})
    ],
    justify='center',
    style={'width': '100%',
           'margin-top': '130px'}
)


main_col_plots = dbc.Col(
    [
        matplot_row,
        families_row,
        PTplot
    ]
)


tabs = dbc.Col(
    [
        title,
        year_slider,
        main_col_plots
    ],
    style={'margin-bottom':'100px'}
)

app.layout = html.Div(
    [
        tabs,
        dcc.Store(id='current-perm')
    ]
)

###################################################### Callbacks ###################################################

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
        s_data = periodTable.plotSimPT(yr, elem)    # Get necessary data for update

        # Update color pattern
        fig['data'][0]['z'] = s_data[0][::-1]

        # Update labels, to account for elements that don't exist at some given year
        new_annots = [e for e in periodTable.annotations if e.text.split('<')[1][2:] in s_data[1]]
        fig['layout']['annotations'] = new_annots

    return fig


@app.callback(
        Output('current-perm','data'),
        [Input('year-slider','value'), Input('opt-button','n_clicks')],
        [State('current-perm','data')]
        )
def select_perm(year, _, current):
    """
    Upon pressing button, select a random index of permutation, and store in dcc.Store
    """
    ctx_trig = dash.callback_context.triggered

    if ctx_trig[0]['prop_id'] == 'opt-button.n_clicks':
        indivs_yr = loadData.opt_permut[2020]
        n = np.random.choice(len(indivs_yr))
        return year, n
    else: return current


@app.callback(
        Output('simmat-plot','figure'),
        [Input('year-slider','value'),Input('opt-button','n_clicks'), Input('current-perm','data')],
        [State('simmat-plot','figure')]
        )
def update_simmat(year, _, store, fig):
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

    if ctx_trig[0]['prop_id'] == 'year-slider.value':   # If only year is updated: don't change permutation
        if store:
            prev_yr, perm_n = store
            indivs_yr = loadData.opt_permut[prev_yr]
            perm = indivs_yr[ perm_n ]
        else:
            perm = np.arange(103)

        fig = simmat.plotSimMat(year, perm)

    if ctx_trig[0]['prop_id'] == 'opt-button.n_clicks': # If opt button was activated: change permutation, to this year's
        # Select a random pre-optimized permutation for this year. 
        _, perm_n = store

        indivs_yr = loadData.opt_permut[year]
        if type(perm_n)==int:   perm = indivs_yr[ perm_n ]
        else:                   perm = range(103)

        fig = simmat.plotSimMat(year, perm)

    return fig



@app.callback(
        Output('show-cost','children'),
        [Input('year-slider','value'), Input('opt-button','n_clicks'), Input('current-perm','data')]
        )
def update_cost(year, _, store):
    """
    Update box showing the cost of permutation, on the current simmat.
    """

    if store:
        perm_year, perm_n = store
        indivs_yr = loadData.opt_permut[perm_year]

        if type(perm_n)==int:   perm = indivs_yr[ perm_n ]
        else:                   perm = np.arange(103)

        S = loadData.simMats[year - 1800].copy()
        P = simmat.symmetrize(S)
        cost = loadData.costPerm(P,perm)

        ctx_trig = dash.callback_context.triggered
        if ctx_trig[0]['prop_id'] in ['year-slider.value', 'opt-button.n_clicks']:
            return [
                html.Div("Sequence optimized in {}".format(perm_year)),
                html.Div("Cost in {} is {:.3f}".format(year, cost))
            ]

    return ["Void container"]



# Select elems to modify closure plot
@app.callback(
        Output('closure-plot','figure'),
        [Input('contain-clos','n_submit')],
        [State('contain-clos','value')],
        )
def update_closure(_,incl_ce,):

    ctx_trig = dash.callback_context.triggered[0]["prop_id"]

    valid_trigs = ['contain-clos.n_submit','non-contain-clos.n_submit'] 

    if ctx_trig in valid_trigs:
        if incl_ce is not None:
            incl_ce = set(incl_ce.split(","))

        #if notincl_ce is not None:
            #notincl_ce = set(notincl_ce.split(","))
    else:
        incl_ce = {'Na'}


    print(incl_ce)#, notincl_ce)
    fig = FEs.HeatmapClosure(2021,
                       FEs.FEs_df,FEs.hist_relev,
                       incl_ce = incl_ce,
                  #     notincl_ce = notincl_ce,
                       thresh_relev=0.)

    return fig


###################################################### Run app server ###################################################    
    
if __name__ == "__main__":
    app.run_server(debug=True,host='0.0.0.0',port=8050)
