
import plotly.graph_objects as go
import numpy as np
from Components import simmat

elemList = simmat.elemList

# Periodic Table Data
import Components.TPs as TPs

# Change this to change the TP to be used
TP = TPs.TPshort

symbol = np.zeros(TP['shape'], 'U4')
for elem in TP['xy'].keys():
    x,y = TP['xy'][elem]
    symbol[x,y] = elem

element = np.vectorize(TP['names'].get)(symbol)
atomic_mass = np.zeros(TP['shape'])
color = np.zeros(TP['shape'])
for elem in TP['xy'].keys():
    x,y = TP['xy'][elem]
    color[x,y] = 1.0



P = simmat.P
def plotSimPT(year,elem):
    """
    Get values of similarity between selected element, and all other elements.
    Format this in a periodic table as specified in TP module.
    """

    Pyr = P[year-1800].copy()

    # Get list of elements that exist this year
    isna = np.isnan(np.diag(Pyr))
    elems_yr = np.array(simmat.elemList)[~isna]

    if elem==None:
        return color, elems_yr


    if elem != '':
        X, Y = TP['xy'][elem]
    else:
        X, Y = 0,0

    # Get list of similarities between our element and all the rest
    sims = np.nan_to_num(Pyr[elemList.index(elem)],0.)
    if sims.sum()==0.:
        sims += 1.    # If all 0., click stops responding

    img = np.zeros(TP['shape'])
    for e in simmat.elemList:
        x, y = TP['xy'][e]
        img[x, y] = sims[elemList.index(e)]

    # Modify entry for our element, so others are not opaque
    img[X, Y] = 0
    img[X, Y] = np.nanmax(img)
    img /= np.nanmax(img)

    return img,elems_yr



# Display element name and atomic mass on hover
hover=[]
for x in range(len(symbol)):
    hover.append([i + '<br>' + 'Atomic Mass: ' + str(j) if i else ''
                      for i, j in zip(element[x], atomic_mass[x])])

# Set Colorscale
colorscale=[
        [0.0, 'rgba(0,0,0,0)'],
        [1.0, 'rgba(49,108,159,0.6)']
        ]


fig = go.Figure()
fig.add_trace(go.Heatmap(
                z=color[::-1],
                colorscale=colorscale,
                #text=hover[::-1],
                hoverinfo='text',
                showscale=False,
                ygap=5,
                xgap=5
                
    ))

annotations = []
for n, row in enumerate(symbol):
    for m, val in enumerate(row):
        annotations.append(go.layout.Annotation(text="<b>{}</b>".format(symbol[::-1][n][m]), x=m, y=n,
                                         xref='x1', yref='y1', showarrow=False))

fig.update_layout(
        annotations=annotations,
        #title='Periodic Table',
        margin=dict(l=40, r=40, t=50, b=0, pad=0),
        xaxis=dict(zeroline=False, showgrid=False,visible=False),
        yaxis=dict(zeroline=False, showgrid=False, scaleanchor="x", visible=False),
        plot_bgcolor='white',
        autosize=True,
)
