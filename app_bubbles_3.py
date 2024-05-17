# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash, dcc, html
import plotly.express as px
import pandas as pd
import plotly.graph_objs as go
import textwrap
import numpy as np


import dash_bootstrap_components as dbc

app = Dash(external_stylesheets=[dbc.themes.COSMO])

df = pd.read_csv('data7.csv')

df['hit_descriptions'] = df['hit_descriptions'].apply(lambda t: "<br>&nbsp;".join(textwrap.wrap(t, 75)))

#df['-log10HC_vs_VLC_p.adj']=-1*(np.log10(df['HC_vs_VLC_p.adj']))
#df['-HC_vs_VLC_log2 fold change']=-1*(df['HC_vs_VLC_log2 fold change'])
df['VLC_avg']=df[['VLC1','VLC2','VLC3']].mean(axis=1)
#df['VLC_SDperc']=(df[['VLC1','VLC2','VLC3']].std(axis=1))/df['VLC_avg']

#df= df.loc[df['RbcL_vs_BST2_log2 fold change']>0]

fig = px.scatter(df, x="HC_vs_VLC_log2 fold change", y="Log2FC(10030/5per)",
                 size="VLC_avg", color="SAGA_vs_BST2_log2FC", hover_data=["assignment", "Prediction", "description_x"],
                 log_y=False, size_max=100, custom_data=['identifier', 'assignment', 'Prediction', 'description_x', 'transmembranes', 'tm_likelihood', 'VLC_avg', 'number_motifs', 'hit_descriptions'])

fig.update_layout(hovermode='closest', legend=dict(title= None), hoverlabel=dict(bgcolor='rgba(255,255,255,0.75)',font=dict(color='black')))

buttonlist_x = []
for col in df.columns:
    buttonlist_x.append(
        dict(
            args=["x",[df[str(col)]] ],
            label=str(col),
            method='restyle'
        )
    )

buttonlist_y = []
for col in df.columns:
    buttonlist_y.append(
        dict(
            args=['y',[df[str(col)]] ],
            label=str(col),
            method='restyle'
        )
    )

buttonlist_size = []
for col in df.columns:
    buttonlist_size.append(
        dict(
            args=['size',[df[str(col)]] ],
            label=str(col),
            method='restyle'
        )
    )

buttonlist_x_scale = (
    dict(
        args=[{'xaxis.type': 'linear'}],
        label="Linear",
        method="relayout"
    ),
    dict(
        args=[{'xaxis.type': 'log'}],
        label="Log<sub>10</sub>",
        method="relayout"
    )
)

buttonlist_y_scale = (
    dict(
        args=[{'yaxis.type': 'linear'}],
        label="Linear",
        method="relayout"
    ),
    dict(
        args=[{'yaxis.type': 'log'}],
        label="Log<sub>10</sub>",
        method="relayout"
    )
)


fig.update_layout(
        updatemenus=[
            go.layout.Updatemenu(
                buttons=buttonlist_x,
                direction="down",
                #pad={"r": 10, "t": 10},
                showactive=True,
                x=0,
                xanchor="left",
                y=1.12,
                yanchor="top"
            ),
            go.layout.Updatemenu(
                buttons=buttonlist_y,
                direction="down",
                #pad={"r": 10, "t": 10},
                showactive=True,
                x=0.25,
                xanchor="left",
                y=1.12,
                yanchor="top"
            ),
            go.layout.Updatemenu(
                buttons=buttonlist_size,
                direction="left",
                #pad={"r": 10, "t": 10},
                showactive=True,
                x=0.5,
                xanchor="center",
                y=1.12,
                yanchor="top"
            ),
            go.layout.Updatemenu(
                buttons=buttonlist_x_scale,
                #type = "buttons",
                direction="down",
                #pad={"r": 10, "t": 10},
                showactive=True,
                x=0.9,
                xanchor="right",
                y=1.12,
                yanchor="top"
            ),
            go.layout.Updatemenu(
                buttons=buttonlist_y_scale,
                direction="down",
                #type = "buttons",
                #pad={"r": 10, "t": 10},
                showactive=True,
                x=1,
                xanchor="right",
                y=1.12,
                yanchor="top"
            ),

        ])
autosize=True




fig.update_layout(
    annotations=[
        dict(text="<b>x-axis</b>", x=0, xref="paper", y=1.19, yref="paper",
             showarrow=False, align="left"),
        dict(text="<b>y-axis</b>", x=0.25, xref="paper", y=1.19, yref="paper",
                             showarrow=False, align = "left"),
        dict(text="size", x=0.54, xref="paper", y=1.06, yref="paper",
                             showarrow=False),
        dict(text="<b>x-axis scale</b>", x=0.9, xref="paper", y=1.19, yref="paper",
             showarrow=False, align="center"),
        dict(text="<b>y-axis scale</b>", x=1, xref="paper", y=1.19, yref="paper",
                             showarrow=False)    ])
  #fig.show()



fig.update_traces(
    hovertemplate =
                "<b>%{customdata[0]}</b><br>" +
                "Manual Assignment: <b>%{customdata[1]}</b><br>"+
                "TargetP prediction: %{customdata[2]}<br>"+
                "HHPred annotation: %{customdata[3]}<br>"+
                "TMHMM prediction: %{customdata[4]}<br>"+
                "TM likelihood (% in): %{customdata[5]}<br>"+
                "VLC protein abundnace: %{customdata[6]:.3E}<br>"+
                "Number of RBMs: %{customdata[7]}<br>"+
                "Chlamy top hit description:<br> %{customdata[8]}<br>",
    mode='markers'
)

#fig.update_yaxes(ticklabelposition="inside top", title=None)

#fig.update_xaxes(title_text='Log<sub>2</sub>FC (VLCO<sub>2</sub> vs. HCO<sub>2</sub>)')
#fig.update_yaxes(title_text='-log<sub>10</sub>[Adjusted <i>P</i> Value]')
#fig.update_yaxes(range=[0, 10])


app.layout = html.Div([
    dcc.Graph(
        id='life-exp-vs-gdp',
        figure=fig
    )
])

if __name__ == '__main__':
    app.run(debug=True)
