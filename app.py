from dash import Dash, dcc, html, Input, Output, State
import plotly.express as px
import pandas as pd
import textwrap
import numpy as np
import dash_bootstrap_components as dbc

app = Dash(external_stylesheets=[dbc.themes.COSMO])

df = pd.read_csv('Chlorella_db_WIP_manual.csv')

# Filtering ribosomal content
strings_to_filter = ['ribosomal', 'ribosome', 'RNA']
columns_to_check = ['HHPred_description', 'kogdefline', 'iprDesc']
df = df[~df[columns_to_check].apply(lambda x: x.str.contains('|'.join(strings_to_filter), case=False, na=False)).any(axis=1)]

# Define plot variables
x_plot = 'PYR_vs_LYS__fraction_log2 fold change'
y_plot = 'Rb_vs_Bst_log2 fold change'
size_plot = 'E439_VLC_avg'
color_plot = 'Log2FC(10030/5per)'

import math

# Define rounded bounds for cleaner sliders
x_min = math.floor(df[x_plot].min())
x_max = math.ceil(df[x_plot].max())

y_min = math.floor(df[y_plot].min())
y_max = math.ceil(df[y_plot].max())

# Define ticks spaced every N units
def generate_marks(start, end, step):
    return {i: str(i) for i in range(start, end + 1, step)}



df['Cr_BLAST_hit_description'] = df['Cr_BLAST_hit_description'].apply(lambda t: "<br>&nbsp;".join(textwrap.wrap(str(t), 80)))

# List of numeric columns suitable for x,y,size,color plotting
numeric_cols = df.select_dtypes(include=['float64', 'int64']).columns.tolist()

# App layout with filters
app.layout = html.Div([

    # Row 1: Axis selectors
    html.Div([
        html.Div([
            html.Label("X-axis:"),
            dcc.Dropdown(
                id='x-axis',
                options=[{'label': col, 'value': col} for col in numeric_cols],
                value='PYR_vs_LYS__fraction_log2 fold change',
                clearable=False
            ),
        ], style={'width': '24%', 'padding': '10px'}),

        html.Div([
            html.Label("Y-axis:"),
            dcc.Dropdown(
                id='y-axis',
                options=[{'label': col, 'value': col} for col in numeric_cols],
                value='Rb_vs_Bst_log2 fold change',
                clearable=False
            ),
        ], style={'width': '24%', 'padding': '10px'}),

        html.Div([
            html.Label("Color:"),
            dcc.Dropdown(
                id='color',
                options=[{'label': col, 'value': col} for col in numeric_cols],
                value='Log2FC(10030/5per)',
                clearable=False
            ),
        ], style={'width': '24%', 'padding': '10px'}),

        html.Div([
            html.Label("Size:"),
            dcc.Dropdown(
                id='size',
                options=[{'label': col, 'value': col} for col in numeric_cols],
                value='E439_VLC_avg',
                clearable=False
            ),
        ], style={'width': '24%', 'padding': '10px'}),
    ], style={'display': 'flex', 'justifyContent': 'space-between'}),

    # Row 2: Search and TargetP2 filters side by side
    html.Div([
        html.Div([
            html.Label("Search Protein ID"),
            dcc.Input(
                id='protein-search',
                type='text',
                placeholder='Enter Protein ID...',
                debounce=True,
                style={'width': '100%'}
            )
        ], style={'width': '48%', 'padding': '10px'}),

        html.Div([
            html.Label("TargetP2 Prediction"),
            dcc.Dropdown(
                id='targetp2-dropdown',
                options=[{'label': val, 'value': val} for val in sorted(df['TargetP2_prediction'].dropna().unique())],
                multi=True,
                placeholder='Filter by TargetP2 prediction...',
                style={'width': '100%'}
            )
        ], style={'width': '48%', 'padding': '10px'}),
    ], style={'display': 'flex', 'justifyContent': 'space-between'}),
    html.Div([
        html.Label("Search Metadata (description, annotation, etc.)"),
        dcc.Input(
            id='metadata-search',
            type='text',
            placeholder='Enter keyword(s) for metadata...',
            debounce=True,
            style={'width': '100%'}
        )
    ], style={'width': '100%', 'padding': '10px'}),

    # Row 3: Sliders for filtering X and Y values with manual inputs
    html.Div([
        html.Div([
            html.Label(id='x-slider-label'),
            html.Div([
                dcc.Input(id='x-min-input', type='number', debounce=True, placeholder='Min', style={'width': '80px'}),
                dcc.Input(id='x-max-input', type='number', debounce=True, placeholder='Max', style={'width': '80px'}),
            ], style={'display': 'flex', 'gap': '10px', 'marginBottom': '5px'}),
            dcc.RangeSlider(
                id='x-range-slider',
                step=0.1,  # allows float values
                tooltip={"placement": "bottom", "always_visible": False},
                allowCross=False
            )
        ], style={'width': '48%', 'padding': '10px'}),

        html.Div([
            html.Label(id='y-slider-label'),
            html.Div([
                dcc.Input(id='y-min-input', type='number', debounce=True, placeholder='Min', style={'width': '80px'}),
                dcc.Input(id='y-max-input', type='number', debounce=True, placeholder='Max', style={'width': '80px'}),
            ], style={'display': 'flex', 'gap': '10px', 'marginBottom': '5px'}),
            dcc.RangeSlider(
                id='y-range-slider',
                step=0.1,
                tooltip={"placement": "bottom", "always_visible": False},
                allowCross=False
            )
        ], style={'width': '48%', 'padding': '10px'}),
    ], style={'display': 'flex', 'justifyContent': 'space-between'}),



    # Graph and point count

    html.Div(id='point-count', style={'padding': '10px', 'fontWeight': 'bold'}),
    html.Div([
        html.Button("Download FASTA", id="download-button", n_clicks=0, className="btn btn-primary"),
        dcc.Download(id="fasta-download")
    ], style={'padding': '10px'}),
    dcc.Graph(id='scatter-plot'),
    dcc.Store(id='filtered-data-store'),
    html.Div(id='sequence-display', style={'whiteSpace': 'pre-wrap', 'padding': '10px', 'border': '1px solid #ccc', 'marginTop': '20px', 'maxHeight': '300px', 'overflowY': 'auto'})


])
@app.callback(
    Output('scatter-plot', 'figure'),
    Output('point-count', 'children'),
    Output('filtered-data-store', 'data'),
    Input('x-axis', 'value'),
    Input('y-axis', 'value'),
    Input('color', 'value'),
    Input('size', 'value'),
    Input('protein-search', 'value'),       # existing ID search
    Input('metadata-search', 'value'),      # <-- new input
    Input('x-range-slider', 'value'),
    Input('y-range-slider', 'value'),
    Input('targetp2-dropdown', 'value')
)
def update_figure(x_col, y_col, color_col, size_col,
                  protein_id_search, metadata_search,
                  x_range, y_range, targetp2_filter):

    # Fallback for ranges if None
    if x_range is None:
        x_range = [df[x_col].min(), df[x_col].max()]
    if y_range is None:
        y_range = [df[y_col].min(), df[y_col].max()]

    x_min, x_max = x_range
    y_min, y_max = y_range

    # Filter DataFrame by slider ranges
    filtered_df = df[
        (df[x_col] >= x_min) & (df[x_col] <= x_max) &
        (df[y_col] >= y_min) & (df[y_col] <= y_max)
    ]

    # Filter by TargetP2 prediction if filter is set
    if targetp2_filter and len(targetp2_filter) > 0:
        filtered_df = filtered_df[filtered_df['TargetP2_prediction'].isin(targetp2_filter)]

    search_columns = [
        'identifier', 'HHPred_description', 'Cr_BLAST_hit',
        'Cr_BLAST_hit_description', 'kogdefline', 'iprDesc',
        'definition', 'pathway', 'prediction', 'TargetP2_prediction',
        'manual_assignment', 'manual_notes'
    ]

    # === Filter by metadata search ===
    search_columns = [
        'identifier', 'HHPred_description', 'Cr_BLAST_hit',
        'Cr_BLAST_hit_description', 'kogdefline', 'iprDesc',
        'definition', 'pathway', 'prediction', 'TargetP2_prediction',
        'manual_assignment', 'manual_notes'
    ]

    if metadata_search and isinstance(metadata_search, str):
        metadata_mask = filtered_df[search_columns].apply(
            lambda row: row.astype(str).str.contains(metadata_search, case=False, na=False)
            , axis=1).any(axis=1)
        filtered_df = filtered_df[metadata_mask]

    if protein_id_search and isinstance(protein_id_search, str):
        highlight_df = filtered_df[filtered_df['identifier'].str.contains(protein_id_search, case=False, na=False)]
        others_df = filtered_df[~filtered_df['identifier'].str.contains(protein_id_search, case=False, na=False)]
    else:
        highlight_df = pd.DataFrame(columns=filtered_df.columns)
        others_df = filtered_df

    # === Add dynamic min/max for color and size ===
    # Apply log scale safely (avoid log(0) or negative)
    def safe_log_scale(sizes):
        return np.log10(sizes.clip(min=1))  # clip at 1 to avoid log(0)

    if not filtered_df.empty:
        color_min = filtered_df[color_col].min()
        color_max = filtered_df[color_col].max()
        size_min = filtered_df[size_col].min()
        size_max = filtered_df[size_col].max()
    else:
        # fallback values if empty filter
        color_min, color_max = 0, 1
        size_min, size_max = 0, 1

    if not filtered_df.empty and size_max > size_min:
        # Normalize size between 0 and 1 (or scale it)
        filtered_df['size_norm'] = (filtered_df[size_col] - size_min) / (size_max - size_min)
        size_col_to_use = 'size_norm'
    else:
        filtered_df['size_norm'] = 0.5  # constant size fallback
        size_col_to_use = 'size_norm'

    # Create scatter plot with Plotly Express
    fig = px.scatter(
        others_df,
        x=x_col,
        y=y_col,
        size=size_col,
        color=color_col,
        color_continuous_scale='plasma',
        size_max=150,
        custom_data=[
            'identifier', 'HHPred_description', 'Cr_BLAST_hit',
            "Cr_BLAST_hit_description", 'transmembrane_likelihood',
            'E439_VLC_avg', 'number_RBMs', 'kogdefline', 'iprDesc',
            'definition', 'pathway', 'prediction', 'TargetP2_prediction',
            'manual_assignment', x_col, y_col, size_col, color_col, 'manual_notes'
        ],
        height=800
    )

    # Highlight matching proteins if search matches exist
    if not highlight_df.empty:
        fig.add_scatter(
            x=highlight_df[x_col],
            y=highlight_df[y_col],
            mode='markers+text',
            text=highlight_df['identifier'],
            textposition='top center',
            marker=dict(size=20, color='red', line=dict(color='black', width=2)),
            name='Search Match',
            customdata=highlight_df[[
                'identifier', 'HHPred_description', 'Cr_BLAST_hit',
                "Cr_BLAST_hit_description", 'transmembrane_likelihood',
                'E439_VLC_avg', 'number_RBMs', 'kogdefline', 'iprDesc',
                'definition', 'pathway', 'prediction', 'TargetP2_prediction',
                'manual_assignment', x_col, y_col, size_col, color_col, 'manual_notes'
            ]].values
        )

    fig.update_layout(
        hovermode='closest',
        legend=dict(title=None),
        hoverlabel=dict(bgcolor='rgba(255,255,255,0.75)', font=dict(color='black'))
    )

    fig.update_traces(
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            "Manual assignment: <b>%{customdata[13]}</b><br>"
            "Manual notes: <b>%{customdata[18]}</b><br>"
            "HHPred annotation: %{customdata[1]}<br>"
            "Chlamy top hit: %{customdata[2]}<br>"
            "Chlamy top hit description:<br> %{customdata[3]}<br>"
            "TM likelihood (% in): %{customdata[4]}<br>"
            "VLC protein abundance: %{customdata[5]:.3E}<br>"
            "Number of RBMs: %{customdata[6]}<br>"
            "KOG annotation: %{customdata[7]}<br>"
            "InterPro annotation: %{customdata[8]}<br>"
            "EC annotation: %{customdata[9]}<br>"
            "GO annotation: %{customdata[10]}<br>"
            "SignalP prediction: %{customdata[11]}<br>"
            "TargetP2 prediction: %{customdata[12]}<br>"
            f"{x_col}: %{{customdata[14]}}<br>"
            f"{y_col}: %{{customdata[15]}}<br>"
            f"{size_col}: %{{customdata[16]}}<br>"
            f"{color_col}: %{{customdata[17]}}<br>"
        )
    )

    num_points = len(filtered_df)
    return fig, f"Number of proteins shown: {num_points}", filtered_df.to_json(date_format='iso', orient='split')


@app.callback(
    Output('sequence-display', 'children'),
    Input('scatter-plot', 'clickData')
)
def display_sequence(clickData):
    if clickData is None:
        return "Click a point to see protein sequence."

    protein_id = clickData['points'][0]['customdata'][0]  # Assuming identifier is first in customdata

    sequence = df.loc[df['identifier'] == protein_id, 'sequence'].values
    if len(sequence) == 0:
        return "Sequence not found."

    seq_text = sequence[0]

    # Format as fasta style
    fasta_text = f">{protein_id}\n{seq_text}"
    return fasta_text

@app.callback(
    Output('x-range-slider', 'min'),
    Output('x-range-slider', 'max'),
    Output('x-range-slider', 'value'),
    Output('x-range-slider', 'marks'),
    Output('x-slider-label', 'children'),
    Output('x-min-input', 'value'),
    Output('x-max-input', 'value'),

    Input('x-axis', 'value'),
    Input('x-min-input', 'value'),
    Input('x-max-input', 'value'),
    State('x-range-slider', 'value')
)
def update_x_slider(x_col, min_val, max_val, current):
    x_min_data = math.floor(df[x_col].min())
    x_max_data = math.ceil(df[x_col].max())
    x_marks = generate_marks(x_min_data, x_max_data, step=2)

    # Determine value to apply to slider
    if min_val is not None and max_val is not None and min_val < max_val:
        x_slider_value = [min_val, max_val]
    else:
        x_slider_value = [x_min_data, x_max_data]

    return (
        x_min_data, x_max_data, x_slider_value,
        x_marks, f"X filter: {x_col}", x_slider_value[0], x_slider_value[1]
    )


@app.callback(
    Output('y-range-slider', 'min'),
    Output('y-range-slider', 'max'),
    Output('y-range-slider', 'value'),
    Output('y-range-slider', 'marks'),
    Output('y-slider-label', 'children'),
    Output('y-min-input', 'value'),
    Output('y-max-input', 'value'),

    Input('y-axis', 'value'),
    Input('y-min-input', 'value'),
    Input('y-max-input', 'value'),
    State('y-range-slider', 'value')
)
def update_y_slider(y_col, min_val, max_val, current):
    y_min_data = math.floor(df[y_col].min())
    y_max_data = math.ceil(df[y_col].max())
    y_marks = generate_marks(y_min_data, y_max_data, step=2)

    # Determine value to apply to slider
    if min_val is not None and max_val is not None and min_val < max_val:
        y_slider_value = [min_val, max_val]
    else:
        y_slider_value = [y_min_data, y_max_data]

    return (
        y_min_data, y_max_data, y_slider_value,
        y_marks, f"Y filter: {y_col}", y_slider_value[0], y_slider_value[1]
    )


@app.callback(
    Output("fasta-download", "data"),
    Input("download-button", "n_clicks"),
    State("filtered-data-store", "data"),
    prevent_initial_call=True
)
def download_fasta(n_clicks, stored_data):
    if stored_data is None:
        return dash.no_update

    filtered_df = pd.read_json(stored_data, orient='split')

    # Ensure required columns exist
    if 'identifier' not in filtered_df or 'sequence' not in filtered_df:
        return dash.no_update

    # Generate FASTA content
    fasta_lines = [
        f">{row['identifier']}\n{row['sequence']}"
        for _, row in filtered_df.iterrows()
        if pd.notnull(row['sequence'])
    ]
    fasta_text = "\n".join(fasta_lines)

    return dict(
        content=fasta_text,
        filename="filtered_proteins.fasta",
        type="text/plain"
    )

if __name__ == '__main__':
    app.run(debug=True)