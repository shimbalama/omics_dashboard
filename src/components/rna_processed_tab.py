from dash import Dash, dcc, html
from dash.dependencies import Input, Output
import pandas as pd
from . import ids
import plotly.express as px
import dash_bio


def render(app: Dash, dfs: dict[str, pd.DataFrame]) -> html.Div:
    @app.callback(
        Output("vp-graph", "figure"),
        Input("plots", "value"),
        Input("effect-size", "value"),
        Input("P-val", "value"),
        Input(ids.PROCESSED_RNA_DATA_DROP, "value"),
    )
    def update_graph(plot, effect_lims, genomic_line, datadset_id):
        """Update rendering of data points upon changing x-value of vertical dashed lines."""
        df = dfs[datadset_id].copy()
        print(df.head(2))
        if plot == "Vol":
            return dash_bio.VolcanoPlot(
                dataframe=df,
                genomewideline_value=float(genomic_line),
                effect_size_line=list(map(float, effect_lims)),
                snp=None,
                gene="gene_symbol",
            )
        elif plot == "MA":
            # Define significance threshold
            alpha = 0.05
            df['color'] = ['red' if pval <= alpha else 'blue' for pval in df['P']]

            return px.scatter(df, y="EFFECTSIZE", x="P", color='color')
        else:
            raise ValueError("wtf")

    return html.Div(
        children=[
            html.H6("Dataset"),
            dcc.Dropdown(
                id=ids.PROCESSED_RNA_DATA_DROP,
                options=list(dfs.keys()),
                value="qlf.APAvsCS",
                multi=False,
            ),
            html.H6("Plot_type"),
            dcc.Dropdown(
                id="plots",
                options=["Vol", "MA"],
                value="Vol",
                multi=False,
            ),
            html.H6("P-val"),
            dcc.Slider(
                id="P-val",
                value=4,
                max=10,
                min=0,
                step=0.01,
                marks={str(num): str(num) for num in range(0, 11, 2)},
            ),
            html.H6("Effect-size"),
            dcc.RangeSlider(
                id="effect-size",
                min=-4,
                max=4,
                value=[-1, 1],
                step=0.01,
                marks={str(num): str(num) for num in range(-4, 5)},
            ),
            html.Div(
                dcc.Graph(
                    id="vp-graph",
                    figure=dash_bio.VolcanoPlot(
                        dataframe=dfs["qlf.APAvsCS"], snp=None, gene="gene_symbol"
                    ),
                )
            ),
        ],
    )
