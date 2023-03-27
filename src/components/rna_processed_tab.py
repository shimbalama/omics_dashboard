from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from . import ids
import plotly.express as px
import dash_bio
from src.read_files import RNASeqData
import pandas as pd


def render(app: Dash, data: dict[str, RNASeqData]) -> html.Div:
    # get comparisons
    @app.callback(
        Output("comp", "options"),
        Input(ids.PROCESSED_RNA_DATA_DROP, "value"),
    )
    def set_comparison_options(selected_comparison):
        selected_data = data[selected_comparison]

        return [{"label": comp, "value": comp} for comp in selected_data.degs]

    # select comparison/s
    @app.callback(
        Output("comp", "value"),
        Input("comp", "options"),
    )
    def select_comparison_values(available_options):
        return available_options

    @app.callback(
        Output("vp-graph", "figure"),
        Input("comp", "value"),
        Input("plots", "value"),
        Input("effect-size", "value"),
        Input("P-val", "value"),
        Input(ids.PROCESSED_RNA_DATA_DROP, "value"),
    )
    def update_graph(
        comp: str,
        plot: str,
        effect_lims: list[str],
        genomic_line: str,
        datadset_id: str,
    ):
        """Update rendering of data points upon changing x-value of vertical dashed lines."""

        selcted_data: RNASeqData = data[datadset_id]
        df: pd.DataFrame = selcted_data.processed_dfs[comp].copy()
        print(df.head())
        if plot == "Vol":
            return dash_bio.VolcanoPlot(
                dataframe=df,
                genomewideline_value=float(genomic_line),
                effect_size_line=list(map(float, effect_lims)),
                snp=None,
                gene="gene_symbol",
                width=1111,
                height=888,
            )
        elif plot == "MA":
            # Define significance threshold
            alpha = 0.05
            df["color"] = ["red" if pval <= alpha else "blue" for pval in df["P"]]

            return px.scatter(df, y="EFFECTSIZE", x="P", color="color")
        else:
            raise ValueError("wtf")

    return html.Div(
        children=[
            html.H6("Dataset"),
            dcc.Dropdown(
                id=ids.PROCESSED_RNA_DATA_DROP,
                options=list(data.keys()),
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
            html.H6("comp"),
            dcc.Dropdown(
                id="comp",
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
                        dataframe=data["BETi"].processed_dfs["qlf.APAvsCS"],
                        snp=None,
                        gene="gene_symbol",
                        width=1111,
                        height=888,
                    ),
                )
            ),
        ],
    )
