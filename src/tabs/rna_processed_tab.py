from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from ..components import ids
import dash_bio
from src.read_files import RNASeqData
import pandas as pd
from src.helpers import make_list_of_dicts

KEY = 'rna_bulk'


def render(app: Dash, data: dict[str, RNASeqData], ids2) -> html.Div:
    # get comparisons
    def draw_volcano(df, genomic_line, effect_lims):
        return dash_bio.VolcanoPlot(
            dataframe=df,
            genomewideline_value=float(genomic_line),
            effect_size_line=list(map(float, effect_lims)),
            snp=None,
            gene="gene_symbol",
            width=1111,
            height=888,
        )

    @app.callback(
        Output(ids.PROCESSED_COMPARISON_DROPDOWN, "options"),
        Input(ids.PROCESSED_RNA_DATA_DROP, "value"),
    )
    def set_comparison_options(experiment: str) -> list[dict[str, str]]:
        """Set the test option by given experiemnt name"""
        return make_list_of_dicts(list(data[KEY][experiment].degs))

    @app.callback(
        Output(ids.PROCESSED_COMPARISON_DROPDOWN, "value"),
        Input(ids.PROCESSED_COMPARISON_DROPDOWN, "options"),
    )
    def select_gene_value(options: list[dict[str, str]]) -> str:
        """Select first test as default value after changing dataset"""
        return options[0]["value"]

    @app.callback(
        Output(ids.VOLCANO, "figure"),
        Input(ids.PROCESSED_COMPARISON_DROPDOWN, "value"),
        Input(ids.EEFECT_SIZE, "value"),
        Input(ids.PVALUE, "value"),
        Input(ids.PROCESSED_RNA_DATA_DROP, "value"),
    )
    def update_graph(
        comp: str,
        effect_lims: list[str],
        genomic_line: str,
        datadset_id: str,
    ):
        """Update rendering of data points upon changing x-value of vertical dashed lines."""
        selcted_data: RNASeqData = data[KEY][datadset_id]
        df: pd.DataFrame = selcted_data.processed_dfs[comp].copy()
        df['gene_symbol'] = df.index
        df = df.reset_index(drop=True)
        return draw_volcano(df, genomic_line, effect_lims)

    experiments = list(data[KEY].keys())
    first_experiment = experiments[0]
    default_comparison = list(data[KEY][first_experiment].processed_dfs.keys())[0]
    default_df = data[KEY][first_experiment].processed_dfs[default_comparison]
    default_df['gene_symbol'] = default_df.index
    default_df = default_df.reset_index(drop=True)
    #print(111122223333,experiments,first_experiment,default_comparison,default_df)
    return html.Div(
        children=[
            html.P(
                "You can use Volcano Plot to interactively "
                "identify clinically meaningful markers in "
                "genomic experiments, i.e., markers that are "
                "statistically significant and have an effect "
                "size greater than some threshold. "
                "Specifically, volcano plots depict the negative "
                "log-base-10 p-values plotted against their "
                "effect size."
            ),
            html.H6("Dataset"),
            dcc.Dropdown(
                id=ids.PROCESSED_RNA_DATA_DROP,
                options=experiments,
                value=first_experiment,
                multi=False,
            ),
            html.H6("test"),
            dcc.Dropdown(
                id=ids.PROCESSED_COMPARISON_DROPDOWN,
                multi=False,
                value=default_comparison,
            ),
            html.H6("P-val"),
            dcc.Slider(
                id=ids.PVALUE,
                value=4,
                max=10,
                min=0,
                step=0.01,
                marks={str(num): str(num) for num in range(0, 11, 2)},
            ),
            html.H6("Effect-size"),
            dcc.RangeSlider(
                id=ids.EEFECT_SIZE,
                min=-4,
                max=4,
                value=[-1, 1],
                step=0.01,
                marks={str(num): str(num) for num in range(-4, 5)},
            ),
            html.Div(
                dcc.Graph(
                    id=ids.VOLCANO,
                    figure=draw_volcano(default_df, 1, [-3, 3]),
                )
            ),
        ],
    )
