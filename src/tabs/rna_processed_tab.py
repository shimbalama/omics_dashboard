from dash import Dash, dcc, html
from dash.dependencies import Input, Output
import dash_bio
from src.parse.read_files import Data
import pandas as pd
from src.helpers import make_list_of_dicts
from icecream import ic
from src.helpers import Params

KEY = "rna_bulk"


def render(
    app: Dash, uninitialised_datasets: dict[str, Data], params: Params
) -> html.Div:
    """Render the layout for the RNA processed tab."""
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
        Output(params.ids.tests_drop, "options"),
        Input(params.ids.data_drop, "value"),
    )
    def set_comparison_options(experiment: str) -> list[dict[str, str]]:
        """Set the test option by given experiemnt name"""
        return make_list_of_dicts(list(datasets[experiment].degs))

    @app.callback(
        Output(params.ids.tests_drop, "value"),
        Input(params.ids.tests_drop, "options"),
    )
    def select_gene_value(options: list[dict[str, str]]) -> str:
        """Select first test as default value after changing dataset"""
        return options[0]["value"]

    @app.callback(
        Output(params.ids.plot, "figure"),
        Input(params.ids.tests_drop, "value"),
        Input(params.ids.slider1, "value"),
        Input(params.ids.slider2, "value"),
        Input(params.ids.data_drop, "value"),
    )
    def update_graph(
        comp: str,
        effect_lims: list[str],
        genomic_line: str,
        datadset_id: str,
    ):
        """Update rendering of data points upon changing x-value of vertical dashed lines."""
        selcted_data: Data = datasets[datadset_id]
        df: pd.DataFrame = selcted_data.processed_dfs[comp].copy()
        df["gene_symbol"] = df.index
        df = df.reset_index(drop=True)
        return draw_volcano(df, genomic_line, effect_lims)

    datasets = {}
    for k, v in uninitialised_datasets.items():
        func, path = v
        dataset = func(path)
        if not dataset.df_FDR.empty:
            datasets[k] = func(path)
    dataset_names = list(datasets.keys())
    first_dataset = datasets[dataset_names[0]]
    default_comparison = list(first_dataset.processed_dfs.keys())[0]
    default_df = first_dataset.processed_dfs[default_comparison]
    default_df["gene_symbol"] = default_df.index
    default_df = default_df.reset_index(drop=True)
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
                id=params.ids.data_drop,
                options=dataset_names,
                value=dataset_names[0],
                multi=False,
            ),
            html.H6("test"),
            dcc.Dropdown(
                id=params.ids.tests_drop,
                multi=False,
                value=default_comparison,
            ),
            html.H6("P-val"),
            dcc.Slider(
                id=params.ids.slider2,
                value=4,
                max=10,
                min=0,
                step=0.01,
                marks={str(num): str(num) for num in range(0, 11, 2)},
                tooltip={"placement": "bottom", "always_visible": False},
            ),
            html.H6("Effect-size"),
            dcc.RangeSlider(
                id=params.ids.slider1,
                min=-4,
                max=4,
                value=[-1, 1],
                step=0.01,
                marks={str(num): str(num) for num in range(-4, 5)},
                tooltip={"placement": "bottom", "always_visible": False},
            ),
            html.Div(
                dcc.Graph(
                    id=params.ids.plot,
                    figure=draw_volcano(default_df, 1, [-3, 3]),
                )
            ),
        ],
    )
