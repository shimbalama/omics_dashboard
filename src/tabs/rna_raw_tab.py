from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from src.read_files import RNASeqData
from src.helpers import (
    make_list_of_dicts,
    draw_box_chart,
    gene_dropdown,
    gene_dropdown_default
)
from ..components import ids

KEY = "rna_bulk"

def render(app: Dash, data: dict[str, dict[str, RNASeqData]]) -> html.Div:
    # see https://dash.plotly.com/basic-callbacks#dash-app-with-chained-callbacks


    gene_dropdown(app, ids.GENE_DROPDOWN, ids.RAW_RNA_DATA_DROP, data[KEY])

    @app.callback(
        Output(ids.COMPARISON_DROPDOWN, "options"),
        Input(ids.RAW_RNA_DATA_DROP, "value"),
    )
    def set_comparison_options(experiment: str) -> list[dict[str, str]]:
        """Populates the comparison selection dropdown with options from teh given dataset"""
        return make_list_of_dicts(list(data[KEY][experiment].comparisons))

    gene_dropdown_default(app, ids.GENE_DROPDOWN)

    @app.callback(
        Output(ids.COMPARISON_DROPDOWN, "value"),
        Input(ids.COMPARISON_DROPDOWN, "options"),
        Input(ids.SELECT_ALL_COMPARISONS_BUTTON, "n_clicks"),
    )
    def select_comparison_values(
        available_comparisons: list[dict[str, str]], _: int
    ) -> list[dict[str, str]]:
        """Default to all available comparisons"""
        return [comp["value"] for comp in available_comparisons]

    @app.callback(
        Output(ids.BOX_CHART, "children"),
        Input(ids.RAW_RNA_DATA_DROP, "value"),
        Input(ids.GENE_DROPDOWN, "value"),
        Input(ids.COMPARISON_DROPDOWN, "value"),
    )
    def update_box_chart(dataset_choice: str, gene: str, comps: list[str]) -> html.Div:
        """Re draws a box and wisker of the CPM data for each set of replicates for eact
        comparison and overlays the respective FDR value"""
        selected_data = data[KEY][dataset_choice]
        filtered: RNASeqData = selected_data.filter(comps)

        return draw_box_chart(gene, filtered)

    default = list(data[KEY].keys())
    return html.Div(
        children=[
            html.H6("Dataset"),
            dcc.Dropdown(
                id=ids.RAW_RNA_DATA_DROP,
                options=default,
                value=default[0],
                multi=False,
            ),
            html.H6("Gene"),
            dcc.Dropdown(
                id=ids.GENE_DROPDOWN,
            ),
            html.H6("Comparison"),
            dcc.Dropdown(
                id=ids.COMPARISON_DROPDOWN,
                multi=True,
            ),
            html.Button(
                className="dropdown-button",
                children=["Select All"],
                id=ids.SELECT_ALL_COMPARISONS_BUTTON,
                n_clicks=0,
            ),
            html.Div(
                draw_box_chart(
                    data[KEY][default[0]].df.columns[0],
                    data[KEY][default[0]],
                )
            ),
        ],
    )
