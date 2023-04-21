from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from src.read_files import Data
from src.helpers import (
    test_dropdown,
    test_dropdown_select_all,
    draw_box_chart,
    gene_dropdown,
    gene_dropdown_default,
    box,
    Params,
)
from ..components import ids

KEY = "rna_bulk"

PARAMs = Params(DIV_ID=ids.BOX_CHART, X="test")


def render(app: Dash, data: dict[str, dict[str, Data]]) -> html.Div:
    # see https://dash.plotly.com/basic-callbacks#dash-app-with-chained-callbacks

    gene_dropdown(app, ids.GENE_DROPDOWN, ids.RAW_RNA_DATA_DROP, data[KEY])
    gene_dropdown_default(app, ids.GENE_DROPDOWN)

    test_dropdown(app, ids.COMPARISON_DROPDOWN, ids.RAW_RNA_DATA_DROP, data[KEY])
    test_dropdown_select_all(
        app,
        ids.COMPARISON_DROPDOWN,
        ids.COMPARISON_DROPDOWN,
        ids.SELECT_ALL_COMPARISONS_BUTTON,
    )
  
    box()
    # @app.callback(
    #     Output(ids.BOX_CHART, "children"),
    #     Input(ids.RAW_RNA_DATA_DROP, "value"),
    #     Input(ids.GENE_DROPDOWN, "value"),
    #     Input(ids.COMPARISON_DROPDOWN, "value"),
    # )
    # def update_box_chart(dataset_choice: str, gene: str, comps: list[str]) -> html.Div:
    #     """Re draws a box and wisker of the CPM data for each set of replicates for eact
    #     test and overlays the respective FDR value"""
    #     selected_data = data[KEY][dataset_choice]
    #     filtered: Data = selected_data.filter(comps)

    #     return draw_box_chart(filtered, gene, PARAMs)

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
            html.H6("test"),
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
                    data[KEY][default[0]].filter(
                        list(data[KEY][default[0]].comparisons)
                    ),
                    data[KEY][default[0]].df.columns[0],
                    PARAMs,
                )
            ),
        ],
    )
