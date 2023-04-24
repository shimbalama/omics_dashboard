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
    get_defaults,
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

    box(
        app,
        ids.RAW_RNA_DATA_DROP,
        ids.GENE_DROPDOWN,
        ids.COMPARISON_DROPDOWN,
        data[KEY],
        PARAMs,
    )

    first_gene, dataset, datasets = get_defaults(data, KEY)
    return html.Div(
        children=[
            html.H6("Dataset"),
            dcc.Dropdown(
                id=ids.RAW_RNA_DATA_DROP,
                options=datasets,
                value=datasets[0],
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
                    dataset,
                    first_gene,
                    PARAMs,
                )
            ),
        ],
    )
