from dash import Dash, dcc, html
from ..components import ids
from src.read_files import ProtData, Data
from src.helpers import (
    draw_box_chart,
    gene_dropdown,
    gene_dropdown_default,
    box,
    get_defaults,
    test_dropdown,
    test_dropdown_select_all,
    Params,
)

KEY = "proteomics"

PARAMs = Params(DIV_ID="iddd222", X="test")


def render(app: Dash, data: dict[str, dict[str, Data]]) -> html.Div:
    gene_dropdown(app, "proteomics_gene_drop", "proteomics_dataset_drop", data[KEY])
    gene_dropdown_default(app, "proteomics_gene_drop")
    test_dropdown(app, "proteomics_test_drop", "proteomics_dataset_drop", data[KEY])
    test_dropdown_select_all(
        app,
        "proteomics_test_drop",
        "proteomics_test_drop",
        "proteomics_select_all",
    )
    box(
        app,
        "proteomics_dataset_drop",
        "proteomics_gene_drop",
        "proteomics_test_drop",
        data[KEY],
        PARAMs,
    )
    first_gene, dataset, datasets = get_defaults(data, KEY)
    return html.Div(
        children=[#this can be func
            html.H6("Dataset"),
            dcc.Dropdown(
                id="proteomics_dataset_drop",
                options=datasets,
                value=datasets[0],
                multi=False,
            ),
            html.P("You can use prot to interactively"),
            html.H6("Gene"),
            dcc.Dropdown(
                id="proteomics_gene_drop",
            ),
            html.H6("test"),
            dcc.Dropdown(
                id="proteomics_test_drop",
                multi=True,
            ),
            html.Button(
                className="dropdown-button",
                children=["Select All"],
                id="proteomics_select_all",
                n_clicks=0,
            ),
            html.Div(draw_box_chart(dataset, first_gene, PARAMs)),
        ],
    )


