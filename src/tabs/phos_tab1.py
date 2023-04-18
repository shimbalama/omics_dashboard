from dash import Dash, dcc, html
from ..components import ids
from src.read_files import ProtData, Data
from src.helpers import (
    draw_box_chart,
    gene_dropdown,
    gene_dropdown_default,
    box,
    Params,
)


KEY = "phosphoproteomics"

PARAMs = Params(DIV_ID="phosphoiddd222", X="gene", COLOUR="ID", LOG=True, Y="abun")


def render(app: Dash, data: dict[str, dict[str, Data]]) -> html.Div:
    # get comparisons
    gene_dropdown(
        app, "phosphoproteomics_gene_drop", "phosphoproteomics_dataset_drop", data[KEY]
    )
    gene_dropdown_default(app, "phosphoproteomics_gene_drop")
    box(
        app,
        "phosphoproteomics_dataset_drop",
        "phosphoproteomics_gene_drop",
        data[KEY],
        PARAMs,
    )
    default = list(data[KEY].keys())
    default_data = data[KEY][default[0]].filter("AAK1")
    return html.Div(
        children=[
            html.H6("Dataset"),
            dcc.Dropdown(
                id="phosphoproteomics_dataset_drop",
                options=default,
                value=default[0],
                multi=False,
            ),
            html.P("You can use prot to interactively " ""),
            html.H6("Gene"),
            dcc.Dropdown(
                id="phosphoproteomics_gene_drop",
            ),
            html.Div(draw_box_chart(default_data, "abun", PARAMs)),
        ],
    )


# , default_data,"abun", "phosphoiddd222", "gene", "ID", True
