from dash import Dash, dcc, html
from ..components import ids
from src.read_files import ProtData, Data
from src.helpers import (
    draw_box_chart,
)

KEY = "phosphoproteomics"


def render(app: Dash, data: dict[str, dict[str, Data]]) -> html.Div:
    # get comparisons

    default = list(data[KEY].keys())
    default_data = data[KEY][default[0]].filter("AAK1")
    print(1111111111, default_data.df.head(3))
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
            html.Div(
                draw_box_chart(
                    "abun", default_data, "phosphoiddd222", "gene", "ID", True
                )
            ),
        ],
    )
