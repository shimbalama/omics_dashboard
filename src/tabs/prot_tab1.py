from dash import Dash, dcc, html
from ..components import ids
from src.read_files import ProtData, Data
from src.helpers import draw_box_chart, gene_dropdown, gene_dropdown_default, box, Params

KEY = "proteomics"

PARAMs = Params(DIV_ID="iddd222", X="test")

def render(app: Dash, data: dict[str, dict[str, Data]]) -> html.Div:
    gene_dropdown(app, "proteomics_gene_drop", "proteomics_dataset_drop", data[KEY])
    gene_dropdown_default(app, "proteomics_gene_drop")
    box(app, "proteomics_dataset_drop", "proteomics_gene_drop", data[KEY], PARAMs)
    default = list(data[KEY].keys())
    return html.Div(
        children=[
            html.H6("Dataset"),
            dcc.Dropdown(
                id="proteomics_dataset_drop",
                options=default,
                value=default[0],
                multi=False,
            ),
            html.P("You can use prot to interactively"),
            html.H6("Gene"),
            dcc.Dropdown(
                id="proteomics_gene_drop",
            ),
            html.Div(draw_box_chart(data[KEY][default[0]], "TTN", PARAMs)),
        ],
    )


# draw_box_chart("TTN", data[KEY][default[0]], "iddd222", "test")
