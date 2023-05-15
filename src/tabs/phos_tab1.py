from dash import Dash, dcc, html
from src.read_files import Data
from src.helpers import (
    draw_box_chart,
    gene_dropdown,
    gene_dropdown_default,
    box,
    get_defaults,
    Params,
    IDs,
)


KEY = "phosphoproteomics"

PARAMs = Params(X="gene", COLOUR="ID", LOG=True, Y="abun")


def render(app: Dash, data: dict[str, dict[str, Data]], ids: IDs) -> html.Div:
    # get comparisons
    gene_dropdown(app, ids, data[KEY])
    gene_dropdown_default(app, ids)
    box(app, ids, data[KEY], PARAMs)

    _, dataset, datasets = get_defaults(data, KEY)
    return html.Div(
        children=[
            html.H6("Dataset"),
            dcc.Dropdown(
                id=ids.data_drop,
                options=datasets,
                value=datasets[0],
                multi=False,
            ),
            html.P("You can use prot to interactively "),
            html.H6("Gene"),
            dcc.Dropdown(
                id=ids.gene_drop,
            ),
            html.H6("test"),
            dcc.Dropdown(
                id=ids.tests_drop,
                multi=True,
            ),
            html.Button(
                className="dropdown-button",
                children=["Select All"],
                id=ids.select_all,
                n_clicks=0,
            ),
            html.Div(draw_box_chart(dataset, "abun", PARAMs, ids.plot)),
        ],
    )


# , default_data,"abun", "phosphoiddd222", "gene", "ID", True
