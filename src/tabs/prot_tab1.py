from dash import Dash, dcc, html

# from ..components import ids
from src.read_files import ProtData, Data
from src.helpers import (
    draw_box_chart,
    gene_dropdown,
    gene_dropdown_default,
    box,
    get_defaults,
    test_dropdown,
    test_dropdown_select_all,
    dropdowns,
    IDs,
    Params,
)

KEY = "proteomics"

PARAMs = Params(X="test")


def render(app: Dash, data: dict[str, dict[str, Data]], ids: IDs) -> html.Div:
    gene_dropdown(app, ids, data[KEY])
    gene_dropdown_default(app, ids)
    test_dropdown(app, ids, data[KEY])
    test_dropdown_select_all(app, ids)
    box(app, ids, data[KEY], PARAMs)
    return html.Div(children=dropdowns(data, KEY, PARAMs, ids))
    # first_gene, dataset, datasets = get_defaults(data, KEY)
    # return html.Div(
    #     children=[
    #         html.H6("Dataset"),
    #         dcc.Dropdown(
    #             id=ids.data_drop,
    #             options=datasets,
    #             value=datasets[0],
    #             multi=False,
    #         ),
    #         html.P("You can use prot to interactively"),
    #         html.H6("Gene"),
    #         dcc.Dropdown(
    #             id=ids.gene_drop,
    #         ),
    #         html.H6("test"),
    #         dcc.Dropdown(
    #             id=ids.tests_drop,
    #             multi=True,
    #         ),
    #         html.Button(
    #             className="dropdown-button",
    #             children=["Select All"],
    #             id=ids.select_all,
    #             n_clicks=0,
    #         ),
    #         html.Div(draw_box_chart(dataset, first_gene, PARAMs, ids.plot)),
    #     ],
    #)
