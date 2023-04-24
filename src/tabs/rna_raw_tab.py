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
    IDs,
    Params,
    dropdowns,
)
from ..components import ids

KEY = "rna_bulk"

PARAMs = Params(X="test")


def render(app: Dash, data: dict[str, dict[str, Data]], ids: IDs) -> html.Div:
    gene_dropdown(app, ids, data[KEY])
    gene_dropdown_default(app, ids)
    test_dropdown(app, ids, data[KEY])
    test_dropdown_select_all(app, ids)
    box(app, ids, data[KEY], PARAMs)

    # first_gene, dataset, datasets = get_defaults(data, KEY)
    return html.Div(children=dropdowns(data, KEY, PARAMs, ids))
