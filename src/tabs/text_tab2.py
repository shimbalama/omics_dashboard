from dash import Dash, html
from pathlib import Path
from src.parse_data.read_files import Data
from src.helpers import Params


def render(app: Dash, data: dict[str, Data], params: Params) -> html.Div:
    # rules = Path(".")
    # print(list(rules.glob('*')))
    # return html.Div(children=[html.P(line) for line in rules.split("\n")])
    return html.Div('ff')