from dash import Dash, html
from pathlib import Path
from src.read_files import RNASeqData


def render(app: Dash, data: dict[str, RNASeqData], ids, params) -> html.Div:
    # rules = Path(".")
    # print(list(rules.glob('*')))
    # return html.Div(children=[html.P(line) for line in rules.split("\n")])
    return html.Div('ff')