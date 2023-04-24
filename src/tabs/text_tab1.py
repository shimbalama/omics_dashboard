from dash import Dash, html

from src.read_files import RNASeqData


def render(app: Dash, data: dict[str, RNASeqData], ids) -> html.Div:
    return html.Div(
        children=[
            html.P("Wecome to the lab Hudson ploting dashbouard! "),
        ]
    )
