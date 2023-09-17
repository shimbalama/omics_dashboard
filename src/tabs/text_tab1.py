from dash import Dash, html, get_asset_url
from random import randint
from src.read_files import RNASeqData


def render(app: Dash, data: dict[str, RNASeqData], ids, params) -> html.Div:
    return html.Div(
        children=[
            html.P("Wecome to the Hudson lab ploting dashboard!"),
            html.P("This dashboard is designed to help you visualize \
                   the data from the Hudson lab experiments. \n"),
            html.Img(src=get_asset_url(f"h{randint(1, 9)}.png")),
        ]
    )
