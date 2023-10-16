from dash import Dash, html, get_asset_url
from random import randint
from src.parse.read_files import Data
from src.helpers import Params


def render(app: Dash, data: dict[str, Data], params: Params) -> html.Div:
    return html.Div(
        children=[
            html.P("Wecome to the Hudson lab ploting dashboard!"),
            html.P("This dashboard is designed to help you visualize \
                   the data from the Hudson lab experiments. \n"),
            html.Img(src=get_asset_url(f"h{randint(1, 10)}.png")),
        ]
    )
