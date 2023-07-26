from dash import Dash, html, get_asset_url

from src.read_files import RNASeqData


def render(app: Dash, data: dict[str, RNASeqData], ids, params) -> html.Div:
    return html.Div(
        children=[
            html.P("Wecome to the lab Hudson ploting dashboard! \n"*4),
            html.Img(src=get_asset_url('h1.png'))
        ]
    )
