from dash import Dash, html

from src.read_files import RNASeqData


def render(app: Dash, data: dict[str, RNASeqData], ids) -> html.Div:
    return html.Div(
        children=[
            html.P(
                "File names and internal column name are used by this code - please follow these conventions "
            ),
        ]
    )
