import plotly.express as px
from dash import Dash, dcc, html
from dash.dependencies import Input, Output

from . import ids

MEDAL_DATA = px.data.medals_long()


def render(app: Dash) -> html.Div:
    @app.callback(
        Output(ids.BOX_CHART, "children"),
        [
            Input('tmp', "value"),
        ],
    )
    def update_box_chart(nations: list[str]) -> html.Div:
        filtered_data = MEDAL_DATA.query("nation in @nations")

        if filtered_data.shape[0] == 0:
            return html.Div("No data selected.", id=ids.BOX_CHART)

        fig = px.box(filtered_data, x="medal", y="count", color="nation", text="nation")

        return html.Div(dcc.Graph(figure=fig), id=ids.BOX_CHART)

    return html.Div(id=ids.BOX_CHART)
