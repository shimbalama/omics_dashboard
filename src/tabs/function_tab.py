from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from src.read_files import Data
import pandas as pd
from src.helpers import make_list_of_dicts
import plotly.graph_objects as go
import plotly.express as px

KEY = "function"


def render(app: Dash, datasets: dict[str, Data], ids2, params) -> html.Div:
    def draw_line(data: Data) -> html.Div:
        """Draws a box and wisker of the CPM data for each set of replicates for eact
        test and overlays the respective FDR value"""
        datas = data.df.describe().to_dict()
        dose = data.dose.to_dict()
        dose["Dose"]["0"] = 0
        
        xs = [dose["Dose"][str(time_point)] for time_point in list(datas.keys())]
        ys = [v.get("50%") for v in datas.values()]
        stds = [v.get("std") for v in datas.values()]
        fig = go.Figure(
            data=go.Scatter(
                x=xs,
                y=ys,
                error_y=dict(
                    type="data",  # value of error bar given in data coordinates
                    array=stds,
                    visible=True,
                ),
            )
        )
        fig.update_xaxes(type="log")
        

        # fig = px.scatter(
        #     data.df, x=xs, y=ys, error_y=stds, error_y_minus=stds, log_x=True
        # )

        # print(7777777, data.df.head(11))
        return fig  # , html.Div(dcc.Graph(figure=fig), id="func_plot")

    @app.callback(
        Output("func_test_drop", "options"),
        Input("func_DATA_DROP", "value"),
    )
    def set_comparison_options(test: str) -> list[dict[str, str]]:
        """Set the test option by given experiemnt name"""
        print(54321, list(datasets[test].test_names))
        return make_list_of_dicts(list(datasets[test].test_names))

    @app.callback(
        Output("func_test_drop", "value"),
        Input("func_test_drop", "options"),
    )
    def select_test_value(options: list[dict[str, str]]) -> str:
        """Select first test as default value after changing dataset"""
        print(123456, options)
        return options[0]["value"]

    @app.callback(
        Output("func_plot", "figure"),
        Input("func_DATA_DROP", "value"),
        Input("func_test_drop", "value"),
        Input("func_metric_drop", "value"),
    )
    def update_graph(datadset_id: str, test: str, metric: str):
        """Update rendering of data points upon changing x-value of vertical dashed lines."""
        print(1111, datadset_id, test, metric)
        selected_data: Data = datasets[datadset_id]
        print(22222)
        filtered_data = selected_data.filter(test, metric)
        print(33333)
        return draw_line(filtered_data)

    dataset_names = list(datasets.keys())
    first_dataset = datasets[dataset_names[0]]
    default_comparison = list(first_dataset.test_names)[0]
    return html.Div(
        children=[
            html.P("blah"),
            html.H6("Dataset"),
            dcc.Dropdown(
                id="func_DATA_DROP",
                options=dataset_names,
                value=dataset_names[0],
                multi=False,
            ),
            html.H6("test"),
            dcc.Dropdown(
                id="func_test_drop",
                multi=False,
                value=default_comparison,
            ),
            html.H6("metric"),
            dcc.Dropdown(
                id="func_metric_drop",
                options=first_dataset.metrics,
                multi=False,
                value=first_dataset.metrics[0],
            ),
            html.Div(
                dcc.Graph(
                    id="func_plot",
                    figure=draw_line(
                        first_dataset.filter(
                            default_comparison, first_dataset.metrics[0]
                        ),
                    ),
                ),
            ),
        ],
    )
