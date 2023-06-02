from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from src.read_files import Data
import pandas as pd
from src.helpers import make_list_of_dicts
import plotly.graph_objects as go
import plotly.express as px
from scipy.special import expit

KEY = "function"


def render(app: Dash, datasets: dict[str, Data], ids2, params) -> html.Div:
    def draw_line(data: Data, test: str) -> html.Div:
        """Draws a box and wisker of the CPM data for each set of replicates for eact
        test and overlays the respective FDR value"""
        arrhythmia = data.arrhythmia.query('Drug == @test')

        print(1111111122, arrhythmia, sep='\n')
        
        
        # conc_with_arrhythmia = data.arrhythmia.get(test)
        # if conc_with_arrhythmia:
        #     print("conc_with_arrhythmia", conc_with_arrhythmia)
        datas = data.df.describe().to_dict()
        print('datas', datas)
        data.dose["0"] = 0

        xs = [data.dose[str(time_point)] for time_point in list(datas.keys())]
        ys = [v.get("50%") for v in datas.values()]
        # # https://www.statology.org/sigmoid-function-python/
        # # https://datagy.io/sigmoid-function-python/
        # ys = expit(xs)[::-1]
        stds = [v.get("std") for v in datas.values()]
        arrhythmia_counts: list[str] = []
        for x in xs:
            hhh = arrhythmia.query('Dose == @x')
            dd: dict[str, int] = dict(zip(hhh['Arrhythmia'], hhh['count']))
            print('ddddddd', dd)
            total = sum(dd.values())
            non_arrhythmic: int = total - dd.get('Y', 0)
            arrhythmia_counts.append(f'{str(non_arrhythmic)}/{str(total)}')
        print('arrhythmia_counts', arrhythmia_counts)
        # fig = go.Figure(
        #     data=go.Scatter(
        #         x=xs,
        #         y=ys,
        #         error_y=dict(
        #             type="data",  # value of error bar given in data coordinates
        #             array=stds,
        #             visible=True,
        #         ),
        #         mode="markers",
        #         text=arrhythmia_counts,
        #         textposition="top center"
        #     )
        # )
        # fig.update_xaxes(type="log")
        # fig.update_yaxes(range=[0, 1.5])

        fig = px.scatter(
            data.df, x=xs, y=ys, error_y=stds, error_y_minus=stds, text=arrhythmia_counts, log_x=True, size_max=60
        )
        fig.update_traces(textposition='middle left')
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
        selected_data: Data = datasets[datadset_id]
        filtered_data = selected_data.filter(test, metric)
        return draw_line(filtered_data, test)

    dataset_names = list(datasets.keys())
    first_dataset = datasets[dataset_names[0]]
    default_comparison = list(first_dataset.test_names)[0]
    print(33333, dataset_names, default_comparison, first_dataset)
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
                        default_comparison,
                    ),
                ),
            ),
        ],
    )
