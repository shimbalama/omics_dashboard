from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from src.read_files import Data
import numpy as np
from src.helpers import make_list_of_dicts
import plotly.graph_objects as go
import plotly.express as px
from scipy.special import expit
from scipy.optimize import curve_fit

KEY = "function"

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


def render(app: Dash, dataset: Data, ids2, params) -> html.Div:
    def draw_line(filtered_data: Data) -> html.Div:
        """Draws a box and wisker of the CPM data for each set of replicates for eact
        test and overlays the respective FDR value"""
        figs = []
        for name, df, dose, arrhythmia in filtered_data.filtered_data:
            print("name, df, dose, arrhythmia", name, df, dose, arrhythmia)

            datas = df.describe().to_dict()
            dose["0"] = 0

            xs = [dose[str(time_point)] for time_point in list(datas.keys())]
            ys = [v.get("50%") for v in datas.values()]
            # # https://www.statology.org/sigmoid-function-python/
            # # https://datagy.io/sigmoid-function-python/
            stds = [v.get("std") for v in datas.values()]
            arrhythmia_counts: list[str] = []
            for x in xs:
                hhh = arrhythmia.query("Dose == @x")
                dd: dict[str, int] = dict(zip(hhh["Arrhythmia"], hhh["count"]))
                total = sum(dd.values())
                non_arrhythmic: int = total - dd.get("Y", 0)
                arrhythmia_counts.append(f"{str(non_arrhythmic)}/{str(total)}")

            fig1 = px.scatter(
                df,
                x=xs,
                y=ys,
                error_y=stds,
                error_y_minus=stds,
                text=arrhythmia_counts,
                size_max=60,
            )
            fig1.update_traces(textposition="top left")
            figs.append(fig1)
        fig1 = figs[0]

        # y = data.df[list(range(1, len(data.df.columns)))].to_numpy()
        # print(y,y.shape, xs, y.shape[0],len(xs), sep='\n')
        # #assert y.shape[0] == len(xs)-1
        # x = np.array([[conc]*y.shape[1] for conc in xs])
        # x, y = x.flatten(), y.flatten()

        # def sigmoid(x, EC50, Hill_slope):
        #     return 1 / (1 + (x / EC50) ** (-Hill_slope))

        # # Fit the sigmoidal function to the data
        # popt, pcov = curve_fit(sigmoid, x, y)
        # EC50 = popt[0]
        # print('EC50', EC50)

        # # Plot the data and the fitted curve
        # x_fit = np.linspace(min(x), max(x), len(x))
        # y_fit = sigmoid(x_fit, *popt)
        # fig2 = px.line(data.df, x=x_fit, y=y_fit)
        # fig3 = go.Figure(data=fig1.data + fig2.data)

        # fig3.update_xaxes(type="log")
        # fig3.update_yaxes(range=[0, 1.5])
        fig1.update_xaxes(type="log")
        max_y = max(ys) if max(ys) > 1.5 else 1.5

        fig1.update_yaxes(range=[-0.1, max_y])

        return fig1

   

    @app.callback(
        Output("func_plot", "figure"),
        Input("func_test_drop", "value"),
        Input("func_metric_drop", "value"),
        Input("func_DATA_DROP", "value"),
        Input("func_condition_DROP", "value"),
    )
    def update_graph(test: str, metric: str, datasets, conditions):
        """Update rendering of data points upon changing x-value of vertical dashed lines."""
        filtered_data: Data = dataset.filter(test, metric, datasets, conditions)
        return draw_line(filtered_data)

    default_comparison = list(dataset.test_names)[0]
    return html.Div(
        children=[
            html.P("blah bla"),
            html.H6("Dataset"),
            dcc.Dropdown(
                id="func_DATA_DROP",
                options=list(dataset.dataset_names),
                value=list(dataset.dataset_names),
                multi=True,
            ),
            html.H6("Condition"),
            dcc.Dropdown(
                id="func_condition_DROP",
                options=list(dataset.condition_names),
                value=list(dataset.condition_names),
                multi=True,
            ),
            html.H6("Drug"),
            dcc.Dropdown(
                id="func_drug_drop",
                options=list(dataset.test_names),
                multi=False,
                value=default_comparison,
            ),
            html.H6("metric"),
            dcc.Dropdown(
                id="func_metric_drop",
                options=dataset.metrics,
                multi=False,
                value=dataset.metrics[0],
            ),
            html.Div(
                dcc.Graph(
                    id="func_plot",
                    figure=draw_line(
                        dataset.filter(
                            default_comparison,
                            dataset.metrics[0],
                            list(dataset.dataset_names),
                            list(dataset.condition_names),
                        ),
                    ),
                ),
            ),
        ],
    )
