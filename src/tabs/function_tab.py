from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from src.read_files import Data
from src.helpers import make_list_of_dicts
import plotly.graph_objects as go
import plotly.express as px


KEY = "function"


CB_COLOURS = [
    "#377eb8",
    "#ff7f00",
    "#4daf4a",
    "#f781bf",
    "#a65628",
    "#984ea3",
    "#999999",
    "#e41a1c",
    "#dede00",
]


def render(app: Dash, dataset: Data, ids2, params) -> html.Div:
    def draw_line(filtered_data: Data) -> html.Div:
        """Draws a box and wisker of the CPM data for each set of replicates for eact
        test and overlays the respective FDR value"""

        fig1 = px.scatter(
            data_frame=filtered_data.plot_data,
            x=filtered_data.drug,
            y=filtered_data.metric,
            error_y="stds",
            error_y_minus="stds",
            text="arrhythmia_counts",
            size_max=60,
            color="name & condition",
        )

        fig1.update_traces(textposition="top left")
        # figs.append(fig1.data[0])
        # print(1111111,figs)

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
        # fig3 = go.Figure(data=figs)
        # fig3.update_xaxes(type="log")
        # fig3.update_yaxes(range=[0, 1.5])
        fig1.update_xaxes(type="log")
        max_y = max(filtered_data.plot_data[filtered_data.metric])
        max_y = max_y if max_y > 1.5 else 1.5

        fig1.update_yaxes(range=[-0.2, max_y + 0.2])

        return fig1

    @app.callback(Output("func_DATA_DROP", "options"), Input("func_drug_drop", "value"))
    def select_gene_value(drug: str) -> str:
        """Select possible datasets after changing drug"""
        return make_list_of_dicts(dataset.possible_dataset_names(drug))

    @app.callback(
        Output("func_condition_DROP", "options"), Input("func_drug_drop", "value")
    )
    def select_gene_value(drug: str) -> str:
        """Select possible conditions after changing drug"""
        return make_list_of_dicts(dataset.possible_condition_names(drug))

    @app.callback(Output("func_DATA_DROP", "value"), Input("func_DATA_DROP", "options"))
    def select_gene_value(options: list[dict[str, str]]) -> str:
        """Select first gene as default value"""
        return [option['value'] for option in options]

    @app.callback(
        Output("func_condition_DROP", "value"), Input("func_condition_DROP", "options")
    )
    def select_gene_value(options: list[dict[str, str]]) -> str:
        """Select first gene as default value"""
        return [option['value'] for option in options]

    @app.callback(
        Output("func_plot", "figure"),
        Input("func_drug_drop", "value"),
        Input("func_metric_drop", "value"),
        Input("func_DATA_DROP", "value"),
        Input("func_condition_DROP", "value"),
    )
    def update_graph(
        test: str, metric: str, datasets: list[str], conditions: list[str]
    ):
        """Update rendering of data points upon changing x-value of vertical dashed lines."""
        filtered_data: Data = dataset.filter(test, conditions, datasets, metric)

        return draw_line(filtered_data)

    default_drug = list(dataset.test_names)[0]
    default_conditions = dataset.possible_condition_names(default_drug)
    default_datasets = dataset.possible_dataset_names(default_drug)
    filtered_dataset = dataset.filter(
        default_drug,
        default_conditions,
        default_datasets,
        dataset.metrics[0],
    )

    return html.Div(
        children=[
            html.P("blah bla"),
            html.H6("Dataset"),
            dcc.Dropdown(
                id="func_DATA_DROP",
                options=default_datasets,
                value=default_datasets,
                multi=True,
            ),
            html.H6("Condition"),
            dcc.Dropdown(
                id="func_condition_DROP",
                options=default_conditions,
                value=default_conditions,
                multi=True,
            ),
            html.H6("Drug"),
            dcc.Dropdown(
                id="func_drug_drop",
                options=list(dataset.test_names),
                multi=False,
                value=default_drug,
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
                    figure=draw_line(filtered_dataset),
                ),
            ),
        ],
    )
