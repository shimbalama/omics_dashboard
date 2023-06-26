from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from src.read_files import Data
from src.helpers import make_list_of_dicts
import plotly.express as px


KEY = "function"

# np.seterr(divide = 'ignore')
# def sigmoid1(x, EC50, hill_slope):
#     return 1 / (1 + np.exp(np.abs(EC50 - x) * hill_slope))

# def sigmoid2(x, EC50, Hill, slope):
#     return 1 / (1 + (EC50 / x) ** Hill) ** slope

# def sigmoid3(x,b,c,d,e):
#     '''This function is basically a copy of the LL.4 function from the R drc package with
#     - b: hill slope
#     - c: min response
#     - d: max response
#     - e: EC50'''
#     return(c+(d-c)/(1+np.exp(b*(np.log(x)-np.log(e)))))


def render(app: Dash, dataset: Data, ids2, params) -> html.Div:
    def draw_line(filtered_data: Data) -> html.Div:
        """Draws a box and wisker of the CPM data for each set of replicates for eact
        test and overlays the respective FDR value"""
        if filtered_data.plot_data.empty:
            return html.Div("No data selected")
        filtered_data.plot_data.to_csv("~/Downloads/plot_data.csv")

        fig = px.scatter(
            data_frame=filtered_data.plot_data,
            x=filtered_data.drug,
            y=filtered_data.metric,
            error_y="stds",
            error_y_minus="stds",
            text="arrhythmia_counts",
            size_max=60,
            color="name & condition",
        )

        poses = ["top left", "top right", "bottom left", "bottom right"] * 3

        for i, _ in enumerate(filtered_data.plot_data.groupby("name & condition")):
            colour = fig.data[i].marker.color
            fig.data[i].update(
                textfont_color=colour,
                textposition=poses[i],
                mode="lines+markers+text",
            )
            fig.data[i].line.color = colour
            fig.data[i].line.shape = "spline"
        fig.update_xaxes(type="log")  # range=[-1.5, max_y + 1])

        # for name, cond_df in filtered_data.plot_data.groupby("name & condition"):

        #     print(33333, cond_df)
        #     cond_df.replace(0, 0.001, inplace=True)
        #     print(444444, cond_df)
        #     colour = fig.data[i].marker.color
        #     x = cond_df[filtered_data.drug].to_numpy()
        #     y = cond_df[filtered_data.metric].to_numpy()
        #     ec50_in_range = False

        #         params, _ = curve_fit(sigmoid1, x, y)
        #         EC50 = params[0]
        #         if min(x) < EC50 < max(x):
        #             ec50_in_range = True
        #             sigmoid = partial(sigmoid1, EC50=EC50, hill_slope=params[1])

        #         params, _ = curve_fit(sigmoid2, x, y)
        #         EC50 = params[0]
        #         if min(x) < EC50 < max(x):
        #             ec50_in_range = True
        #             sigmoid = partial(sigmoid2, EC50=EC50, Hill=params[1], slope=params[2])
        #     print(11111111, EC50, ec50_in_range)
        #     # Evaluate sigmoid curve at desired concentration points
        #     curve_concentrations = np.logspace(
        #         np.log10(x.min()), np.log10(x.max()), 100
        #     )
        #     curve_response = sigmoid(curve_concentrations)
        #     print(9999999999,curve_concentrations, curve_response)
        #     # Add sigmoid curve to the plot

        #     fig.add_trace(px.line(x=curve_concentrations, y=curve_response).data[0])
        #     fig.data[0].line.color = colour

        #     if not ec50_in_range:
        #         EC50 = 0.1

        #     max_stds = max(cond_df["stds"]) TODO add this back in from provided EC50s don't del
        #     fig.add_shape(
        #         type="line",
        #         x0=EC50,
        #         y0=min(curve_response) - max_stds,
        #         x1=EC50,
        #         y1=max(curve_response) + max_stds,
        #         line=dict(color=colour, dash="dash"),
        #     )
        #     fig.add_annotation(
        #         x=np.log10(EC50),
        #         y=max(curve_response) - (i/5)-0.1,
        #         text=f"EC50 = {EC50:.2f}",
        #         showarrow=True,
        #         arrowhead=1,
        #         arrowcolor=colour,
        #         bordercolor=colour,
        #         ax=30,
        #         ay=-40,
        #     )

        #         # fig.add_vline(
        #         #     x=EC50,
        #         #     line_width=3,
        #         #     line_dash="dash",
        #         #     line_color=colour,
        #         #     name=f"EC50: {EC50:.2f}",
        #         # )

        #     i += 1

        # max_y = max(filtered_data.plot_data[filtered_data.metric])
        # max_y = max_y if max_y > 1.5 else 1.5
        # max_x = max(filtered_data.plot_data[filtered_data.drug])
        # max_x = max_x if max_x > 1.5 else 1.5
        # fig.update_yaxes(range=[-0.2, max_y + 0.2])
        # fig.update_xaxes()
        # fig.data[i].update(
        #     textfont_color=fig1.data[i].marker.color,
        #     textposition=poses[i],
        #     mode="text",
        #     showlegend=False,
        # )

        # Create a DataFrame for plotting
        # df = pd.DataFrame(
        #     [xs, ys, names, ec50s],
        # ).T
        # df.columns = [filtered_data.drug, filtered_data.metric, "name", "EC50"]

        # # Plot the data and the fitted sigmoid curve
        # fig2 = px.line(
        #     df,
        #     x=filtered_data.drug,
        #     y=filtered_data.metric,
        #     line_shape="spline",
        #     color="name",
        # )
        # print(1111111, fig1.data[-1], fig1.data[-1].marker.color, len(fig1.data))
        # new_trace = fig1.data[-1]

        return fig

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
        return [option["value"] for option in options]

    @app.callback(
        Output("func_condition_DROP", "value"), Input("func_condition_DROP", "options")
    )
    def select_gene_value(options: list[dict[str, str]]) -> str:
        """Select first gene as default value"""
        return [option["value"] for option in options]

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
