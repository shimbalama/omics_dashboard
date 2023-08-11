from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from src.read_files import Data
from src.helpers import make_list_of_dicts
import plotly.express as px
import numpy as np
import pandas as pd
import plotly.graph_objects as go

KEY = "function"


def render(app: Dash, dataset: Data, ids2, params) -> go.Figure:
    def draw_line(filtered_data: Data) -> html.Div:
        """Draws a box and wisker of the CPM data for each set of replicates for eact
        test and overlays the respective FDR value"""

        if filtered_data.plot_df.empty:
            # Data for the letters
            letters = {
                "N": [(1, 4), (1, 1), (2, 4), (2, 1)],
                "O": [(3, 1), (3, 4), (4, 4), (4, 1), (3, 1)],
                "D": [(6, 1), (6, 4), (7, 4), (7, 2), (6, 1)],
                "A": [(8, 1), (8, 4), (9, 4), (9, 1), (8, 3), (9, 3)],
                "T": [(10, 4), (11, 4), (10.5, 4), (10.5, 1)],
                "A2": [(12, 1), (12, 4), (13, 4), (13, 1), (12, 3), (13, 3)],
            }

            # Create empty figure
            fig = go.Figure()

            # Add lines for each letter
            for letter, coords in letters.items():
                x, y = zip(*coords)
                fig.add_trace(go.Scatter(x=x, y=y, mode="lines+markers"))

            # Hide axes for cleaner look
            fig.update_xaxes(showticklabels=False, zeroline=False)
            fig.update_yaxes(showticklabels=False, zeroline=False)

            # Set layout
            fig.update_layout(width=666, height=666)

            return fig

        fig = px.scatter(
            data_frame=filtered_data.plot_df,
            x=filtered_data.drug,
            y=filtered_data.metric,
            error_y="stds",
            error_y_minus="stds",
            text="arrhythmia_counts",
            size_max=60,
            color="name & condition",
        )

        poses = ["top left", "top right", "bottom left", "bottom right"] * 3

        colors = {}
        for i, name_condition in enumerate(
            filtered_data.plot_df.groupby("name & condition")
        ):
            colour = fig.data[i].marker.color
            colors[name_condition[0].split("_")[1][0]] = colour
            fig.data[i].update(
                textfont_color=colour,
                textposition=poses[i],
                mode="markers+text",
            )  # lines+
            # fig.data[i].line.color = colour
            # fig.data[i].line.shape = "spline"
        fig.update_xaxes(type="log")

        ec50 = filtered_data.ec50.query(
            "Drug == @filtered_data.drug & Metric == @filtered_data.metric"
        )
        assert ec50.shape[0] <= len(colors)
        for experiment, EC50 in zip(list(ec50["Experiment"]), list(ec50["EC50 (uM)"])):
            colour = colors.get(str(experiment), "black")
            # EC50 = ec50.loc[experiment][]
            if EC50 > 0:
                fig.add_shape(
                    type="line",
                    x0=EC50,
                    y0=min(filtered_data.plot_df[filtered_data.metric]) - 0.1,
                    x1=EC50,
                    y1=max(filtered_data.plot_df[filtered_data.metric]) + 0.1,
                    line=dict(color=colour, dash="dash"),
                )
                x = np.log10(EC50)
            else:
                x = 0.1

            fig.add_annotation(
                x=x,
                y=max(filtered_data.plot_df[filtered_data.metric]) - 0.1,
                text=f"EC50 = {EC50:.2f}",
                showarrow=True,
                arrowhead=1,
                arrowcolor=colour,
                bordercolor=colour,
                ax=30,
                ay=-40,
            )

        # Set layout
        fig.update_layout(
            width=666,
            height=666,
            legend=dict(
                orientation="h",
                entrywidth=333,
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1,
            ),
        )
        fig.update_layout(modebar_orientation='v')
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
