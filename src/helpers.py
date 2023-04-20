import numpy as np
from time import time
import pandas as pd
from dash import dcc, html
import plotly.express as px
from src.read_files import Data
from .components import ids
from dash.dependencies import Input, Output

from dataclasses import dataclass


@dataclass
class Params:
    """Parameters for box plotting"""

    DIV_ID: str
    X: str
    COLOUR: str | None = None
    LOG: bool = False
    Y: str | None = None


def make_list_of_dicts(values: list[str]) -> list[dict[str, str]]:
    """Convert a list of strs into a list where those strings are values in dicts
    against keys label and value, for use in callbacks"""
    return [{"label": val, "value": val} for val in sorted(values)]


def get_y_range(number_of_comparisons, interline=0.03):
    # Specify in what y_range to plot for each pair of columns
    y_range = np.zeros([number_of_comparisons, 2])
    for i in range(number_of_comparisons):
        y_range[i] = [0.93 + i * interline, 0.94 + i * interline]
    return y_range


def add_FDR_brackets(
    fig, FDR, i, column_pair, y_range, text_height=1.01, color="black"
):
    """Adds notations giving the significance level between two box plot data (t-test two-sided test)

    Parameters:
    ----------
    fig: figure
        Plotly boxplot figure.
    FDR: float
        False Discovery Rate (FDR) value for the test.
    i: int
        Index of the y_range to use for the box plot data.
    column_pair: list
        List of two column indices to compare.
        e.g.: [0, 1] compares column 0 with column 1.
    y_range: list
        List of y-ranges for the box plot data.
    text_height: float, optional
        Height of the text annotation above the box plot, default is 1.01.
    color: str, optional
        Color of the lines and text annotation, default is "black".

    Returns:
    -------
    fig: figure
        Figure with the added notation.
    """

    if FDR >= 0.05:
        symbol = "ns"
    elif FDR >= 0.01:
        symbol = "*"
    elif FDR >= 0.001:
        symbol = "**"
    else:
        symbol = "***"
    # Vertical line
    fig.add_shape(
        type="line",
        xref="x",
        yref="y" + " domain",
        x0=column_pair[0],
        y0=y_range[i][0],
        x1=column_pair[0],
        y1=y_range[i][1],
        line=dict(
            color=color,
            width=2,
        ),
    )
    # Horizontal line
    fig.add_shape(
        type="line",
        xref="x",
        yref="y" + " domain",
        x0=column_pair[0],
        y0=y_range[i][1],
        x1=column_pair[1],
        y1=y_range[i][1],
        line=dict(
            color=color,
            width=2,
        ),
    )
    # Vertical line
    fig.add_shape(
        type="line",
        xref="x",
        yref="y" + " domain",
        x0=column_pair[1],
        y0=y_range[i][0],
        x1=column_pair[1],
        y1=y_range[i][1],
        line=dict(
            color=color,
            width=2,
        ),
    )
    ## add text at the correct x, y coordinates
    ## for bars, there is a direct mapping from the bar number to 0, 1, 2...
    fig.add_annotation(
        dict(
            font=dict(color=color, size=14),
            x=(column_pair[0] + column_pair[1]) / 2,
            y=y_range[i][1] * text_height,
            showarrow=False,
            text=symbol,
            textangle=0,
            xref="x",
            yref="y" + " domain",
        )
    )

    return fig


def rubbish(name: str) -> bool:
    return name.startswith((".", "~"))


def draw_box_chart(data: Data, y_gene: str, params: type) -> html.Div:
    """Draws a box and wisker of the CPM data for each set of replicates for eact
    test and overlays the respective FDR value"""

    fig = px.box(
        data.df,
        x=params.X,
        y=y_gene,
        points="all",
        width=1300,
        height=900,
        color=params.COLOUR,
        log_y=params.LOG,
        title=f"Boxplot for CPMs",
        labels={y_gene: "CPM"},
        facet_row_spacing=0.75,
    )
    #s = time()

    FDRs = data.df_FDR[y_gene].to_dict()
    #print(f"time1:{s-time()}")
    # Get the unique values of the 'test' column
    tests = np.unique(data.df["test"])

    # Initialize a dictionary to store the medians
    CPMs = {
        test: np.median(data.df.loc[data.df["test"] == test, y_gene]) for test in tests
    }

    
    #CPMs = data.df[["test", y_gene]].groupby("test").agg("median").to_dict()
    #print(f"time2:{s-time()}")
    # print(1111, data.df.columns)
    unique_comparisons = data.df["test"].unique()
    #print(f"time3:{s-time()}")
    y_range = get_y_range(len(unique_comparisons))
    #print(f"time4:{s-time()}")
    # print(1212121212, data.point_of_reference, unique_comparisons)
    if data.point_of_reference in unique_comparisons:
        #print(f"time5:{s-time()}")
        # print(2222222)
        point_of_reference_index = [
            i
            for i, test in enumerate(unique_comparisons)
            if test == data.point_of_reference
        ]
        #print(f"time6:{s-time()}")
        assert len(point_of_reference_index) == 1
        #print(f"time7:{s-time()}")
        for i, test in enumerate(unique_comparisons):
            #print(f"time8:{s-time()}")
            # print(11111111111, i, test, y_gene, data.point_of_reference)
            if test == data.point_of_reference:
                assert i == point_of_reference_index[0]
                continue
            #print(f"time9:{s-time()}")
            FDR = FDRs.get(test, 0.0)
            #print(f"time10:{s-time()}")
            y = CPMs[test]
            #print(f"time11:{s-time()}")
            fig.add_annotation(
                x=test,
                y=y,
                text=f"{FDR:.1e}",
                yshift=10,
                showarrow=False,
            )
            #print(f"time12:{s-time()}")
            # print(i, test, y, FDR)
            fig = add_FDR_brackets(
                fig, FDR, i, [i, point_of_reference_index[0]], y_range
            )
            #print(f"time13:{s-time()}")
        fig.update_layout(margin=dict(t=i * 33))
    # if "test" in set(data.df.columns):
    #     unique_comparisons = data.df["test"].unique()
    #     y_range = get_y_range(len(unique_comparisons))
    #     if data.point_of_reference in unique_comparisons:
    #         point_of_reference_index = list(unique_comparisons).index(
    #             data.point_of_reference
    #         )

    #         # Compute FDR and median for all tests at once using vectorized operations
    #         FDR_values = data.df_FDR.loc[
    #             data.df_FDR["test"].isin(unique_comparisons), y_gene
    #         ]
    #         print(FDR_values)
    #         median_values = (
    #             data.df.loc[data.df["test"].isin(unique_comparisons)]
    #             .groupby("test")[y_gene]
    #             .median()
    #         )

    #         for i, test in enumerate(unique_comparisons):
    #             if test == data.point_of_reference:
    #                 continue

    #             # Get the FDR and median values for this test from the precomputed arrays
    #             if FDR_values.empty:
    #                 print(98788766)  # wtf???
    #                 continue
    #             FDR = FDR_values.loc[data.df_FDR["test"] == test].values[0]
    #             y = median_values.loc[test]

    #             # Add the annotation and FDR bracket to the plot
    #             fig.add_annotation(
    #                 x=test,
    #                 y=y,
    #                 text=f"{FDR:.1e}",
    #                 yshift=10,
    #                 showarrow=False,
    #             )
    #             fig = add_FDR_brackets(
    #                 fig, FDR, i, [i, point_of_reference_index], y_range
    #             )

    #         # Update the layout with the new top margin
    #         fig.update_layout(margin=dict(t=i * 33))

    return html.Div(dcc.Graph(figure=fig), id=params.DIV_ID)


def gene_dropdown(app, cb_out: str, cb_in: str, data_set: Data):
    @app.callback(Output(cb_out, "options"), Input(cb_in, "value"))
    def set_gene_options(experiment: str) -> list[dict[str, str]]:
        """Populates the gene selection dropdown with options from teh given dataset"""
        return make_list_of_dicts(list(data_set[experiment].df.columns))


def gene_dropdown_default(app, cb_id: str):
    @app.callback(Output(cb_id, "value"), Input(cb_id, "options"))
    def select_gene_value(gene_options: list[dict[str, str]]) -> str:
        """Select first gene as default value"""
        return gene_options[0]["value"]


def box(app, cb_in: str, cb_in2: str, data_set: Data, params: type):
    @app.callback(
        Output(params.DIV_ID, "children"), Input(cb_in, "value"), Input(cb_in2, "value")
    )
    def update_box_chart(experiment: str, gene: str) -> html.Div:
        """Re draws a box and wisker of the CPM data for each set of replicates for eact
        test and overlays the respective FDR value"""
        selected_data = data_set[experiment]
        filtered: Data = selected_data.filter(gene)
        y_param = params.Y if params.Y else gene
        return draw_box_chart(filtered, y_param, params)
