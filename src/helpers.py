import numpy as np
from dash import dcc, html
import plotly.express as px
from src.read_files import Data
from dataclasses import dataclass
from functools import partial
import pandas as pd


@dataclass
class Params:
    """Parameters for box plotting"""

    X: str
    COLOUR: str | None = None
    LOG: bool = False
    Y: str | None = None


@dataclass
class IDs:
    """div ids"""

    name: str

    @property
    def data_drop(self) -> str:
        return f"{self.name}_dataset-dropdown"

    @property
    def gene_drop(self) -> str:
        return f"{self.name}_gene-dropdown"

    @property
    def plot(self) -> str:
        return f"{self.name}_plot"

    @property
    def tests_drop(self) -> str:
        return f"{self.name}_tests-dropdown"

    @property
    def select_all(self) -> str:
        return f"{self.name}_select-all"

    @property
    def slider1(self) -> str:
        return f"{self.name}_slider1"

    @property
    def slider2(self) -> str:
        return f"{self.name}_slider2"


def make_list_of_dicts(values: list[str]) -> list[dict[str, str]]:
    """Convert a list of strs into a list where those strings are values in dicts
    against keys label and value, for use in callbacks"""
    return [{"label": val, "value": val} for val in sorted(set(values))]


def get_y_range(number_of_comparisons: int, interline: float = 0.03) -> np.ndarray:
    """Specify in what y_range to plot for each pair of columns"""
    y_range = np.zeros([number_of_comparisons, 2])
    for i in range(number_of_comparisons):
        y_range[i] = [0.93 + i * interline, 0.94 + i * interline]
    return y_range


def draw_box_chart(data: Data, y_gene: str, params: type, plot_id: str) -> html.Div:
    """Draws a box and wisker of the CPM data for each set of replicates for eact
    test and overlays the respective FDR value"""

    df: pd.DataFrame = data.plot_df
    fig = px.box(
        df,
        x=params.X,
        y=y_gene,
        points="all",
        width=1000,
        height=700,
        color=params.COLOUR,  # TODO
        log_y=params.LOG,
        title=f"Boxplot for CPMs",
        labels={y_gene: "CPM"},
        facet_row_spacing=0.75,
    )

    fig = make_brackets(fig, data)

    return html.Div(dcc.Graph(figure=fig), id=plot_id)


def make_brackets(fig, data):
    def get_point_of_reference():
        point_of_reference_index = [
            i
            for i, test in enumerate(data.test_names)
            if test == data.point_of_reference
        ]
        assert len(point_of_reference_index) == 1
        return point_of_reference_index[0]

    # https://community.plotly.com/t/grouped-bar-charts-with-corresponding-line-chart/19562/4
    def bracket_subplots(fig, data):
        point_of_reference = get_point_of_reference()
        # add_FDR_brackets_partial = partial(
        #     add_FDR_brackets, point_of_reference, y_range, 1.01, "black"
        # )
        for i, test in enumerate(data.test_names):
            if test == data.point_of_reference:
                assert i == point_of_reference
                continue
            FDR = data.get_FDR(test)
            y = data.get_median_CPMs(test)

            fig.add_annotation(
                x=test,
                y=y,
                text=f"{FDR:.1e}",
                yshift=10,
                showarrow=False,
            )
            bracket = Bracket(
                x_POR=point_of_reference, x_other=i, y_range=y_range, y_pos=i, FDR=FDR
            )
            fig = add_FDR_brackets(fig, bracket)
        fig.update_layout(margin=dict(t=i * 33))

        return fig

    def bracket_subplots2(i, fig, data, prot_pos):
        point_of_reference = get_point_of_reference()

        number_plots_per_subplot = len(data.test_names)
        bargap = 0.2
        gap_width = 1 - bargap
        bar_width = (
            gap_width / number_plots_per_subplot
        )  # boxes_per_sub_plot make var TODO

        gpEdges = i - gap_width / 2

        sub_poses = [
            gpEdges + bar_width / 2 + bar_width * n
            for n in range(number_plots_per_subplot)
        ]

        for j, test in enumerate(data.test_names):
            if test == data.point_of_reference:
                assert j == point_of_reference
                continue

            bracket = Bracket(
                x_POR=sub_poses[point_of_reference],
                x_other=sub_poses[j],
                y_range=y_range,
                y_pos=j,
                FDR=data.get_FDR(test, prot_pos),
            )
            fig = add_FDR_brackets(fig, bracket)
        fig.update_layout(margin=dict(t=i * 33))

        return fig

    if data.point_of_reference in data.test_names:
        y_range = get_y_range(len(data.test_names))
        if isinstance(data.df_FDR, pd.DataFrame):
            for i, prot_pos in enumerate(data.df_FDR.columns):
                if prot_pos != "test":
                    fig = bracket_subplots2(i, fig, data, prot_pos)
        else:
            fig = bracket_subplots(fig, data)

    return fig


@dataclass(slots=True, frozen=True)
class Bracket:
    """Bracket for FDR notation"""

    x_POR: float
    x_other: float
    y_range: np.ndarray
    y_pos: int
    FDR: float
    text_height: float = 1.01
    color: str = "black"

    @property
    def stars(self):
        if self.FDR >= 0.05:
            return "ns"
        elif self.FDR >= 0.01:
            return "*"
        elif self.FDR >= 0.001:
            return "**"
        else:
            return "***"

    @property
    def y0(self):
        return self.y_range[self.y_pos][0]

    @property
    def y1(self):
        return self.y_range[self.y_pos][1]

    @property
    def star_x(self):
        return (self.x_POR + self.x_other) / 2


def add_FDR_brackets(fig, bracket: Bracket):
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

    # Vertical line
    # Vertical line at value
    fig.add_shape(
        type="line",
        xref="x",
        yref="y" + " domain",
        x0=bracket.x_POR,
        y0=bracket.y0,
        x1=bracket.x_POR,
        y1=bracket.y1,
        line=dict(
            color=bracket.color,
            width=2,
        ),
    )
    # Horizontal line
    fig.add_shape(
        type="line",
        xref="x",
        yref="y" + " domain",
        x0=bracket.x_POR,
        y0=bracket.y1,
        x1=bracket.x_other,
        y1=bracket.y1,
        line=dict(
            color=bracket.color,
            width=2,
        ),
    )
    # Vertical line at POR
    fig.add_shape(
        type="line",
        xref="x",
        yref="y" + " domain",
        x0=bracket.x_other,
        y0=bracket.y0,
        x1=bracket.x_other,
        y1=bracket.y1,
        line=dict(
            color=bracket.color,
            width=2,
        ),
    )
    ## add text at the correct x, y coordinates
    ## for bars, there is a direct mapping from the bar number to 0, 1, 2...
    fig.add_annotation(
        dict(
            font=dict(color=bracket.color, size=14),
            x=bracket.star_x,
            y=bracket.y1 * bracket.text_height,
            showarrow=False,
            text=bracket.stars,
            textangle=0,
            xref="x",
            yref="y" + " domain",
        )
    )

    return fig


# def add_FDR_brackets(bracket: Bracket):
#     """Adds notations giving the significance level between two box plot data (t-test two-sided test)

#     Parameters:
#     ----------
#     fig: figure
#         Plotly boxplot figure.
#     FDR: float
#         False Discovery Rate (FDR) value for the test.
#     i: int
#         Index of the y_range to use for the box plot data.
#     column_pair: list
#         List of two column indices to compare.
#         e.g.: [0, 1] compares column 0 with column 1.
#     y_range: list
#         List of y-ranges for the box plot data.
#     text_height: float, optional
#         Height of the text annotation above the box plot, default is 1.01.
#     color: str, optional
#         Color of the lines and text annotation, default is "black".

#     Returns:
#     -------
#     fig: figure
#         Figure with the added notation.
#     """

#     # Vertical line
#     fig.add_shape(
#         type="line",
#         xref="x",
#         yref="y" + " domain",
#         x0=i,
#         y0=y_range[i][0],
#         x1=i,
#         y1=y_range[i][1],
#         line=dict(
#             color=color,
#             width=2,
#         ),
#     )
#     # Horizontal line
#     fig.add_shape(
#         type="line",
#         xref="x",
#         yref="y" + " domain",
#         x0=i,
#         y0=y_range[i][1],
#         x1=point_of_reference_index,
#         y1=y_range[i][1],
#         line=dict(
#             color=color,
#             width=2,
#         ),
#     )
#     # Vertical line
#     fig.add_shape(
#         type="line",
#         xref="x",
#         yref="y" + " domain",
#         x0=point_of_reference_index,
#         y0=y_range[i][0],
#         x1=point_of_reference_index,
#         y1=y_range[i][1],
#         line=dict(
#             color=color,
#             width=2,
#         ),
#     )
#     ## add text at the correct x, y coordinates
#     ## for bars, there is a direct mapping from the bar number to 0, 1, 2...
#     fig.add_annotation(
#         dict(
#             font=dict(color=color, size=14),
#             x=(i + point_of_reference_index) / 2,
#             y=y_range[i][1] * text_height,
#             showarrow=False,
#             text=symbol,
#             textangle=0,
#             xref="x",
#             yref="y" + " domain",
#         )
#     )

#     return fig


def rubbish(name: str) -> bool:
    """checks if file name starts with . or ~"""
    return name.startswith((".", "~"))
