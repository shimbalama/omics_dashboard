import numpy as np
from dash import dcc, html
import plotly.express as px
from src.read_files import Data
from dataclasses import dataclass

@dataclass
class Params:
    """Parameters for box plotting"""

    X: str
    COLOUR: str | None = None
    LOG: bool = False
    Y: str | None = None


@dataclass
class IDs:
    """Parameters for box div ids"""

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


# class TestIDs(IDs):
# class VolcanoIDs(IDs):



def make_list_of_dicts(values: list[str]) -> list[dict[str, str]]:
    """Convert a list of strs into a list where those strings are values in dicts
    against keys label and value, for use in callbacks"""
    return [{"label": val, "value": val} for val in sorted(set(values))]


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


def draw_box_chart(data: Data, y_gene: str, params: type, plot_id: str) -> html.Div:
    """Draws a box and wisker of the CPM data for each set of replicates for eact
    test and overlays the respective FDR value"""
    df = data.pandas_df
    fig = px.box(
        df,
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
    if "test" in set(df.columns):
        FDRs = data.df_FDR[y_gene].to_dict()
        median_CPMs = {
            test: np.median(df.loc[df["test"] == test, y_gene])
            for test in np.unique(df["test"])
        }
        unique_comparisons = df["test"].unique()
        y_range = get_y_range(len(unique_comparisons))
        if data.point_of_reference in unique_comparisons:
            point_of_reference_index = [
                i
                for i, test in enumerate(unique_comparisons)
                if test == data.point_of_reference
            ]
            assert len(point_of_reference_index) == 1
            for i, test in enumerate(unique_comparisons):
                if test == data.point_of_reference:
                    assert i == point_of_reference_index[0]
                    continue
                FDR = FDRs.get(test, 0.0)
                y = median_CPMs[test]
                fig.add_annotation(
                    x=test,
                    y=y,
                    text=f"{FDR:.1e}",
                    yshift=10,
                    showarrow=False,
                )
                fig = add_FDR_brackets(
                    fig, FDR, i, [i, point_of_reference_index[0]], y_range
                )
            fig.update_layout(margin=dict(t=i * 33))

    return html.Div(dcc.Graph(figure=fig), id=plot_id)
