import numpy as np
from dash import dcc, html
import plotly.express as px
from src.parse.read_files import Data
from dataclasses import dataclass, field
import pandas as pd
import plotly.graph_objects as go


@dataclass
class IDs:  # TODO more needed
    """div ids"""

    name: str

    @property
    def data_drop(self) -> str:
        return f"{self.name}_dataset-dropdown"

    @property
    def gene_drop(self) -> str:
        return f"{self.name}_gene-dropdown"
    
    @property
    def condition_drop(self) -> str:
        return f"{self.name}_cond-dropdown"
    
    @property
    def metric_drop(self) -> str:
        return f"{self.name}_metric-dropdown"

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
    def stats(self) -> str:
        return f"{self.name}_stats"

    @property
    def stats_out(self) -> str:
        return f"{self.name}_stats-out"

    @property
    def log(self) -> str:
        return f"{self.name}_log"

    @property
    def log_out(self) -> str:
        return f"{self.name}_log-out"

    @property
    def csv(self) -> str:
        return f"{self.name}_csv"

    @property
    def csv_out(self) -> str:
        return f"{self.name}_csv-out"

    @property
    def blurb(self) -> str:
        return f"{self.name}_blurb"
    
    @property
    def phospho_drop(self) -> str:
        return f"{self.name}_phospho-drop"
    
    @property
    def phospho_select_all(self) -> str:
        return f"{self.name}_phospho-select-all"

    @property
    def slider1(self) -> str:
        return f"{self.name}_slider1"

    @property
    def slider2(self) -> str:
        return f"{self.name}_slider2"


@dataclass
class Params:
    """Parameters for box plotting"""

    name: str
    ids: IDs = field(init=False)
    X: str | None = None
    x_axis_title: str | None = None
    y_axis_title: str | None = None
    Y: str | None = None
    COLOUR: str | None = None
    class_name: str = "control-tab"
    legend_title: str | None = None

    def __post_init__(self):
        self.ids = IDs(name=self.name)


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


def draw_box_chart(data: Data, params: Params) -> html.Div:
    """Draws a box and wisker of the CPM data for each set of replicates for eact
    test and overlays the respective FDR value"""

    fig = px.box(
        data.plot_df,
        x=params.X,
        y=data.gene,
        points="all",
        width=666,
        height=666,
        color=params.COLOUR,  # TODO
        log_y=data.stats["log"],
        facet_row_spacing=0.75,
    )  # title=f"Boxplot for CPMs",labels={y_gene: "CPM"},
    if not params.name == "Phosphoproteomics":
        fig.update_xaxes(categoryorder="array", categoryarray=data.ordered_test_names)
    fig = make_brackets(fig, data)
    # Set layout TODO this code is duplicated in func
    fig.update_layout(
        modebar_orientation="v",
        xaxis_title=params.x_axis_title,
        yaxis_title=params.y_axis_title,
        legend_title=params.legend_title,
        legend=dict(
            orientation="h",
            entrywidth=333,
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
        ),
    )

    return html.Div(dcc.Graph(figure=fig), id=params.ids.plot)


def get_sub_poses(number_plots_per_subplot: int, i: int) -> list[float]:
    """Get the positions of the subplots for each test"""

    bargap = 0.2
    gap_width = 1 - bargap
    bar_width = gap_width / number_plots_per_subplot

    gpEdges = i - gap_width / 2

    return [
        gpEdges + bar_width / 2 + bar_width * n for n in range(number_plots_per_subplot)
    ]


def make_brackets(fig: go.Figure, data: Data) -> go.Figure:
    """Adds notations giving the significance level between two box plot data"""
    if data.point_of_reference in data.test_names:
        y_range = get_y_range(len(data.test_names))
        point_of_reference = get_point_of_reference(data)
        if isinstance(data.df_FDR, pd.DataFrame):
            if len(data.df_FDR.columns) == 2 and "test" in data.df_FDR.columns:
                # issue #27 - only bracket if 1 pos... (from James)
                for i, phos_pos in enumerate(data.df_FDR.columns):
                    if phos_pos != "test":
                        sub_poses = get_sub_poses(len(data.test_names), i)
                        add_bracket_per_test(
                            fig,
                            data,
                            point_of_reference,
                            y_range,
                            add_annotation=False,
                            prot_pos=phos_pos,
                            sub_poses=sub_poses,
                        )
        else:
            if len(data.df_FDR) + 1 != len(data.ordered_test_names):
                # issue 13
                data.stats["brackets"] = False
            add_bracket_per_test(fig, data, point_of_reference, y_range)

    return fig


def add_bracket_per_test(
    fig: go.Figure,
    data: Data,
    point_of_reference: int,
    y_range: np.ndarray,
    add_annotation: bool = True,
    prot_pos: str | None = None,
    sub_poses: str | None = None,
) -> go.Figure:  # TODO fix this mess
    """Adds notations giving the significance level between two box plot data"""
    for i, test in enumerate(data.ordered_test_names):
        if test == data.point_of_reference:
            assert i == point_of_reference
            fig = add_annotations_per_test(fig, data, "", test)
            continue
        FDR = data.get_FDR(test, prot_pos) if prot_pos else data.get_FDR(test)

        if data.stats["brackets"]:
            POR = sub_poses[point_of_reference] if sub_poses else point_of_reference
            x = sub_poses[i] if sub_poses else i
            bracket = Bracket(x_POR=POR, x_other=x, y_range=y_range, y_pos=i, FDR=FDR)
            fig = add_FDR_brackets(fig, bracket)
            fig.update_layout(margin=dict(t=i * 15))

        if add_annotation:  # need this or cols with all nan excluded
            if data.stats["brackets"]:
                FDR = f"{FDR:.1e}"
            else:
                FDR = ""
            fig = add_annotations_per_test(fig, data, FDR, test)

    return fig


def add_annotations_per_test(
    fig: go.Figure, data: Data, FDR: str, test: str
) -> go.Figure:
    """Adds notations giving the significance level between two box plot data
    or add blank annotations if all data nan (required to keep it in plot)"""

    fig.add_annotation(
        x=test,
        y=data.get_median_CPMs(test),
        text=FDR,
        yshift=10,
        showarrow=False,
        font=dict(color="black", size=10),
    )

    return fig


def get_point_of_reference(data: Data) -> int:
    """Get the index of the point of reference in the list of tests"""
    point_of_reference_index = [
        i
        for i, test in enumerate(data.ordered_test_names)
        if test == data.point_of_reference
    ]
    assert len(point_of_reference_index) == 1
    return point_of_reference_index[0]


@dataclass(slots=True, frozen=True)
class Bracket:
    """Bracket for FDR notation

    Parameters:
    ----------
    TODO
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

    x_POR: float
    x_other: float
    y_range: np.ndarray
    y_pos: int
    FDR: float
    text_height: float = 1.01
    color: str = "black"

    @property
    def stars(self):
        if self.FDR == 0.0:
            return "ns"
        elif self.FDR >= 0.05:
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


def add_FDR_brackets(fig: go.Figure, bracket: Bracket) -> go.Figure:
    """Adds notations giving the significance level between two box plot
    data (t-test two-sided test)

    Parameters:
    ----------
    fig: figure
        Plotly boxplot figure.
    """
    for line in [
        [bracket.x_POR, bracket.y0, bracket.x_POR, bracket.y1],
        [bracket.x_POR, bracket.y1, bracket.x_other, bracket.y1],
        [
            bracket.x_other,
            bracket.y0,
            bracket.x_other,
            bracket.y1,
        ],
    ]:
        x0, y0, x1, y1 = line
        fig.add_shape(
            type="line",
            xref="x",
            yref="y" + " domain",
            x0=x0,
            y0=y0,
            x1=x1,
            y1=y1,
            line=dict(
                color=bracket.color,
                width=2,
            ),
        )

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


def rubbish(name: str) -> bool:
    """checks if file name starts with . or ~"""
    return name.startswith((".", "~"))
