import numpy as np

import pandas as pd
from dash import dcc, html
import plotly.express as px
from src.read_files import RNASeqData, Data
from .components import ids


def make_list_of_dicts(values: list[str]) -> list[dict[str, str]]:
    """Convert a list of strs into a list where those strings are values in dicts
    against keys label and value, for use in callbacks"""
    return [{"label": val, "value": val} for val in values]


def get_y_range(number_of_comparisons, interline=0.03):
    # Specify in what y_range to plot for each pair of columns
    y_range = np.zeros([number_of_comparisons, 2])
    for i in range(number_of_comparisons):
        y_range[i] = [0.9 + i * interline, 0.91 + i * interline]
    return y_range


def add_FDR_brackets(
    fig, FDR, i, column_pair, y_range, text_height=1.01, color="black"
):
    """Adds notations giving the significance level between two box plot data (t-test two-sided comparison)

    Parameters:
    ----------
    fig: figure
        Plotly boxplot figure.
    FDR: float
        False Discovery Rate (FDR) value for the comparison.
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


def draw_box_chart(
    gene: str, data: Data, div_id: str = ids.BOX_CHART, x: str = "comparison", colour =None, log= False
) -> html.Div:
    """Draws a box and wisker of the CPM data for each set of replicates for eact
    comparison and overlays the respective FDR value"""
    # rna_seq_data: RNASeqData = data[BULK][dataset_choice]

    # fig = px.box(df3, y="abun", x="prot_loc", color="sample", log_y=True, points="all",
    #       hover_data=df3.columns)
    try:
        fig = px.box(
        data.df,
        x=x,
        y=gene,
        points="all",
        width=999,
        height=666,
        color=colour,
        log_y=log,
        title=f"Boxplot for {gene} CPMs",
        labels={"comparison": "Comparison type", gene: "CPM"},
        facet_row_spacing=0.75,
        hover_data=data.df.columns,
    )
    except Exception as e:
        print(data.df.head(), gene, e)
    

    # unique_comparisons = rna_seq_data.raw_df.comparison.unique()
    # y_range = get_y_range(len(unique_comparisons))
    # cytokine_storm = "None"
    # cytokine_storm_FDR = 0.0  # fix this shit...
    # for i, comp in enumerate(unique_comparisons):
    #     if comp == rna_seq_data.point_of_reference:
    #         cytokine_storm = i
    #         DEG_df: pd.DataFrame | None = rna_seq_data.processed_dfs.get(comp)
    #         if DEG_df is not None:
    #             cytokine_storm_FDR = float(DEG_df.query("gene_id == @gene").FDR.iloc[0])
    #         break
    # for i, comp in enumerate(unique_comparisons):
    #     if cytokine_storm == "None":
    #         break  # if removed can't plot
    #     if i == cytokine_storm:
    #         continue
    #     DEG_df: pd.DataFrame | None = rna_seq_data.processed_dfs.get(comp)
    #     if DEG_df is not None:
    #         FDR = float(DEG_df.query("gene_id == @gene").FDR.iloc[0])
    #     else:
    #         FDR = cytokine_storm_FDR
    #         if "CTRL" not in comp.upper():
    #             print(f"log.warn: CTRL not in comp : {comp}")
    #     y = rna_seq_data.raw_df.query("comparison == @comp")[gene].median()
    #     fig.add_annotation(
    #         x=comp,
    #         y=y,
    #         text=f"{FDR:.1e}",
    #         yshift=10,
    #         showarrow=False,
    #     )
    #     fig = add_FDR_brackets(fig, FDR, i, [i, cytokine_storm], y_range)

    return html.Div(dcc.Graph(figure=fig), id=div_id)
