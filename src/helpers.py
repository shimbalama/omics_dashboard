import numpy as np

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
    fig,
    FDR,
    i,
    column_pair,
    y_range,
    text_height=1.01,
    color="black"
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
