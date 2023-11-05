from dash import Dash, html, dcc
from ..tabs import box_tabs, rna_processed_tab, text_tab1, text_tab2, function_tab
from src.parse.read_files import Data
from src.helpers import Params
from typing import Callable

__version__ = '1.0.01'

def tab_layout(name: str, subtabs_id: str, *args) -> dcc.Tab:
    """Defines tab/subtab structure"""
    return dcc.Tab(
        label=name,
        value=name,
        children=[dcc.Tabs(id=subtabs_id, value=subtabs_id, children=args)],
    )


def sub_tab_layout(
    app: Dash,
    data: dict[str, Data],
    params: Params,
    render: Callable[[Dash, dict[str, Data], Params], html.Div],
) -> dcc.Tab:
    """Defines base unit of UI"""

    return dcc.Tab(
        label=params.name,
        value=params.name,
        children=html.Div(
            className=params.class_name,
            children=[
                html.Div(
                    className=params.class_name,
                    children=[
                        html.H1(app.title),
                        html.Hr(),
                        html.Div(
                            className=params.class_name,
                            children=[render(app, data, params)],
                        ),
                    ],
                )
            ],
        ),
    )


def create_layout(app: Dash, data: dict[str, dict[str, Data]]) -> html.Div:
    """Create the layout for Dash application.

    Args:
        app (Dash): The Dash application object.
        data (Dict[str, Dict[str, 'Data']]): A dictionary containing the data needed for the layout.

    Returns:
        html.Div: A Dash html.Div object that represents the layout of the application.
    """

    return html.Div(
        className="div-app",
        id="div-app",
        children=[
            html.Div(
                html.Div(
                    id="vp-control-tabs",
                    className="control-tabs",
                    children=[
                        dcc.Tabs(
                            id="vp-tabs",
                            value="what-is",
                            children=[
                                tab_layout(
                                    "User guide",
                                    "subtabs_id0",
                                    sub_tab_layout(
                                        app,
                                        data,
                                        Params(name="Introduction"),
                                        text_tab1.render,
                                    ),
                                    sub_tab_layout(
                                        app,
                                        data,
                                        Params(name="Input_files"),
                                        text_tab2.render,
                                    ),
                                ),
                                tab_layout(
                                    "RNA",
                                    "subtabs_id1",
                                    sub_tab_layout(
                                        app,
                                        data["rna_bulk"],
                                        Params(
                                            name="CPM",
                                            X="test",
                                            y_axis_title="CPM",
                                            x_axis_title="Condition",
                                        ),
                                        box_tabs.render,
                                    ),
                                    sub_tab_layout(
                                        app,
                                        data["rna_bulk"],
                                        Params(name="Volcano", X="test"),
                                        rna_processed_tab.render,
                                    ),
                                ),
                                tab_layout(
                                    "Protein",
                                    "subtabs_id4",
                                    sub_tab_layout(
                                        app,
                                        data["proteomics"],
                                        Params(
                                            name="Proteins",
                                            X="test",
                                            y_axis_title="Protein abundance",
                                            x_axis_title="Condition",
                                        ),
                                        box_tabs.render,
                                    ),
                                    sub_tab_layout(
                                        app,
                                        data["phosphoproteomics"],
                                        Params(
                                            name="Phosphoproteomics",
                                            X="gene",
                                            COLOUR="test",
                                            Y="abun",
                                            x_axis_title="Phosphopeptide",
                                            y_axis_title="Phosphopeptide abundance",
                                            legend_title="Condition",
                                        ),
                                        box_tabs.render,
                                    ),
                                ),
                                tab_layout(
                                    "Other",
                                    "subtabs_id5",
                                    sub_tab_layout(
                                        app,
                                        data["function"]["test"],
                                        Params(name="Function", X="test"),
                                        function_tab.render,
                                    ),
                                ),
                            ],
                        )
                    ],
                ),
            ),
            html.Div(f"Version: {__version__}"),
        ],
    )
