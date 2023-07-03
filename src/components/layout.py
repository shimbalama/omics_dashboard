from dash import Dash, html, dcc

from ..tabs import (
    box_tabs,
    rna_processed_tab,
    text_tab1,
    text_tab2,
    function_tab
)
from src.read_files import Data

from src.helpers import IDs, Params


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
    name: str,
    rendering,
    class_name: str = "control-tab",
):
    '''Defines base unit of UI'''
    ids = IDs(name)
    params = (
        Params(X="gene", COLOUR="test", LOG=True, Y="abun")
        if name == "Phosphoproteomics"
        else Params(X="test")
    )
    return dcc.Tab(
        label=name,
        value=name,
        children=html.Div(
            className=class_name,
            children=[
                html.Div(
                    className=class_name,
                    children=[
                        html.H1(app.title),
                        html.Hr(),
                        html.Div(
                            className=class_name,
                            children=[rendering.render(app, data, ids, params)],
                        ),
                    ],
                )
            ],
        ),
    )


def create_layout(app: Dash, data: dict[str, dict[str, Data]]) -> html.Div:
    return html.Div(
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
                            sub_tab_layout(app, data, "Introduction", text_tab1),
                            sub_tab_layout(app, data, "Input_files", text_tab2),
                        ),
                        tab_layout(
                            "RNA",
                            "subtabs_id1",
                            sub_tab_layout(app, data["rna_bulk"], "CPM", box_tabs),
                            sub_tab_layout(
                                app, data["rna_bulk"], "Volcano", rna_processed_tab
                            ),
                            sub_tab_layout(app, data, "scRNA", text_tab1),
                        ),
                        tab_layout(
                            "Protein",
                            "subtabs_id4",
                            sub_tab_layout(
                                app, data["proteomics"], "Proteins", box_tabs
                            ),
                            sub_tab_layout(
                                app,
                                data["phosphoproteomics"],
                                "Phosphoproteomics",
                                box_tabs,
                            ),
                        ),
                        tab_layout(
                            "Other",
                            "subtabs_id5",
                            sub_tab_layout(
                                app, data["function"]['test'], "Function", function_tab
                            ),
                        ),
                    ],
                )
            ],
        ),
    )
