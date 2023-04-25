from dash import Dash, html, dcc

from ..tabs import (
    box_tabs,
    rna_processed_tab,
    text_tab1,
    text_tab2,
)

from src.helpers import IDs, Params


def tab_layout(name, subtabs_id, *args):
    return dcc.Tab(
        label=name,
        value=name,
        children=[dcc.Tabs(id=subtabs_id, value=subtabs_id, children=args)],
    )


def sub_tab_layout(app, data, name, rendering, class_name="control-tab"):
    ids = IDs(name)
    params = (
        Params(X="gene", COLOUR="ID", LOG=True, Y="abun")
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


def create_layout(app: Dash, data) -> html.Div:
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
                            "RNAseq (bulk)",
                            "subtabs_id1",
                            sub_tab_layout(app, data["rna_bulk"], "CPM", box_tabs),
                            sub_tab_layout(
                                app, data["rna_bulk"], "Volcano", rna_processed_tab
                            ),
                        ),
                        tab_layout("RNA_single_cell", "subtabs_id2"),
                        tab_layout(
                            "Phospho/proteomics",
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
                        tab_layout("Function", "subtabs_id5"),
                    ],
                )
            ],
        ),
    )
