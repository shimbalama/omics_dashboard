from dash import Dash, html, dcc

from . import rna_raw_tab, rna_processed_tab, text_tab1, text_tab2


def tab_layout(name, subtabs_id, *args):
    return dcc.Tab(
        label=name,
        value=name,
        children=[dcc.Tabs(id=subtabs_id, value=subtabs_id, children=args)],
    )


def sub_tab_layout(app, data, name, rendering, class_name="control-tab"):
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
                            children=[rendering.render(app, data)],
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
                            sub_tab_layout(app, data, "Input files", text_tab2),
                        ),
                        tab_layout(
                            "RNAseq (bulk)",
                            "subtabs_id1",
                            sub_tab_layout(app, data, "CPM and FDR", rna_raw_tab),
                            sub_tab_layout(
                                app, data, "Volcano of DEG", rna_processed_tab
                            ),
                        ),
                        tab_layout("RNA_single_cell", "subtabs_id2"),
                        tab_layout("Protein", "subtabs_id3"),
                        tab_layout("Phosphoproteomics", "subtabs_id4"),
                        tab_layout("Function", "subtabs_id5"),
                    ],
                )
            ],
        ),
    )
