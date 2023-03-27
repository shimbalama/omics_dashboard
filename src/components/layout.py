from dash import Dash, html, dcc

from . import rna_raw_tab, box_chart, rna_processed_tab  # bar_chart, nation_dropdown,


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
                        dcc.Tab(
                            label="RNAseq (bulk)",
                            value="tab1",
                            children=[
                                dcc.Tabs(
                                    id="subtabs",
                                    value="subtab1",
                                    children=[
                                        dcc.Tab(
                                            label="raw",
                                            value="raw",
                                            children=html.Div(
                                                className="control-tab",
                                                children=[
                                                    html.Div(
                                                        className="app-controls-block",
                                                        children=[
                                                            html.H1(app.title),
                                                            html.Hr(),
                                                            html.Div(
                                                                className="dropdown-container",
                                                                children=[
                                                                    rna_raw_tab.render(
                                                                        app, data
                                                                    )
                                                                ],
                                                            ),
                                                        ],
                                                    )
                                                ],
                                            ),
                                        ),
                                        dcc.Tab(
                                            label="processed",
                                            value="processed",
                                            children=html.Div(
                                                className="control-tab",
                                                children=[
                                                    html.Div(
                                                        className="app-controls-block",
                                                        children=[
                                                            html.H1(app.title),
                                                            html.Hr(),
                                                            html.Div(
                                                                className="dropdown-container",
                                                                children=[
                                                                    rna_processed_tab.render(
                                                                        app, data
                                                                    )
                                                                ],
                                                            ),
                                                        ],
                                                    )
                                                ],
                                            ),
                                        ),
                                    ],
                                )
                            ],
                        ),
                        dcc.Tab(
                            label="RNA_single_cell",
                            value="RNA_single_cell",
                            children=html.Div(
                                className="control-tab",
                                children=[
                                    html.Div(
                                        className="app-controls-block",
                                        children=[
                                            html.Div(
                                                className="app-controls-name",
                                                children="Dataset: ",
                                            )
                                        ],
                                    )
                                ],
                            ),
                        ),
                        dcc.Tab(
                            label="Protein",
                            value="Protein",
                            children=html.Div(
                                className="control-tab",
                                children=[
                                    html.Div(
                                        className="app-controls-block",
                                        children=[
                                            html.Div(
                                                className="app-controls-name",
                                                children="Dataset: ",
                                            )
                                        ],
                                    )
                                ],
                            ),
                        ),
                        dcc.Tab(
                            label="Phosphoproteomics",
                            value="Phosphoproteomics",
                            children=html.Div(
                                className="control-tab",
                                children=[
                                    html.Div(
                                        className="app-controls-block",
                                        children=[
                                            html.Div(
                                                className="app-controls-name",
                                                children="Dataset: ",
                                            )
                                        ],
                                    )
                                ],
                            ),
                        ),
                        # dcc.Tab(
                        #     label="data",
                        #     value="what-is",
                        #     children=html.Div(
                        #         className="control-tab",
                        #         children=[
                        #             html.H1(app.title),
                        #             html.Hr(),
                        #             html.Div(
                        #                 className="dropdown-container",
                        #                 children=[nation_dropdown.render(app)],
                        #             ),
                        #             bar_chart.render(app),
                        #         ],
                        #     ),
                        # ),
                    ],
                )
            ],
        ),
    )


# def create_layout(app: Dash) -> html.Div:
#     return html.Div(
#         className="app-div",
#         children=[
#             html.H1(app.title),
#             html.Hr(),
#             html.Div(
#                 className="dropdown-container",
#                 children=[
#                     nation_dropdown.render(app),
#                 ],
#             ),
#             bar_chart.render(app),
#         ],
#     )
