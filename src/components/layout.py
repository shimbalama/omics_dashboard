from dash import Dash, html, dcc

from . import bar_chart, nation_dropdown




def create_layout(app: Dash) -> html.Div:
    return html.Div(
        html.Div(
            id='vp-control-tabs', className='control-tabs', children=[
                dcc.Tabs(id='vp-tabs', value='what-is', children=[
                    dcc.Tab(label='RNA_raw',value='RNA_raw', children=html.Div(className='control-tab', 
                            children=[
                                html.Div(className='app-controls-block', children=[
                                    html.Div(className='app-controls-name',children='Dataset: ')
                                    ]
                                )
                            ]
                        )
                    ),
                    dcc.Tab(label='RNA_processed',value='RNA_processed', children=html.Div(className='control-tab', 
                            children=[
                                html.Div(className='app-controls-block', children=[
                                    html.Div(className='app-controls-name',children='Dataset: ')
                                    ]
                                )
                            ]
                        )
                    ),
                    dcc.Tab(label='RNA_single_cell',value='RNA_single_cell', children=html.Div(className='control-tab', 
                            children=[
                                html.Div(className='app-controls-block', children=[
                                    html.Div(className='app-controls-name',children='Dataset: ')
                                    ]
                                )
                            ]
                        )
                    ),
                    dcc.Tab(label='Protein',value='Protein', children=html.Div(className='control-tab', 
                            children=[
                                html.Div(className='app-controls-block', children=[
                                    html.Div(className='app-controls-name',children='Dataset: ')
                                    ]
                                )
                            ]
                        )
                    ),
                    dcc.Tab(label='Phosphoproteomics',value='Phosphoproteomics', children=html.Div(className='control-tab', 
                            children=[
                                html.Div(className='app-controls-block', children=[
                                    html.Div(className='app-controls-name',children='Dataset: ')
                                    ]
                                )
                            ]
                        )
                    ),
                    dcc.Tab(label='data', value='what-is', children=html.Div(className='control-tab', children=[
                        html.H1(app.title),
                        html.Hr(),
                        html.Div(
                            className="dropdown-container",
                            children=[
                                nation_dropdown.render(app)
                            ]
                        ),
                        bar_chart.render(app),
                                ]
                            )
                        )
                    ]
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