# from dash import Dash, dcc, html
from src.read_files import Data

# from dash.dependencies import Input, Output
from typing import Any
from src.helpers import draw_box_chart, Params, IDs, make_list_of_dicts
import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, callback, callback_context
import dash_daq as daq


def download(app, ids: IDs, datasets: dict[str, Data]):
    @callback(
        Output(ids.csv_out, "data"),
        Input(ids.csv, "n_clicks"),
        Input(ids.data_drop, "value"),
        Input(ids.gene_drop, "value"),
        Input(ids.tests_drop, "value"),
        Input("phos_site_drop", "value"),
        prevent_initial_call=True,
    )
    def func(n_clicks, dataset: str, gene: str, tests: list[str], phos):
        if 'csv' in callback_context.triggered_id:
            selected_data = datasets[dataset]
            if callback_context.triggered_id == 'Phosphoproteomics_csv':
                filtered: Data = selected_data.filter(phos, tests, {})
            else:
                filtered: Data = selected_data.filter(gene, tests, {})
            return dcc.send_data_frame(filtered.plot_df.to_csv, "mydf.csv")


def tog_stat(app, ids: IDs):
    @callback(
        Output(ids.stats_out, "value"),
        Input(ids.stats, "value"),
    )
    def update_output(value):
        return value


def tog_log(app, ids: IDs):
    @callback(
        Output(ids.log_out, "value"),
        Input(ids.log, "value"),
    )
    def update_output(value):
        return value


def gene_dropdown(app, ids: IDs, datasets: dict[str, Data]):
    @app.callback(Output(ids.gene_drop, "options"), Input(ids.data_drop, "value"))
    def set_gene_options(dataset: str) -> list[dict[str, str]]:
        """Populates the gene selection dropdown with options from teh given dataset"""
        return make_list_of_dicts(list(datasets[dataset].df.columns))


def gene_dropdown_default(app, ids: IDs):
    @app.callback(Output(ids.gene_drop, "value"), Input(ids.gene_drop, "options"))
    def select_gene_value(gene_options: list[dict[str, str]]) -> str:
        """Select first gene as default value"""
        return gene_options[0]["value"]


def box(app, ids: IDs, datasets: dict[str, Data], params: type):
    @app.callback(
        Output(ids.plot, "children"),
        Input(ids.data_drop, "value"),
        Input(ids.gene_drop, "value"),
        Input(ids.tests_drop, "value"),
        Input(ids.stats_out, "value"),
        Input(ids.log_out, "value"),
    )
    def update_box_chart(
        dataset: str, gene: str, tests: list[str], stats: bool, log: bool
    ) -> html.Div:
        """Re draws a box and wisker of the CPM data for each set of replicates for eact
        test and overlays the respective FDR value"""
        selected_data = datasets[dataset]
        stats_d = {"brackets": stats, "log": log}
        filtered: Data = selected_data.filter(gene, tests, stats_d)
        return draw_box_chart(filtered, params, ids.plot)


def test_dropdown(app, ids: IDs, datasets: dict[str, Data]):
    @app.callback(
        Output(ids.tests_drop, "options"),
        Input(ids.data_drop, "value"),
    )
    def set_comparison_options(dataset: str) -> list[dict[str, str]]:
        """Populates the test selection dropdown with options from teh given dataset"""
        return make_list_of_dicts(list(datasets[dataset].test_names))


def test_dropdown_select_all(app, ids: IDs):
    @app.callback(
        Output(ids.tests_drop, "value"),
        Input(ids.tests_drop, "options"),
        Input(ids.select_all, "n_clicks"),
    )
    def select_comparison_values(
        available_comparisons: list[dict[str, str]], _: int
    ) -> list[dict[str, str]]:
        """Default to all available comparisons"""
        return [comp["value"] for comp in available_comparisons]


def phospho_site_dropdown(app, ids: IDs, datasets: dict[str, Data]):
    @app.callback(
        Output("phos_site_drop", "options"),
        Input(ids.data_drop, "value"),
        Input(ids.gene_drop, "value"),
        Input(ids.tests_drop, "value"),
    )
    def set_comparison_options(
        dataset: str, gene: str, tests: list[str]
    ) -> list[dict[str, str]]:
        """Populates the test selection dropdown with options from teh given dataset"""
        selected_data = datasets[dataset]
        filtered: Data = selected_data.filter(
            [col for col in selected_data.df_FDR.columns if col.startswith(f"{gene}_")],
            tests,
            True,
        )
        return make_list_of_dicts(
            [phos for phos in list(filtered.df_FDR.columns) if phos != "test"]
        )


def phospho_site_dropdown_select_all(app, ids: IDs):
    @app.callback(
        Output("phos_site_drop", "value"),
        Input("phos_site_drop", "options"),
        Input("phos_site_drop_select_all", "n_clicks"),
    )
    def select_comparison_values(
        available_comparisons: list[dict[str, str]], _: int
    ) -> list[dict[str, str]]:
        """Default to all available comparisons"""
        return [comp["value"] for comp in available_comparisons]


def phospho_site_box(app, ids: IDs, datasets: dict[str, Data], params: type):
    @app.callback(
        Output(ids.plot, "children"),
        Input(ids.data_drop, "value"),
        Input(ids.gene_drop, "value"),
        Input(ids.tests_drop, "value"),
        Input("phos_site_drop", "value"),
        Input(ids.stats_out, "value"),
        Input(ids.log_out, "value"),
    )
    def update_box_chart(
        dataset: str,
        gene: str,
        tests: list[str],
        phos_sites: list[str],
        stats: bool,
        log: bool,
    ) -> html.Div:
        """Re draws a box and wisker of the CPM data for each set of replicates for eact
        test and overlays the respective FDR value"""
        stats_d = {"brackets": stats, "log": log}
        selected_data = datasets[dataset]
        filtered: Data = selected_data.filter(phos_sites, tests, stats_d)

        return draw_box_chart(filtered, params, ids.plot)


def get_defaults(datasets: dict[str, Data], params) -> tuple[str, Data, list[str]]:
    dataset_names = list(datasets.keys())
    first_dataset = datasets[dataset_names[0]]
    if params.name == "Phosphoproteomics":
        first_gene = sorted(first_dataset.df.columns)[0]
        dataset = first_dataset.filter(
            [sorted(first_dataset.df_FDR.columns)[0]],
            list(first_dataset.test_names),
            {"brackets": True, "log": True},
        )
    else:
        first_gene = first_dataset.df.columns[0]
        dataset = first_dataset.filter(
            first_gene, list(first_dataset.test_names), {"brackets": True, "log": False}
        )

    return dataset, dataset_names


def dropdowns(datasets: dict[str, Data], params: Params, ids: IDs) -> list[Any]:
    dataset, dataset_names = get_defaults(datasets, params)

    children = [
        dbc.Row(
            [
                dbc.Col(html.H6("Dataset"), width=4),
                dbc.Col(html.H6("Gene"), width=3),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Dropdown(
                        id=ids.data_drop,
                        options=dataset_names,
                        value=dataset_names[0],
                        multi=False,
                    ),
                    width=4,
                ),
                dbc.Col(
                    dcc.Dropdown(
                        id=ids.gene_drop,
                    ),
                    width=3,
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(html.H6("Test"), width=7),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Dropdown(
                        id=ids.tests_drop,
                        multi=True,
                    ),
                    width=7,
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Button(
                        className="dropdown-button",
                        children=["Select All"],
                        id=ids.select_all,
                        n_clicks=0,
                    ),
                    width=3,
                ),
                dbc.Col(
                    daq.ToggleSwitch(id=ids.stats_out, value=True, label="Stats"),
                    width=1,
                ),
                dbc.Col(
                    daq.ToggleSwitch(id=ids.log_out, value=True, label="Log10"),
                    width=1,
                ),
                dbc.Col(
                    [
                        html.Button("Download CSV", id=ids.csv, n_clicks=0),
                        dcc.Download(id=ids.csv_out),
                    ],
                    width=2,
                ),
            ],
        ),
    ]
    if params.name == "Phosphoproteomics":
        children += [
            dbc.Row(
                [
                    dbc.Col(html.H6("Phospho Site"), width=3),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dcc.Dropdown(
                            id="phos_site_drop",
                            multi=True,
                            persistence=True,
                            persistence_type="session",
                        ),
                        width=7,
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        html.Button(
                            className="dropdown-button",
                            children=["Select All"],
                            id="phos_site_drop_select_all",
                            n_clicks=0,
                        ),
                        width=5,
                    )
                ]
            ),
            html.Div(draw_box_chart(dataset, params, ids.plot)),
        ]
    else:
        children += [html.Div(draw_box_chart(dataset, params, ids.plot))]

    return children


def render(app: Dash, datasets: dict[str, Data], ids: IDs, params: Params) -> html.Div:
    initialised_datasets = {}
    for k, v in datasets.items():
        func, path = v
        initialised_datasets[k] = func(path)
    gene_dropdown(app, ids, initialised_datasets)
    gene_dropdown_default(app, ids)
    test_dropdown(app, ids, initialised_datasets)
    test_dropdown_select_all(app, ids)
    tog_stat(app, ids)
    tog_log(app, ids)
    download(app, ids, initialised_datasets)
    if params.name == "Phosphoproteomics":
        phospho_site_dropdown(app, ids, initialised_datasets)
        phospho_site_dropdown_select_all(app, ids)
        phospho_site_box(app, ids, initialised_datasets, params)
    else:
        box(app, ids, initialised_datasets, params)
    return html.Div(children=dropdowns(initialised_datasets, params, ids))
