# from dash import Dash, dcc, html
from src.parse_data.read_files import Data
from typing import Any
from src.helpers import draw_box_chart, Params, make_list_of_dicts
import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, callback, callback_context
import dash_daq as daq
from multiprocessing import Pool
from icecream import ic


def download(app: Dash, datasets: dict[str, Data], params: Params) -> Any:
    '''Download the data as a csv file'''
    @callback(
        Output(params.ids.csv_out, "data"),
        Input(params.ids.csv, "n_clicks"),
        Input(params.ids.data_drop, "value"),
        Input(params.ids.gene_drop, "value"),
        Input(params.ids.tests_drop, "value"),
        Input(params.ids.phospho_drop, "value"),
        prevent_initial_call=True,
    )
    def func(n_clicks, dataset: str, gene: str, tests: list[str], phos):
        if "csv" in callback_context.triggered_id:
            selected_data = datasets[dataset]
            if callback_context.triggered_id == "Phosphoproteomics_csv":
                filtered: Data = selected_data.filter(phos, tests, {})
            else:
                filtered: Data = selected_data.filter(gene, tests, {})
            return dcc.send_data_frame(filtered.plot_df.to_csv, "mydf.csv")


def tog_stat(app, params: Params) -> Any:
    '''Toggle the stats on and off'''
    @callback(
        Output(params.ids.stats_out, "value"),
        Input(params.ids.stats, "value"),
    )
    def update_output(value):
        return value


def tog_log(app, params: Params) -> Any:
    '''Toggle the log on and off'''
    @callback(
        Output(params.ids.log_out, "value"),
        Input(params.ids.log, "value"),
    )
    def update_output(value):
        return value


def gene_dropdown(app, datasets: dict[str, Data], params: Params) -> list[dict[str, str]]:
    '''Populate the gene dropdown with options from the given dataset'''
    @app.callback(Output(params.ids.gene_drop, "options"), Input(params.ids.data_drop, "value"))
    def set_gene_options(dataset: str) -> list[dict[str, str]]:
        """Populates the gene selection dropdown with options from teh given dataset"""
        return make_list_of_dicts(list(datasets[dataset].df.columns))


def gene_dropdown_default(app, params: Params) -> str:
    '''Select first gene as default value'''
    @app.callback(Output(params.ids.gene_drop, "value"), Input(params.ids.gene_drop, "options"))
    def select_gene_value(gene_options: list[dict[str, str]]) -> str:
        """Select first gene as default value"""
        return gene_options[0]["value"]


def box(app, datasets: dict[str, Data], params: Params) -> html.Div:
    '''Draw a box and wisker of the CPM data for each set of replicates for each'''
    @app.callback(
        Output(params.ids.plot, "children"),
        Input(params.ids.data_drop, "value"),
        Input(params.ids.gene_drop, "value"),
        Input(params.ids.tests_drop, "value"),
        Input(params.ids.stats_out, "value"),
        Input(params.ids.log_out, "value"),
    )
    def update_box_chart(
        dataset: str, gene: str, tests: list[str], stats: bool, log: bool
    ) -> html.Div:
        """Re draws a box and wisker of the CPM data for each set of replicates for eact
        test and overlays the respective FDR value"""
        selected_data = datasets[dataset]
        stats_d = {"brackets": stats, "log": log}
        filtered: Data = selected_data.filter(gene, tests, stats_d)
        return draw_box_chart(filtered, params)


def test_dropdown(app, datasets: dict[str, Data], params: Params) -> list[dict[str, str]]:
    '''Populate the test selection dropdown with options from teh given dataset'''
    @app.callback(
        Output(params.ids.tests_drop, "options"),
        Input(params.ids.data_drop, "value"),
    )
    def set_comparison_options(dataset: str) -> list[dict[str, str]]:
        """Populates the test selection dropdown with options from teh given dataset"""
        return make_list_of_dicts(list(datasets[dataset].test_names))


def test_dropdown_select_all(app, params: Params) -> list[dict[str, str]]:
    '''Default to all available comparisons'''
    @app.callback(
        Output(params.ids.tests_drop, "value"),
        Input(params.ids.tests_drop, "options"),
        Input(params.ids.select_all, "n_clicks"),
    )
    def select_comparison_values(
        available_comparisons: list[dict[str, str]], _: int
    ) -> list[dict[str, str]]:
        """Default to all available comparisons"""
        return [comp["value"] for comp in available_comparisons]


def phospho_site_dropdown(app, datasets: dict[str, Data], params: Params) -> list[dict[str, str]]:
    '''Populate the phospho site dropdown with options from teh given dataset'''''
    @app.callback(
        Output(params.ids.phospho_drop, "options"),
        Input(params.ids.data_drop, "value"),
        Input(params.ids.gene_drop, "value"),
        Input(params.ids.tests_drop, "value"),
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


def phospho_site_dropdown_select_all(app, params: Params) -> list[dict[str, str]]:
    '''Select all phospho sites'''
    @app.callback(
        Output(params.ids.phospho_drop, "value"),
        Input(params.ids.phospho_drop, "options"),
        Input(params.ids.phospho_select_all, "n_clicks"),
    )
    def select_comparison_values(
        available_comparisons: list[dict[str, str]], _: int
    ) -> list[dict[str, str]]:
        """Default to all available comparisons"""
        return [comp["value"] for comp in available_comparisons]


def phospho_site_box(app, datasets: dict[str, Data], params: Params) -> html.Div:
    '''Draw multiple box and wisker of the CPM data for each set of replicates for each'''
    @app.callback(
        Output(params.ids.plot, "children"),
        Input(params.ids.data_drop, "value"),
        Input(params.ids.gene_drop, "value"),
        Input(params.ids.tests_drop, "value"),
        Input(params.ids.phospho_drop, "value"),
        Input(params.ids.stats_out, "value"),
        Input(params.ids.log_out, "value"),
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
        return draw_box_chart(filtered, params)


def get_defaults(datasets: dict[str, Data], params) -> tuple[str, Data, list[str]]:
    '''Get the default dataset, gene and test'''
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


def blurby(app, initialised_datasets, params: Params) -> Any:
    '''Write the blurb'''
    @app.callback(Output(params.ids.blurb, "children"), Input(params.ids.data_drop, "value"))
    def write_blurb(dataset: str) -> html.Div:
        ic(dataset)
        selected_data = initialised_datasets[dataset]
        children = []
        for line in selected_data.blurb:
            children += [html.P("\t" + line)]
        return html.Div(children=children, id=params.ids.blurb)


def dropdowns(datasets: dict[str, Data], params: Params) -> list[Any]:
    '''Create the dropdowns'''
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
                        id=params.ids.data_drop,
                        options=dataset_names,
                        value=dataset_names[0],
                        multi=False,
                    ),
                    width=4,
                ),
                dbc.Col(
                    dcc.Dropdown(
                        id=params.ids.gene_drop,
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
                        id=params.ids.tests_drop,
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
                        id=params.ids.select_all,
                        n_clicks=0,
                    ),
                    width=3,
                ),
                dbc.Col(
                    daq.ToggleSwitch(id=params.ids.stats_out, value=True, label="Stats"),
                    width=1,
                ),
                dbc.Col(
                    daq.ToggleSwitch(id=params.ids.log_out, value=True, label="Log10"),
                    width=1,
                ),
                dbc.Col(
                    [
                        html.Button("Download CSV", id=params.ids.csv, n_clicks=0),
                        dcc.Download(id=params.ids.csv_out),
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
                            id=params.ids.phospho_drop,
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
                            id=params.ids.phospho_select_all,
                            n_clicks=0,
                        ),
                        width=5,
                    )
                ]
            ),
            html.Div(draw_box_chart(dataset, params)),
        ]
    else:
        children += [html.Div(draw_box_chart(dataset, params))]
    children += [html.Div(dataset.blurb, id=params.ids.blurb)]
    return children


def render(app: Dash, datasets: dict[str, Data], params: Params) -> html.Div:
    '''Render the tab'''
    initialised_datasets = {}

    with Pool(processes=3) as pool:
        for k, result in pool.imap_unordered(call_funk, datasets.items()):
            initialised_datasets[k] = result

    gene_dropdown(app, initialised_datasets, params)
    gene_dropdown_default(app, params)
    test_dropdown(app, initialised_datasets, params)
    test_dropdown_select_all(app, params)
    tog_stat(app, params)
    tog_log(app, params)
    download(app, initialised_datasets, params)
    if params.name == "Phosphoproteomics":
        phospho_site_dropdown(app, initialised_datasets, params)
        phospho_site_dropdown_select_all(app, params)
        phospho_site_box(app, initialised_datasets, params)
    else:
        box(app, initialised_datasets, params)
    blurby(app, initialised_datasets, params)

    return html.Div(children=dropdowns(initialised_datasets, params))


def call_funk(tup) -> tuple[str, Data]:
    '''Call the function and return the result'''
    k, v = tup
    func, path = v
    return k, func(path)
