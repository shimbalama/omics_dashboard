from dash import Dash, dcc, html
from src.read_files import Data
from dash.dependencies import Input, Output
from typing import Any
from src.helpers import draw_box_chart, Params, IDs, make_list_of_dicts

def gene_dropdown(app, ids: IDs, data_set: dict[str, Data]):
    @app.callback(Output(ids.gene_drop, "options"), Input(ids.data_drop, "value"))
    def set_gene_options(experiment: str) -> list[dict[str, str]]:
        """Populates the gene selection dropdown with options from teh given dataset"""
        return make_list_of_dicts(list(data_set[experiment].df.columns))


def gene_dropdown_default(app, ids: IDs):
    @app.callback(Output(ids.gene_drop, "value"), Input(ids.gene_drop, "options"))
    def select_gene_value(gene_options: list[dict[str, str]]) -> str:
        """Select first gene as default value"""
        return gene_options[0]["value"]


def box(app, ids: IDs, data_set: Data, params: type):
    @app.callback(
        Output(ids.plot, "children"),
        Input(ids.data_drop, "value"),
        Input(ids.gene_drop, "value"),
        Input(ids.tests_drop, "value"),
    )
    def update_box_chart(experiment: str, gene: str, tests: list[str]) -> html.Div:
        """Re draws a box and wisker of the CPM data for each set of replicates for eact
        test and overlays the respective FDR value"""
        selected_data = data_set[experiment]
        filtered: Data = selected_data.filter(gene, tests)
        y_param = params.Y if params.Y else gene
        return draw_box_chart(filtered, y_param, params, ids.plot)


def test_dropdown(app, ids: IDs, data_set: Data):
    @app.callback(
        Output(ids.tests_drop, "options"),
        Input(ids.data_drop, "value"),
    )
    def set_comparison_options(experiment: str) -> list[dict[str, str]]:
        """Populates the test selection dropdown with options from teh given dataset"""
        return make_list_of_dicts(list(data_set[experiment].test_names))


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


def get_defaults(data: Data) -> tuple[str, Data]:
    datasets = list(data.keys())
    first_dataset = data[datasets[0]]
    first_gene = first_dataset.df.columns[0]
    dataset = first_dataset.filter(first_gene, list(first_dataset.test_names))

    return first_gene, dataset, datasets


def dropdowns(data: Data, params: Params, ids: IDs) -> list[Any]:
    first_gene, dataset, datasets = get_defaults(data)
    y_param = params.Y if params.Y else first_gene
    children = [
        html.H6("Dataset"),
        dcc.Dropdown(
            id=ids.data_drop,
            options=datasets,
            value=datasets[0],
            multi=False,
        ),
        html.H6("Gene"),
        dcc.Dropdown(
            id=ids.gene_drop,
        ),
        html.H6("test"),
        dcc.Dropdown(
            id=ids.tests_drop,
            multi=True,
        ),
        html.Button(
            className="dropdown-button",
            children=["Select All"],
            id=ids.select_all,
            n_clicks=0,
        ),
        html.Div(draw_box_chart(dataset, y_param, params, ids.plot)),
    ]

    return children


def render(app: Dash, datasets: dict[str, Data], ids: IDs, params: Params) -> html.Div:
    gene_dropdown(app, ids, datasets)
    gene_dropdown_default(app, ids)
    test_dropdown(app, ids, datasets)
    test_dropdown_select_all(app, ids)
    box(app, ids, datasets, params)
    return html.Div(children=dropdowns(datasets, params, ids))
