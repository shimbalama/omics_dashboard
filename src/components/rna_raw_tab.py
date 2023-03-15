from dash import Dash, dcc, html
from dash.dependencies import Input, Output
import pandas as pd
from . import ids


def render(app: Dash, dfs: dict[str, pd.DataFrame]) -> html.Div:
    # see https://dash.plotly.com/basic-callbacks#dash-app-with-chained-callbacks

    @app.callback(Output(ids.GENE_DROPDOWN, "options"), Input(ids.DATASET_DROPDOWN, "value"))
    def set_comparison_options(selected_comparison):
        a =  [{"label": gene, "value": gene} for gene in dfs[selected_comparison].columns]
        print(a[:3], len(a))
        return a
    
    @app.callback(Output(ids.COMPARISON_DROPDOWN, "options"), Input(ids.DATASET_DROPDOWN, "value"))
    def set_comparison_options(selected_comparison):
        seen = set([])
        b= []
        for comp in dfs[selected_comparison].comparison:
            if comp not in seen:
                b.append({"label": comp, "value": comp})
            seen.add(comp)
        print( b)
        return b
    
    @app.callback(Output(ids.GENE_DROPDOWN, "value"), Input(ids.GENE_DROPDOWN, "options"))
    def set_gene_value(available_options):
        return available_options[0]["value"]

    @app.callback(Output(ids.COMPARISON_DROPDOWN, "value"), Input(ids.COMPARISON_DROPDOWN, "options"))
    def set_cities_value(available_options):
        return available_options
    # @app.callback(
    #             Output('vp-dataset-div', 'title'),
    #     [
    #         Input(ids.DATASET_DROPDOWN, 'value')
    #     ]
    # )
    # def select_df(dataset_name: str) -> pd.DataFrame:
    #     return dfs[dataset_name]

    # @app.callback(
    #             Output('vp-dataset-div', 'title'),
    #     [
    #         Input(ids.DATASET_DROPDOWN, 'value')
    #     ]
    # )
    # def select_gene() -> pd.DataFrame:
    #     return dfs[dataset_name]

    @app.callback(
        Output(ids.COMPARISON_DROPDOWN, "value"),
        Input(ids.COMPARISON_DROPDOWN, "value"),
        Input(ids.SELECT_ALL_COMPARISONS_BUTTON, "n_clicks")
    )
    def select_all_comparisons() -> list[str]:
        return list(set(df.comparison))

    return html.Div(
        children=[
            html.H6("Dataset"),
            dcc.Dropdown(
                id=ids.DATASET_DROPDOWN,
                options=list(dfs.keys()),
                value="CSexp2_shRNA",
                multi=False,
            ),
            html.H6("Gene"),
            dcc.Dropdown(
                id=ids.GENE_DROPDOWN,
            ),
            html.H6("Comparison"),
            dcc.Dropdown(
                id=ids.COMPARISON_DROPDOWN,
                multi=True,
            ),
            html.Button(
                className="dropdown-button",
                children=["Select All"],
                id=ids.SELECT_ALL_COMPARISONS_BUTTON,
                n_clicks=0,
            ),
        ]
    )
