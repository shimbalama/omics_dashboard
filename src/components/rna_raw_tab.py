from dash import Dash, dcc, html
from dash.dependencies import Input, Output
import pandas as pd
from . import ids
import plotly.express as px


def render(app: Dash, dfs: dict[str, pd.DataFrame]) -> html.Div:
    # see https://dash.plotly.com/basic-callbacks#dash-app-with-chained-callbacks

    # from dataset dropdown get 3 things
    # 1: get df
    # @app.callback(Output(ids.BOX_CHART, "value"), Input(ids.RAW_RNA_DATA_DROP, "value"))
    # def get_df(selected_comparison):
    #     return selected_comparison

    # 2: get genes
    @app.callback(
        Output(ids.GENE_DROPDOWN, "options"), Input(ids.RAW_RNA_DATA_DROP, "value")
    )
    def set_gene_options(selected_comparison):
        a = [
            {"label": gene, "value": gene} for gene in dfs[selected_comparison].columns
        ]
        print(a[:3], len(a))
        return a

    # 3: get comparisons
    @app.callback(
        Output(ids.COMPARISON_DROPDOWN, "options"), Input(ids.RAW_RNA_DATA_DROP, "value")
    )
    def set_comparison_options(selected_comparison):
        seen = set([])
        b = []
        for comp in dfs[selected_comparison].comparison:
            if comp not in seen:
                b.append({"label": comp, "value": comp})
            seen.add(comp)
        print(b)
        return b

    # select gene
    @app.callback(
        Output(ids.GENE_DROPDOWN, "value"), Input(ids.GENE_DROPDOWN, "options")
    )
    def select_gene_value(available_options):
        return available_options[0]["value"]

    # select comparison/s
    @app.callback(
        Output(ids.COMPARISON_DROPDOWN, "value"),
        Input(ids.COMPARISON_DROPDOWN, "options"),
    )
    def select_comparison_values(available_options):
        return available_options

    @app.callback(
        Output(ids.BOX_CHART, "children"),
        Input(ids.RAW_RNA_DATA_DROP, "value"),
        Input(ids.GENE_DROPDOWN, "value"),
        Input(ids.COMPARISON_DROPDOWN, "value"),
    )
    def update_box_chart(dataset_choice: str, gene: str, comps: list[str]) -> html.Div:
        df_filtered = dfs[dataset_choice].query('comparison in @comps')
        fig = px.box(df_filtered, x="comparison", y=gene)
        return html.Div(dcc.Graph(figure=fig), id=ids.BOX_CHART)

    # @app.callback(
    #     Output(ids.COMPARISON_DROPDOWN, "value"),
    #     Input(ids.COMPARISON_DROPDOWN, "value"),
    #     Input(ids.SELECT_ALL_COMPARISONS_BUTTON, "n_clicks")
    # )
    # def select_all_comparisons() -> list[str]:
    #     return list(set(df.comparison))

    return html.Div(
        children=[
            html.H6("Dataset"),
            dcc.Dropdown(
                id=ids.RAW_RNA_DATA_DROP,
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
            html.Div(id=ids.BOX_CHART),
        ],
    )
