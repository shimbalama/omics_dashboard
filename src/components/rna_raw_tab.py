from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from . import ids
from src.read_files import RNASeqData
import plotly.express as px


def render(app: Dash, data: dict[str, RNASeqData]) -> html.Div:
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
        selected_data = data[selected_comparison]

        return [{"label": gene, "value": gene} for gene in selected_data.raw_df.columns]

    # 3: get comparisons
    @app.callback(
        Output(ids.COMPARISON_DROPDOWN, "options"),
        Input(ids.RAW_RNA_DATA_DROP, "value"),
    )
    def set_comparison_options(selected_comparison):
        selected_data = data[selected_comparison]
        return [{"label": comp, "value": comp} for comp in selected_data.comparisons]

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
        selected_data = data[dataset_choice]
        print(1212121212, dataset_choice, gene, comps)
        # TODO add FDR
        df_filtered = selected_data.raw_df.query("comparison in @comps")
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
                options=list(data.keys()),
                value="BETi",
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
