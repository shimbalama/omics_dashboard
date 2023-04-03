from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from . import ids
from src.read_files import RNASeqData
import plotly.express as px
import pandas as pd

from src.helpers import make_list_of_dicts


def render(app: Dash, data: dict[str, RNASeqData]) -> html.Div:
    # see https://dash.plotly.com/basic-callbacks#dash-app-with-chained-callbacks

    def draw_box_chart(df: pd.DataFrame, gene: str, dataset_choice: str) -> html.Div:
        """Draws a box and wisker of the CPM data for each set of replicates for eact
        comparison and overlays the respective FDR value"""
        def draw_line(fig, x_start, x_end, y, height, plot_type='line', ref='paper', line_vals = {'color':"magenta",'width':3}):
            # Top line
            fig.add_shape(type=plot_type,
                xref=ref, yref=ref,
                x0=x_start, y0=y+height,
                x1=x_end, y1=y+height,
                line=line_vals)

            # Left line
            fig.add_shape(type=plot_type,
                xref=ref, yref=ref,
                x0=x_start, y0=y,
                x1=x_start, y1=y+height,
                line=line_vals)

            # Right line
            fig.add_shape(type=plot_type,
                xref=ref, yref=ref,
                x0=x_end, y0=y+height,
                x1=x_end, y1=y,
                line=line_vals)
            
            return fig
        
        fig = px.box(
            df,
            x="comparison",
            y=gene,
            points="all",
            width=1111,
            height=888,
            title=f"Boxplot for {gene} CPMs",
            labels={"comparison": "Comparison type", gene: "CPM"},
        )

        for i, comp in enumerate(df.comparison.unique(), start=1):
            DEG_df: pd.DataFrame | None = data[dataset_choice].processed_dfs.get(comp)
            if DEG_df is not None:
                FDR = float(DEG_df.query("gene_id == @gene").FDR.iloc[0])
            else:
                FDR = 0.0
            y = df.query("comparison == @comp")[gene].median()
            fig.add_annotation(
                x=comp,
                y=y,
                text=f"{FDR:.1e}",
                yshift=10,
                showarrow=False,
            )
            i /= 9
            print(111111,i, y)
            fig = draw_line(fig, i, i+0.2, 0.6, 0.6+i, plot_type='line', ref='paper', line_vals = {'color':"magenta",'width':3})
        fig.add_shape(type="line",
            xref="paper", yref="paper",
            x0=0.1, y0=0.8,
            x1=0.43, y1=0.8,
            line=dict(
                color="magenta",
                width=3,
            ),
                    )

        return html.Div(dcc.Graph(figure=fig), id=ids.BOX_CHART)

    @app.callback(
        Output(ids.GENE_DROPDOWN, "options"), Input(ids.RAW_RNA_DATA_DROP, "value")
    )
    def set_gene_options(experiment: str) -> list[dict[str, str]]:
        """Populates the gene selection dropdown with options from teh given dataset"""
        return make_list_of_dicts(list(data[experiment].raw_df.columns))

    @app.callback(
        Output(ids.COMPARISON_DROPDOWN, "options"),
        Input(ids.RAW_RNA_DATA_DROP, "value"),
    )
    def set_comparison_options(experiment: str) -> list[dict[str, str]]:
        """Populates the comparison selection dropdown with options from teh given dataset"""
        return make_list_of_dicts(list(data[experiment].comparisons))

    @app.callback(
        Output(ids.GENE_DROPDOWN, "value"), Input(ids.GENE_DROPDOWN, "options")
    )
    def select_gene_value(gene_options: list[dict[str, str]]) -> str:
        """Select first gene as default value"""
        return gene_options[0]["value"]

    @app.callback(
        Output(ids.COMPARISON_DROPDOWN, "value"),
        Input(ids.COMPARISON_DROPDOWN, "options"),
        Input(ids.SELECT_ALL_COMPARISONS_BUTTON, "n_clicks"),
    )
    def select_comparison_values(
        available_comparisons: list[dict[str, str]], _: int
    ) -> list[dict[str, str]]:
        """Default to all available comparisons"""
        return [comp["value"] for comp in available_comparisons]

    @app.callback(
        Output(ids.BOX_CHART, "children"),
        Input(ids.RAW_RNA_DATA_DROP, "value"),
        Input(ids.GENE_DROPDOWN, "value"),
        Input(ids.COMPARISON_DROPDOWN, "value"),
    )
    def update_box_chart(dataset_choice: str, gene: str, comps: list[str]) -> html.Div:
        """Re draws a box and wisker of the CPM data for each set of replicates for eact
        comparison and overlays the respective FDR value"""
        selected_data = data[dataset_choice]
        df_filtered = selected_data.raw_df.query("comparison in @comps")

        return draw_box_chart(df_filtered, gene, dataset_choice)

    default = list(data.keys())
    return html.Div(
        children=[
            html.H6("Dataset"),
            dcc.Dropdown(
                id=ids.RAW_RNA_DATA_DROP,
                options=default,
                value=default[0],
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
            html.Div(
                draw_box_chart(
                    data[default[0]].raw_df,
                    data[default[0]].raw_df.columns[0],
                    default[0],
                )
            ),
        ],
    )
