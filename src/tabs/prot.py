from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from ..components import ids
import dash_bio
from src.read_files import RNASeqData
import pandas as pd
from src.helpers import make_list_of_dicts

KEY = 'prot'


def render(app: Dash, data: dict[str, RNASeqData]) -> html.Div:
    # get comparisons
    def draw_box_chart(df: pd.DataFrame, gene: str, dataset_choice: str) -> html.Div:
        """Draws a box and wisker of the CPM data for each set of replicates for eact
        comparison and overlays the respective FDR value"""
        #RNA_seq_data: RNASeqData = data[KEY][dataset_choice]
        
        fig = px.box(
            df,
            x="comparison",
            y=gene,
            points="all",
            width=999,
            height=666,
            title=f"Boxplot for {gene} CPMs",
            labels={"comparison": "Comparison type", gene: "CPM"},
            facet_row_spacing=0.75
        )

    
    return html.Div(
        children=[
            html.P(
                "You can use Volcano Plot to interactively "
                ""
            ),
            html.H6("Comparison"),
            dcc.Dropdown(
                id=ids.PROCESSED_COMPARISON_DROPDOWN,
                multi=False,
                value=default_comparison,
            ),
            html.H6("P-val"),
            dcc.Slider(
                id=ids.PVALUE,
                value=4,
                max=10,
                min=0,
                step=0.01,
                marks={str(num): str(num) for num in range(0, 11, 2)},
            ),
            html.H6("Effect-size"),
            dcc.RangeSlider(
                id=ids.EEFECT_SIZE,
                min=-4,
                max=4,
                value=[-1, 1],
                step=0.01,
                marks={str(num): str(num) for num in range(-4, 5)},
            ),
            html.Div(
                draw_box_chart(
                    data[BULK][default[0]].raw_df,
                    data[BULK][default[0]].raw_df.columns[0],
                    default[0],
                )
            ),
        ],
    )
