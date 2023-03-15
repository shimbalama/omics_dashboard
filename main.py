from dash import Dash
from dash_bootstrap_components.themes import BOOTSTRAP
from src.read_files import load_raw_rna_files
from src.components.layout import create_layout
from pathlib import Path

DATA_PATH = Path("./data").absolute()


def main() -> None:
    dfs = load_raw_rna_files(DATA_PATH)#dict of functions? 2do?
    app = Dash(external_stylesheets=[BOOTSTRAP])
    app.title = "Omics dashboard"
    app.layout = create_layout(app, dfs)
    app.run()


if __name__ == "__main__":
    main()
