from dash import Dash
from dash_bootstrap_components.themes import BOOTSTRAP
from src.read_files import RNASeqData
from src.components.layout import create_layout
from pathlib import Path

DATA_PATH = Path("./data").absolute()
RNASEQ = ['BETi']


def main() -> None:
    data = {name: RNASeqData(name, DATA_PATH) for name in RNASEQ}
    app = Dash(external_stylesheets=[BOOTSTRAP])
    app.title = "Omics dashboard"
    app.layout = create_layout(app, data)
    app.run()


if __name__ == "__main__":
    main()
