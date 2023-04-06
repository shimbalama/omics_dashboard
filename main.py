from dash import Dash
from dash_bootstrap_components.themes import BOOTSTRAP
from src.read_files import RNASeqData
from src.components.layout import create_layout
from pathlib import Path

DATA_PATH = Path("./data").absolute()
bulk = DATA_PATH / "rna" / "bulk"
RNASEQ = [path for path in bulk.glob("*") if not path.name.startswith(".")]
print(RNASEQ)

def main() -> None:
    data = {path.name: RNASeqData(path) for path in RNASEQ}
    app = Dash(external_stylesheets=[BOOTSTRAP])
    app.title = "Omics dashboard"
    app.layout = create_layout(app, data)
    app.run()

if __name__ == "__main__":
    main()
