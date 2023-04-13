from dash import Dash
from dash_bootstrap_components.themes import BOOTSTRAP
from src.read_files import RNASeqData
from src.components.layout import create_layout
from pathlib import Path
from src.parse_data.phosphoprot import load_phospho_data, load_prot_data


DATA_PATH = Path("./data").absolute()
bulk = DATA_PATH / "rna" / "bulk"
RNASEQ = [path for path in bulk.glob("*") if not path.name.startswith(".")]


def main() -> None:
    
    data = {"rna_bulk": {path.name: RNASeqData(path) for path in RNASEQ},
            "rna_single": {},
            'prot': load_prot_data(),
            'phospho': load_phospho_data()}
    app = Dash(external_stylesheets=[BOOTSTRAP])
    app.title = "Omics dashboard"
    app.layout = create_layout(app, data)
    app.run()


if __name__ == "__main__":
    main()
