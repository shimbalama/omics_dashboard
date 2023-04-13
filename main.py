from dash import Dash
from dash_bootstrap_components.themes import BOOTSTRAP
from src.read_files import RNASeqData
from src.components.layout import create_layout
from pathlib import Path
from src.parse_data.phosphoprot import load_phospho_data, load_prot_data


DATA_PATH = Path("./data").absolute()
bulk = DATA_PATH / "rna" / "bulk"
RNASEQ = [path for path in bulk.glob("*") if not path.name.startswith(".")]
print(RNASEQ)

def main() -> None:
    print(111, load_prot_data())
    print(2222,load_phospho_data())
    data = {path.name: RNASeqData(path) for path in RNASEQ}
    app = Dash(external_stylesheets=[BOOTSTRAP])
    app.title = "Omics dashboard"
    app.layout = create_layout(app, data)
    app.run()

if __name__ == "__main__":
    main()
