from dash import Dash

# from multiprocessing import Pool
from dash_bootstrap_components.themes import BOOTSTRAP
from src.read_files import (
    load_RNAseq_data,
    load_phospho_data,
    load_prot_data,
    load_function_data,
)
from src.components.layout import create_layout
from pathlib import Path
from src.helpers import rubbish
from time import time
import json

import dash_auth

# Keep this out of source code repository - save in a file or a database
passwords = Path("auth/passwords.json")
with open("auth/passwords.json", "r") as f:
    VALID_USERNAME_PASSWORD_PAIRS = json.load(f)
DATA_PATH = Path("./data").absolute()

READERS = {
    "rna_bulk": load_RNAseq_data,
    "phosphoproteomics": load_phospho_data,
    "proteomics": load_prot_data,
    "function": load_function_data,
}


def read_all(path: Path):  # might have to chnage htis to read on demand??
    return path.name, {
        sub_path.name: (READERS.get(path.name), sub_path)
        for sub_path in path.glob("*")
        if not rubbish(sub_path.name)
    }


def main() -> None:
    start = time()
    data_folders = [path for path in DATA_PATH.glob("*") if not rubbish(path.name)]

    # with Pool(processes=8) as pool:
    #    data = dict(pool.imap_unordered(read_all, data_folders))
    data = dict(read_all(gg) for gg in data_folders)
    print(f"data loaded in {time()-start} seconds")
    app = Dash(external_stylesheets=[BOOTSTRAP])
    auth = dash_auth.BasicAuth(app, VALID_USERNAME_PASSWORD_PAIRS)
    app.title = "Omics dashboard"
    app.layout = create_layout(app, data)
    app.run_server(debug=False, host="0.0.0.0", port=8080)


if __name__ == "__main__":
    main()
