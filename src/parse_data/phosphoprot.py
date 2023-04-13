from functools import partial, reduce
from typing import Callable

import pandas as pd


# TODO read all data with MP pool


class SchemaProt:
    COMBINED_PROT_DESCRIPTION = "Description"
    GENE = "gene"
    CAT = "category"


class SchemaPhos:
    UNQIUE_ID = "ID"
    POS_IN_PROT = "Positions in Proteins"
    MOD_IN_PROT = "Modifications in Proteins"
    COMBINED_PROT_DESCRIPTION = "Master Protein Descriptions"
    GENE = "gene"
    CAT = "category"


Preprocessor = Callable[[type, pd.DataFrame], pd.DataFrame]


def compose(schema: type, *functions: Preprocessor) -> Preprocessor:
    """Helper function to call all df functions sequencially"""
    partially_filled_funcs = [partial(func, schema) for func in functions]
    return reduce(
        lambda func1, func2: lambda x: func2(func1(x)), partially_filled_funcs
    )


def make_gene_col(schema: type, df: pd.DataFrame) -> pd.DataFrame:
    df[schema.GENE] = df[schema.COMBINED_PROT_DESCRIPTION].apply(
        lambda x: [val for val in x.strip().split(" ") if "GN=" in val]
        .pop()
        .replace("GN=", "")
        if "GN=" in x
        else x.split(" ")[-1]
    )
    return df


def remove_unwanted_cols(schema: type, df: pd.DataFrame) -> pd.DataFrame:
    abundances: list[str] = [
        col
        for col in list(df.columns)
        if col.startswith(("Abundance: F", "Abundance Ratio Adj. P-Value:"))
    ]
    misc = [schema.GENE]
    return df[misc + abundances]


def add_categories(schema: type, df: pd.DataFrame) -> pd.DataFrame:
    df = df.T
    index_str = "index_str"
    df[index_str] = df.index
    df[schema.CAT] = df[index_str].apply(
        lambda x: "P_value" if "P-Value:" in x else str(x).split(",")[-1]
    )
    df[schema.CAT] = df[schema.CAT].astype("category")
    del df[index_str]

    return df.T


def load_prot_data() -> pd.DataFrame:
    # load the data from the CSV file
    data = pd.read_excel(
        "/Users/liam/code/data/Finalised_files_for_dashboard/Proteomics/RawoutputPD_FibrosisStim.xlsx",
    )
    pipe = [make_gene_col, remove_unwanted_cols, add_categories]
    preprocessor = compose(SchemaProt, *pipe)
    return preprocessor(data)


def create_unqiue_id(schema: type, df: pd.DataFrame) -> pd.DataFrame:
    def choose_pos(gene: str, pos_extact: str | float, pos_range: str) -> str:
        if not isinstance(pos_extact, float):
            poses = pos_extact.strip().split(" ")[-1]
            pos: str = "_".join(bit.split("(")[0][1:] for bit in poses.split(";"))
        else:
            pos: str = pos_range.split(" ")[-1][1:-1]
        return "_".join([gene, pos])

    df[schema.UNQIUE_ID] = df.apply(
        lambda x: choose_pos(
            x[schema.GENE],
            x[schema.MOD_IN_PROT],
            x[schema.POS_IN_PROT],
        ),
        axis=1,
    )
    return df


def load_phospho_data() -> pd.DataFrame:
    # load the data from the CSV file
    df = pd.read_excel(
        "/Users/liam/code/data/Finalised_files_for_dashboard/Phosphoproteomics/Raw_Phosphoproteomics_ET1-His.xlsx",
    )
    pipe = [make_gene_col, create_unqiue_id, remove_unwanted_cols, add_categories]
    preprocessor = compose(SchemaPhos, *pipe)

    return preprocessor(df)
