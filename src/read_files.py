from pathlib import Path
from dataclasses import dataclass, field
import pandas as pd
from typing import Protocol
import polars as pl

import numpy as np

# from .helpers import rubbish
from functools import partial, reduce
from typing import Callable


def rubbish(name: str) -> bool:
    return name.startswith((".", "~"))


class Data(Protocol):
    def filter():
        ...


@dataclass(slots=True, frozen=True)
class ProtData:
    """Data pertaining to a chromosomal region
    Parameters
    ----------
    features_to_drop : str or list, default=None
        Variable(s) to be dropped from the dataframe


    Methods
    -------
    find_hom_pol_positions:
        Returns all zero based genomic positions that are in or
        adjacent to homoploymers of given length

    Properties
    -------
    homopolymer_lengths:
        lengths of homolymers
    """

    df: pd.DataFrame = field(repr=False)
    df_FDR: pd.DataFrame = field(repr=False)

    def filter(self, gene2: str):
        return ProtData(self.df, self.df_FDR)


@dataclass(slots=True, frozen=True)
class PhosphoProtData:
    """Data pertaining to a chromosomal region
    Parameters
    ----------
    features_to_drop : str or list, default=None
        Variable(s) to be dropped from the dataframe


    Methods
    -------
    find_hom_pol_positions:
        Returns all zero based genomic positions that are in or
        adjacent to homoploymers of given length

    Properties
    -------
    homopolymer_lengths:
        lengths of homolymers
    """

    df: pd.DataFrame = field(repr=False)
    df_FDR: pd.DataFrame = field(repr=False)

    def filter(self, gene: str):
        dd = self.df.copy().T
        dd["gene"] = dd.index
        dd = dd.query('gene == @gene | gene == "category"')
        del dd["gene"]
        dd = dd.set_index("ID")
        dd = dd.T
        df3 = dd.melt(
            value_vars=list(dd.columns)[:-1],
            id_vars="ID",
            var_name="gene",
            value_name="abun",
        )
        return PhosphoProtData(df3, self.df_FDR)


@dataclass(slots=True, frozen=True)
class RNASeqData:
    """Data pertaining to a chromosomal region
    Parameters
    ----------
    features_to_drop : str or list, default=None
        Variable(s) to be dropped from the dataframe

    Methods
    -------
    find_hom_pol_positions:
        Returns all zero based genomic positions that are in or
        adjacent to homoploymers of given length

    Properties
    -------
    homopolymer_lengths:
        lengths of homolymers
    """

    # path: Path
    # name: str = field(init=False)
    df: pd.DataFrame = field(repr=False)
    df_FDR: pd.DataFrame = field(repr=False)
    processed_dfs: dict[str, pd.DataFrame] = field(repr=False)

    def filter(self, comparisons: list[str] = None):
        df = self.df.copy(deep=True)
        print(99999, comparisons)
        filtered_data: pd.DataFrame = df.query("comparison in @comparisons")
        return RNASeqData(
            filtered_data.copy(deep=True),
            self.df_FDR.copy(deep=True),
            self.processed_dfs,
        )

    @property
    def comparisons(self):
        return set(self.df["comparison"])

    @property
    def degs(self):
        return set(self.processed_dfs.keys())

    @property
    def point_of_reference(self):
        """frgrg"""

        ref_df = self.df.query('point_of_ref == "yes"')
        point_of_ref = set(ref_df.comparison)
        if len(point_of_ref) == 1:
            return point_of_ref.pop()
        elif len(point_of_ref) > 1:
            raise ValueError("Multiple comparisons have POR tag!")
        else:
            return None

    @property
    def mean_count(self) -> pd.DataFrame:
        return self.df.groupby("comparison").agg("mean")


def load_processed_rna_files(path: Path) -> dict[str, pd.DataFrame]:
    """reads processed tsv files"""
    data_dict: dict[str, pd.DataFrame] = {}

    def read_individual(csv: Path) -> pd.DataFrame:
        if csv.suffix.endswith(".csv"):
            delimitor = ","
        elif csv.suffix.endswith(".tsv") or csv.suffix.endswith(".txt"):
            delimitor = "\t"
        else:
            raise ValueError("delimitor unknown!")
        df = pl.read_csv(csv, separator=delimitor).to_pandas(
            use_pyarrow_extension_array=True
        )
        df.columns = [
            "gene_id",
            "gene_symbol",
            "gene_biotype",
            "EFFECTSIZE",
            "logCPM",
            "F",
            "P",
            f"{csv.stem}_FDR",
        ]

        df = df.reset_index(drop=True)
        return df

    DEGs = path / "DEGs"
    for csv in DEGs.glob("*"):
        # name = str(csv.name).strip().split("_")[0]
        data_dict[csv.stem] = read_individual(csv)
    if not data_dict:
        raise FileNotFoundError(f"No DEG data found for {self.path}")
    return data_dict


def read_CPM(path: Path) -> pd.DataFrame:
    # df = pd.read_csv(, index_col=1)
    fin = list(path.glob("*.csv"))
    if fin:
        fin = fin.pop()
        df = pl.read_csv(str(fin)).to_pandas(use_pyarrow_extension_array=True)
        df = df.replace("NA", np.NaN)
        df["index"] = df.gene_symbol.fillna(df.gene_id)
        df = df.set_index("index")
        df = df.T.iloc[3:]
        names: list[str] = list(df.index)
        df["comparison"] = [name.split("_")[0] for name in names]
        df["point_of_ref"] = ["yes" if "_POR_" in name else "no" for name in names]
        return df
    else:
        print(f"log.warn {path} has no data!!!")
        raise FileNotFoundError()


# def merge_FDR(fdrs):
#     fdrs2 = [df[["gene_id", f"{name}_FDR"]] for name, df in fdrs.items()]
#     final_df = reduce(
#         lambda left, right: pd.merge(left, right, on=["gene_id"], how="outer"), fdrs2
#     )
#     final_df.columns = ["gene_id"] + list(fdrs.keys())
#     return final_df


def merge_FDR(fdrs):
    genes = list(fdrs.values())[0]["gene_id"]
    return pd.concat([genes] + [df[f"{name}_FDR"] for name, df in fdrs.items()], axis=1)


def load_RNAseq_data(path: Path) -> RNASeqData:
    CPM = read_CPM(path)
    FDR = load_processed_rna_files(path)
    merged_FDRs = merge_FDR(FDR)
    CPM.to_csv("~/Downloads/CPM3333.csv")
    merged_FDRs.to_csv("~/Downloads/FDR44444.csv")
    return RNASeqData(CPM, merged_FDRs, FDR)


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


def remove_unwanted_cols(
    schema: type, df: pd.DataFrame, misc: list[str]
) -> pd.DataFrame:
    abundances: list[str] = [
        col
        for col in list(df.columns)
        if col.startswith(("Abundance: F", "Abundance Ratio Adj. P-Value:"))
    ]
    return df[misc + abundances]


def add_categories(schema: type, df: pd.DataFrame) -> pd.DataFrame:
    df = df.set_index(schema.GENE)
    df = df.T
    index_str = "index_str"
    df[index_str] = df.index
    df[schema.CAT] = df[index_str].apply(
        lambda x: "P_value" if "P-Value:" in x else str(x).split(",")[-1]
    )
    df[schema.CAT] = df[schema.CAT].astype("category")
    del df[index_str]

    return df.T


def read_excel(path: Path) -> Data:
    # load the data
    files = [fin for fin in list(path.glob("*xlsx")) if not rubbish(fin.name)]
    assert len(files) == 1
    return pd.read_excel(files[0])


def split_dfs(df):
    df = df.copy(deep=True).T
    df_FDR = df.query('category == "P_value"')
    df = df.query('category != "P_value"')

    return df, df_FDR


def load_prot_data(path: Path) -> ProtData:
    df = read_excel(path)
    pipe = [
        make_gene_col,
        partial(remove_unwanted_cols, misc=[SchemaProt.GENE]),
        add_categories,
    ]
    preprocessor = compose(SchemaProt, *pipe)
    df = preprocessor(df)
    df, df_FDR = split_dfs(df)
    df.to_csv("~/Downloads/prot2222.csv")
    return ProtData(df, df_FDR)


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


def load_phospho_data(path: Path) -> PhosphoProtData:
    df = read_excel(path)
    pipe = [
        make_gene_col,
        create_unqiue_id,
        partial(remove_unwanted_cols, misc=[SchemaPhos.GENE, SchemaPhos.UNQIUE_ID]),
        add_categories,
    ]
    preprocessor = compose(SchemaPhos, *pipe)
    df = preprocessor(df)
    df, df_FDR = split_dfs(df)
    df.to_csv("~/Downloads/phos111.csv")
    return PhosphoProtData(df, df_FDR)


# load_prot_data(Path("/Users/liam/code/omics_dashboard/data/proteomics/FibrosisStim/"))
# load_phospho_data(
#     Path("/Users/liam/code/omics_dashboard/data/phosphoproteomics/ET1-His")
# )

# load_RNAseq_data(Path('/Users/liam/code/omics_dashboard/data/rna_bulk/Reid_unpub_Cytokine_Storm_vs_BETi'))
