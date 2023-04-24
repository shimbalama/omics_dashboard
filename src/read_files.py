from pathlib import Path
from dataclasses import dataclass, field
import pandas as pd
from typing import Protocol
import polars as pl
from time import time
import numpy as np
import numexpr as ne
import polars as pl
from collections import Counter

# from .helpers import rubbish
from functools import partial, reduce
from typing import Callable


def rubbish(name: str) -> bool:
    return name.startswith((".", "~"))


class Data(Protocol):
    def filter():
        ...

    def point_of_reference():
        ...

    def pandas_df():
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

    name: str
    df: pl.DataFrame = field(repr=False)
    df_FDR: pd.DataFrame = field(repr=False)

    def filter(self, gene: str, tests: list[str]):
        df = self.df.filter(self.df["test"].is_in(tests))
        df = df.select([gene, 'test'])
        return ProtData(self.name, df, self.df_FDR)

    @property
    def point_of_reference(self):
        return "CTRL"

    @property
    def pandas_df(self):
        return self.df.to_pandas()
    
    @property
    def test_names(self):
        return set(self.df["test"])


@dataclass(slots=True, frozen=True)
class PhosphoProtData:#wait for more data before tightening bolts here... 
    #need to figure if want FDR bracket and if filter needs 'tests'
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

    name: str
    df: pd.DataFrame = field(repr=False)
    df_FDR: pd.DataFrame = field(repr=False)

    def filter(self, gene: str, tests: list[str]):
        df = self.df.copy().T
        df["gene"] = df.index
        df = df.query('gene == @gene | gene == "test"')
        del df["gene"]
        df = df.set_index("ID")
        df = df.T
        df_melted = df.melt(
            value_vars=list(df.columns)[:-1],
            id_vars="ID",
            var_name="gene",
            value_name="abun",
        )
        return PhosphoProtData(self.name, df_melted, self.df_FDR)

    @property
    def point_of_reference(self):
        return "CTRL"

    @property
    def pandas_df(self):
        return self.df
    
    @property
    def test_names(self):
        return set(self.df["test"])


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

    name: str
    df: pl.DataFrame = field(repr=False)
    df_FDR: pd.DataFrame = field(repr=False)
    processed_dfs: dict[str, pd.DataFrame] = field(repr=False)
    point_of_reference: str = field(init=False)

    def __post_init__(self):
        """set attribute point_of_reference"""

        object.__setattr__(self, "point_of_reference", self.find_point_of_reference())

    def filter(self, gene: str, tests: list[str] = None):
        df = self.df.filter(self.df["test"].is_in(tests))
        df = df.select([gene, 'test' , 'point_of_ref'])
        return RNASeqData(
            self.name,
            df,
            self.df_FDR,
            self.processed_dfs,
        )

    # TODO - post init to check that indexes are equal. found examples of missign and duplicates gene names

    @property
    def pandas_df(self):
        return self.df.to_pandas()

    @property
    def test_names(self):
        return set(self.df["test"])

    @property
    def degs(self):
        return set(self.processed_dfs.keys())

    def find_point_of_reference(self):
        """frgrg"""
        
        ref_df = self.df.filter(pl.col("point_of_ref") == "yes")
        point_of_ref = set(ref_df["test"])
        if len(point_of_ref) == 1:
            return point_of_ref.pop()
        elif len(point_of_ref) > 1:
            raise ValueError("Multiple comparisons have POR tag!")
        else:
            return None

    # @property
    # def mean_count(self) -> pd.DataFrame:
    #     return self.df.groupby("test").agg("mean")


def make_index(df):
    df.gene_symbol = df.gene_symbol.replace("NA", np.NaN)
    df.gene_symbol = df["gene_symbol"].mask(
        df["gene_symbol"].duplicated(keep=False), other=np.NaN
    )
    df["index2"] = df.gene_symbol.fillna(df.gene_id)
    df = df.set_index("index2")
    idx = df.index
    if any(idx.duplicated()):
        print("log.warn!", "dupppssss!!!!", df[idx.duplicated()])
    return df


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
            csv.stem,
        ]
        df = make_index(df)
        return df

    DEGs = path / "DEGs"
    for csv in DEGs.glob("*"):
        data_dict[csv.stem] = read_individual(csv)
    if not data_dict:
        raise FileNotFoundError(f"No DEG data found for {self.path}")
    return data_dict


def read_CPM(path: Path) -> pd.DataFrame:
    fin = list(path.glob("*.csv"))
    if fin:
        fin = fin.pop()
        df = pl.read_csv(str(fin)).to_pandas(use_pyarrow_extension_array=True)
        df = make_index(df)
        df = df.T.iloc[3:]
        names: list[str] = list(df.index)
        df["test"] = [name.split("_")[0] for name in names]
        df["point_of_ref"] = ["yes" if "_POR_" in name else "no" for name in names]
        return df
    else:
        print(f"log.warn {path} has no data!!!")
        raise FileNotFoundError()


def merge_FDR(fdrs):

    df = pd.concat([df[name] for name, df in fdrs.items()], axis=1).T
    df["test"] = df.index

    return df


def load_RNAseq_data(path: Path) -> RNASeqData:
    CPM = read_CPM(path)
    CPM = CPM.query('test != "description"')
    CPM = CPM.astype({col: "Float32" for col in CPM.columns[:-2]})
    FDR = load_processed_rna_files(path)
    merged_FDRs = merge_FDR(FDR)
    CPM = pl.from_pandas(CPM)
    # CPM.to_csv("~/Downloads/CPM3333.csv")
    # merged_FDRs.to_csv("~/Downloads/FDR44444.csv")
    return RNASeqData(path.stem, CPM, merged_FDRs, FDR)


class SchemaProt:
    COMBINED_PROT_DESCRIPTION = "Description"
    GENE = "gene"
    CAT = "test"


class SchemaPhos:
    UNQIUE_ID = "ID"
    POS_IN_PROT = "Positions in Proteins"
    MOD_IN_PROT = "Modifications in Proteins"
    COMBINED_PROT_DESCRIPTION = "Master Protein Descriptions"
    GENE = "gene"
    CAT = "test"


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


def make_gene_col_unique(schema: type, df: pd.DataFrame) -> pd.DataFrame:
    non_unique_genes = list(df[schema.GENE])
    counts = Counter(non_unique_genes)
    unique_genes = []
    for i, gene in enumerate(non_unique_genes):
        if counts[gene] > 1:
            unique_genes.append(
                f"{gene}{counts[gene]-non_unique_genes[:i].count(gene)}"
            )
        else:
            unique_genes.append(gene)

    # genes = list(df[schema.GENE])
    # unique_genes = [ f"{a}{c[a]}" for c in [Counter()] for a in genes if [c.update([a])] ]
    df[schema.GENE] = unique_genes
    return df


def remove_unwanted_cols(
    schema: type, df: pd.DataFrame, misc: list[str]
) -> pd.DataFrame:
    abundances: list[str] = [
        col
        for col in list(df.columns)
        if col.startswith(("Abundance: F", "Abundance Ratio Adj. P-Value:"))
    ]
    df = df[misc + abundances]

    return df


def add_categories(schema: type, df: pd.DataFrame) -> pd.DataFrame:
    df = df.set_index(schema.GENE)
    df = df.T
    index_str = "index_str"
    df[index_str] = df.index
    df[schema.CAT] = df[index_str].apply(
        lambda x: "P_value" if "P-Value:" in x else str(x).split(",")[-1].strip()
    )
    df[schema.CAT] = df[schema.CAT].astype("category")
    del df[index_str]

    return df.T


def read_excel(path: Path) -> Data:
    # load the data
    files = [fin for fin in list(path.glob("*xlsx")) if not rubbish(fin.name)]
    assert len(files) == 1
    return pd.read_excel(files[0])
    # pls caused too many dramas ie The current offset in the file is 3909947 bytes.
    # can't handle values like 0.001000


def split_dfs(df):
    df = df.T
    df_FDR = df.query('test == "P_value"')
    df_FDR["index_cpy"] = df_FDR.index
    df_FDR["test"] = df_FDR.index_cpy.apply(
        lambda x: x.strip()
        .split(") /")[0]
        .replace("Abundance Ratio Adj. P-Value: (", "")
    )
    del df_FDR["index_cpy"]
    df = df.query('test != "P_value"')

    return df, df_FDR


def load_prot_data(path: Path) -> ProtData:
    df = read_excel(path)
    pipe = [
        make_gene_col,
        make_gene_col_unique,
        partial(remove_unwanted_cols, misc=[SchemaProt.GENE]),
        add_categories,
    ]
    preprocessor = compose(SchemaProt, *pipe)
    df = preprocessor(df)

    df, df_FDR = split_dfs(df)
    # df.to_csv("~/Downloads/prot2222.csv")
    return ProtData(path.stem, pl.from_pandas(df), df_FDR)


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
        make_gene_col,  # make_gene_col_unique, can't as are mainnly dups..
        create_unqiue_id,
        partial(remove_unwanted_cols, misc=[SchemaPhos.GENE, SchemaPhos.UNQIUE_ID]),
        add_categories,
    ]
    preprocessor = compose(SchemaPhos, *pipe)
    df = preprocessor(df)
    df, df_FDR = split_dfs(df)
    # df.to_csv("~/Downloads/phos111.csv")
    return PhosphoProtData(path.stem, df, df_FDR)


# load_prot_data(Path("/Users/liam/code/omics_dashboard/data/proteomics/FibrosisStim/"))
# load_phospho_data(
#     Path("/Users/liam/code/omics_dashboard/data/phosphoproteomics/ET1-His")
# )

# load_RNAseq_data(Path('/Users/liam/code/omics_dashboard/data/rna_bulk/Reid_unpub_Cytokine_Storm_vs_BETi'))
