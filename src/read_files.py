from pathlib import Path
from dataclasses import dataclass, field
import pandas as pd
from typing import Protocol

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

    def filter(self, gene2: str):
        # print(11111, gene2, self.df.head(), self.df.shape, sep='\n')
        # rows = [gene2, 'gene']
        # filtered_data: pd.DataFrame = self.df.query("gene in @rows").T
        # print('filtered_data', filtered_data, filtered_data.columns, sep='\n')
        return ProtData(self.df.T)
    
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

    def filter(self, gene: str):
        filtered_data: pd.DataFrame = self.df.query("gene == @gene").T.tail(-1)
        # print('filtered_data', filtered_data)
        # filtered_data = filtered_data.melt(
        #     value_vars=list(filtered_data.columns)[:-1],
        #     id_vars="sample",
        #     var_name="prot_loc",
        #     value_name="abun",
        # )
        return ProtData(filtered_data)


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
    processed_dfs: dict[str, pd.DataFrame] = field(repr=False)

    def filter(
        self,
        comparisons: list[str] = None):
        filtered_data: pd.DataFrame = self.df.query(
            "comparison in @comparisons"
        )
        return RNASeqData(filtered_data, self.processed_dfs)
    
    # def __post_init__(self):
    #     """set attribute homopolymer_positions"""
    #     object.__setattr__(self, "name", self.path.name)
    #     object.__setattr__(self, "df", self.read_individual())
    #     object.__setattr__(self, "processed_dfs", self.load_processed_rna_files())


    #TODO, these properties are the same, just dif name in dif files!! they not same

    @property
    def comparisons(self):
        return set(self.df.comparison)
    
    @property
    def degs(self):
        return set(self.processed_dfs.keys())
    
    @property
    def point_of_reference(self):
        '''frgrg'''

        ref_df = self.df.query('point_of_ref == "yes"')
        point_of_ref = set(ref_df.comparison)
        if len(point_of_ref) == 1:
            return point_of_ref.pop()
        elif len(point_of_ref) > 1:
            raise ValueError('Multiple comparisons have POR tag!')
        else:
            return None

    @property
    def mean_count(self) -> pd.DataFrame:
        return self.df.groupby("comparison").agg("mean")

def load_processed_rna_files(path: Path) -> dict[str, pd.DataFrame]:
    """reads processed tsv files"""
    data_dict: dict[str, pd.DataFrame] = {}

    def read_individual(csv: Path) -> pd.DataFrame:
        if csv.suffix.endswith('.csv'):
            delimitor = ','
        elif csv.suffix.endswith('.tsv') or csv.suffix.endswith('.txt'):
            delimitor = '\t'
        else:
            raise ValueError('delimitor unknown!')
        df = pd.read_csv(csv, sep=delimitor)
        df.columns =[
                "gene_id",
                "gene_symbol",
                "gene_biotype",
                "EFFECTSIZE",
                "logCPM",
                "F",
                "P",
                "FDR",
            ]
        
        df = df.reset_index(drop=True)
        return df
    DEGs = path / 'DEGs'
    for csv in DEGs.glob("*"):
        name = str(csv.name).strip().split("_")[0]
        data_dict[name] = read_individual(csv)
    if not data_dict:
        raise FileNotFoundError(f'No DEG data found for {self.path}')

    return data_dict

def read_CPM(path: Path) -> pd.DataFrame:
        # df = pd.read_csv(, index_col=1)
        fin = list(path.glob("*.csv"))
        if fin:
            fin = fin.pop()
            df = pd.read_csv(fin, index_col=0).T.iloc[3:]
            names: list[str] = list(df.index)
            df["comparison"] = [name.split("_")[0] for name in names]
            df["point_of_ref"] = ['yes' if '_POR_' in name else 'no' for name in names]
            return df
        else:
            print(f'log.warn {path} has no data!!!')
            raise FileNotFoundError()
        
def load_RNAseq_data(path: Path) -> RNASeqData:
    return RNASeqData(read_CPM(path), load_processed_rna_files(path))


from functools import partial, reduce
from typing import Callable
from pathlib import Path
import pandas as pd


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


def read_excel(path: Path) -> ProtData:
    # load the data
    files = list(path.glob("*xlsx"))
    assert len(files) == 1
    return pd.read_excel(files[0])

def load_prot_data(path: Path) -> pd.DataFrame:
    df = read_excel(path)
    pipe = [make_gene_col, remove_unwanted_cols, add_categories]
    preprocessor = compose(SchemaProt, *pipe)
    df = preprocessor(df)
    df.to_csv('~/Downloads/prot2222.csv')
    return ProtData(df)

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

def load_phospho_data(path: Path) -> ProtData:
    df = read_excel(path)
    pipe = [make_gene_col, create_unqiue_id, remove_unwanted_cols, add_categories]
    preprocessor = compose(SchemaPhos, *pipe)
    df = preprocessor(df)
    df.to_csv('~/Downloads/phos111.csv')
    return PhosphoProtData(df)
