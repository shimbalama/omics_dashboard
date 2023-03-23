from pathlib import Path
import pandas as pd
from dataclasses import dataclass, field
from typing import Optional


@dataclass(slots=True, frozen=True)
class RNASeqData:
    """Data pertaining to a chromosomal region
    Parameters
    ----------
    features_to_drop : str or list, default=None
        Variable(s) to be dropped from the dataframe
    chrom: str, default=None
        the chromosome number

    Methods
    -------
    find_hom_pol_positions:
        Returns all zero based genomic positions that are in or
        adjacent to homoploymers of given length
    get_hom_pol_lengths:
        Adds the length of the associated homolymer to to each zero based position
        in homopolymer_positions
    Properties
    -------
    homopolymer_lengths:
        lengths of homolymers
    """
    name: str
    path: Path
    raw: pd.DataFrame = field(init=False, repr=False)
    processed: dict[str, pd.DataFrame] = field(init=False, repr=False)

    # def filter(
    #     self,
    #     gene: str,
    #     comparisons: Optional[list[str]] = None) -> DataSource:
    #     filtered_data: pd.DataFrame = self._data[gene].query(
    #         "comparison in @comparisons"
    #     )
    #     return DataSource(filtered_data)
    def __post_init__(self):
        """set attribute homopolymer_positions"""

        object.__setattr__(self, "raw", self.read_individual())
        object.__setattr__(self, "processed", self.load_processed_rna_files())
    
    def read_individual(self) -> pd.DataFrame:
        #df = pd.read_csv(, index_col=1)
        fin = list(self.path.glob(f'*{self.name}*.csv')).pop()
        
        df = pd.read_csv(fin, index_col=1).T.iloc[3:]
        names: list[str] = list(df.index)
        mal_formed_names: list[str] = [
            name for name in names if len(name.split("_")) != 2
        ]
        if mal_formed_names:
            raise ValueError(f"all names should have 1 _, but {mal_formed_names}")
        df["comparison"] = [name.split("_")[0] for name in names]
        print(1, df.head(2), df.dtypes, sep="\n")
        return df

    def load_processed_rna_files(self) -> dict[str, pd.DataFrame]:
        """reads processed tsv files"""
        data_dict: dict[str, pd.DataFrame] = {}

        def read_individual(csv: Path) -> pd.DataFrame:
            df = pd.read_csv(csv, sep="\t", index_col=0)
            df['EFFECTSIZE'] = df.logFC.copy()
            df['P'] = df.PValue.copy()
            df = df.reset_index(drop=True)
            print(1222, df.head(), df.dtypes, sep="\n")
            return df

        for csv in self.path.glob("*.txt"):
            name = str(csv.name).strip().split("_")[0]
            data_dict[name] = read_individual(csv)

        return data_dict



