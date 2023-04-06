from pathlib import Path
from dataclasses import dataclass, field
import pandas as pd



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

    path: Path
    name: str = field(init=False)
    raw_df: pd.DataFrame = field(init=False, repr=False)
    #raw_mean_df: pd.DataFrame = field(init=False, repr=False)
    processed_dfs: dict[str, pd.DataFrame] = field(init=False, repr=False)

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
        object.__setattr__(self, "name", self.path.name)
        object.__setattr__(self, "raw_df", self.read_individual())
        object.__setattr__(self, "processed_dfs", self.load_processed_rna_files())

    def read_individual(self) -> pd.DataFrame:
        # df = pd.read_csv(, index_col=1)
        print(2222222222,self.path, list(self.path.glob('*')), list(self.path.glob("*.csv")), sep='\n')
        fin = list(self.path.glob("*.csv"))
        if fin:
            fin = fin.pop()
            df = pd.read_csv(fin, index_col=0).T.iloc[3:]
            names: list[str] = list(df.index)
            # mal_formed_names: list[str] = [
            #     name for name in names if len(name.split("_")) != 2
            # ]
            # if mal_formed_names:
            #     raise ValueError(f"all names should have 1 _, but {mal_formed_names}")
            df["comparison"] = [name.split("_")[0] for name in names]
            return df
        else:
            print(f'log.warn {self.path} has no data!!!')
            raise FileNotFoundError()

    #TODO, these properties are the same, just dif name in dif files!!
    @property
    def comparisons(self):
        return set(self.raw_df.comparison)
    
    @property
    def degs(self):
        return set(self.processed_dfs.keys())

    @property
    def mean_count(self) -> pd.DataFrame:
        return self.raw_df.groupby("comparison").agg("mean")

    def load_processed_rna_files(self) -> dict[str, pd.DataFrame]:
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
        DEGs = self.path / 'DEGs'
        for csv in DEGs.glob("*"):
            name = str(csv.name).strip().split("_")[0]
            data_dict[name] = read_individual(csv)
        if not data_dict:
            raise FileNotFoundError(f'No DEG data found for {self.path}')

        return data_dict
