from pathlib import Path
from dataclasses import dataclass, field
import pandas as pd


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

    def filter(self, gene: str):
        filtered_data: pd.DataFrame = self.df.query("gene == @gene").T.tail(-1)
        print('filtered_data', filtered_data)
        filtered_data = filtered_data.melt(
            value_vars=list(filtered_data.columns)[:-1],
            id_vars="sample",
            var_name="prot_loc",
            value_name="abun",
        )
        return ProtData(filtered_data)
