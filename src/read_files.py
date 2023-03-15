from pathlib import Path
import pandas as pd
from dataclasses import dataclass
from typing import Optional



def load_raw_rna_files(data_path: Path) -> dict[str, pd.DataFrame]:
    """Reads raw RNA csv files"""
    data_dict: dict[str, pd.DataFrame] = {}

    def read_individual(file_path: Path) -> pd.DataFrame:
        df = pd.read_csv(file_path, index_col=1).T.iloc[3:]
        names: list[str] = list(df.index)
        mal_formed_names: list[str] = [
            name for name in names if len(name.split("_")) != 2
        ]
        if mal_formed_names:
            raise ValueError(f"all names should have 1 _, but {mal_formed_names}")
        df["comparison"] = [name.split("_")[0] for name in names]
        return df

    for csv in ["CSexp2_BETi", "CSexp2_shRNA"]:
        data_dict[csv] = read_individual(data_path / f"{csv}_CPM_after_filtering.csv")
    return data_dict

# @dataclass
# class DataSource:
#     _data: pd.DataFrame

#     def filter(
#         self,
#         gene: str,
#         comparisons: Optional[list[str]] = None) -> DataSource:
#         filtered_data: pd.DataFrame = self._data[gene].query(
#             "comparison in @comparisons"
#         )
#         return DataSource(filtered_data)