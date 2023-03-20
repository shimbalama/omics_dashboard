from pathlib import Path
import pandas as pd
from dataclasses import dataclass
from typing import Optional


def load_raw_rna_files(data_path: Path) -> dict[str, pd.DataFrame]:
    """Reads raw RNA csv files"""
    data_dict: dict[str, pd.DataFrame] = {}

    def read_individual(file_path: Path) -> pd.DataFrame:
        #df = pd.read_csv(file_path, index_col=1)
        df = pd.read_csv(file_path, index_col=1).T.iloc[3:]
        names: list[str] = list(df.index)
        mal_formed_names: list[str] = [
            name for name in names if len(name.split("_")) != 2
        ]
        if mal_formed_names:
            raise ValueError(f"all names should have 1 _, but {mal_formed_names}")
        df["comparison"] = [name.split("_")[0] for name in names]
        #print(1, df.head(), df.dtypes, sep="\n")
        return df

    for csv in ["CSexp2_BETi", "CSexp2_shRNA"]:
        data_dict[csv] = read_individual(data_path / f"{csv}_CPM_after_filtering.csv")
    return data_dict


def load_processed_rna_files(data_path: Path) -> dict[str, pd.DataFrame]:
    """reads processed tsv files"""
    data_dict: dict[str, pd.DataFrame] = {}

    def read_individual(csv: Path) -> pd.DataFrame:
        df = pd.read_csv(csv, sep="\t", index_col=0)
        df['EFFECTSIZE'] = df.logFC.copy()
        df['P'] = df.PValue.copy()
        df = df.reset_index(drop=True)
        return df

    for csv in data_path.glob("*.txt"):
        name = str(csv.name).strip().split("_")[0]
        data_dict[name] = read_individual(csv)

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
#         return DataSource(filtered_data) 239 ES6
