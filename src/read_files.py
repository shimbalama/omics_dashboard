from pathlib import Path
from dataclasses import dataclass, field
import pandas as pd
from typing import Protocol
import polars as pl
import numpy as np
from collections import Counter, defaultdict
from abc import ABC, abstractmethod, abstractproperty

from functools import partial, reduce
from typing import Callable


def rubbish(name: str) -> bool:
    return name.startswith((".", "~"))


class Data(Protocol):
    def filter():
        """Returns a filtred version of the data"""
        ...

    def point_of_reference():  # still?TODO
        ...

    def plot_df():
        ...


@dataclass(slots=True, frozen=True)
class Data2(ABC):
    name: str
    gene: str
    ordered_test_names: list[str]
    stats: dict[str, bool]

    @abstractmethod
    def filter(self):  # prot, phos
        """Returns a filtred version of the data"""
        pass

    @abstractproperty
    def plot_df(self):  # prot, phos
        pass

    @abstractproperty
    def test_names(self):  # prot, phos
        pass

    @abstractproperty
    def point_of_reference(self):  # still?TODO #prot, phos
        pass

    def get_FDR(self, test, prot_pos=None):
        try:
            return self.df_FDR[test]
        except KeyError:
            return 0.0

    def get_median_CPMs(self, test, prot_pos=None):
        return np.median(self.plot_df.query("test == @test")[self.gene])


@dataclass(slots=True, frozen=True)
class ProtData(Data2):
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

    df: pl.DataFrame = field(repr=False)
    df_FDR: pd.DataFrame = field(repr=False)
    # stats: bool = True
    # gene: str | None = None

    def filter(self, gene: str, tests: list[str], stats: dict[str, bool]) -> "ProtData":
        df = self.df.filter(self.df["test"].is_in(tests))
        df = df.select([gene, "test"])
        return ProtData(
            name=self.name,
            df=df,
            df_FDR=self.df_FDR[gene],
            gene=gene,
            ordered_test_names=[test for test in tests if test != "ID"],
            stats=stats,
        )

    @property
    def point_of_reference(self):
        return "CTRL"

    @property
    def plot_df(self):
        return self.df.to_pandas()

    @property
    def test_names(self):
        return sorted(np.unique(self.df[Schema.CAT]))


@dataclass(slots=True, frozen=True)
class PhosphoProtData(Data2):  # wait for more data before tightening bolts here...
    # need to figure if want FDR bracket and if filter needs 'tests'
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
    ###TODO DataFrame columns are not unique, some columns will be omitted.UserWarning
    # name: str
    df: pd.DataFrame = field(repr=False)
    df_FDR: pd.DataFrame = field(repr=False)
    # stats: bool = True
    # gene: str | None = None

    def filter(self, gene: list[str], tests: list[str], stats: dict[str, bool]):  # , phos_sites: list[str]):
        # gene is list of phospho sites and gene name from from them
        phos_sites = gene[:]
        gene = gene.pop().split("_")[0]
        if "ID" not in tests:
            tests.append("ID")
        df = self.df.copy()
        df = df.query("test == @tests").T
        df[Schema.GENE] = df.index
        df = df.query('gene == @gene | gene == "test"')
        df = df.query('ID == @phos_sites | ID == "ID"')

        del df[Schema.GENE]

        df = df.set_index(Schema.UNQIUE_ID)
        df = df.T
        df_melted = df.melt(
            value_vars=list(df.columns)[:-1],
            id_vars=Schema.UNQIUE_ID,
            var_name=Schema.GENE,
            value_name=gene,
        )
        df_melted.columns = ["test", Schema.GENE, gene]
        for phos_pos, pos_df in df_melted.groupby("gene"):
            if not any(pos_df.query("test == @self.point_of_reference")[gene].notna()):
                try:
                    self.df_FDR[
                        phos_pos
                    ] = 0.0  # tiny P values due to no CTRL data need to be zeroed
                except IndexError:
                    print(
                        "phos_pos, pos_df", phos_pos, pos_df, df_melted, sep="\n\n"
                    )  # TODO - this is a hack, need to figure out why this happens. maybe cuz 846 has no data?

        return PhosphoProtData(
            name=self.name,
            df=df_melted,
            df_FDR=self.df_FDR[phos_sites + ["test"]],
            gene=gene,
            ordered_test_names=[test for test in tests if test != "ID"],
            stats=stats,
        )

    def get_FDR(self, test, prot_pos=None):
        try:
            return self.df_FDR.query("test == @test")[prot_pos].values[0]
        except KeyError:
            return 0.0

    def get_median_CPMs(self, test, prot_pos=None):
        return np.median(self.plot_df.query("gene == @prot_pos").fillna(0.0)[self.gene])

    @property
    def point_of_reference(self):
        return "CTRL"

    @property
    def plot_df(self):
        return self.df

    @property
    def test_names(self):
        return sorted(
            set([test for test in list(self.df["test"]) if test != Schema.UNQIUE_ID])
        )


@dataclass(slots=True, frozen=True)
class RNASeqData(Data2):
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

    # name: str
    df: pl.DataFrame = field(repr=False)
    df_FDR: pd.DataFrame = field(repr=False)
    processed_dfs: dict[str, pd.DataFrame] = field(repr=False)
    # stats: bool = True
    # gene: str | None = None

    def filter(self, gene: str, tests: list[str], stats: dict[str, bool]): #tests: list[str] = None,
        df: pl.DataFrame = self.df.filter(self.df["test"].is_in(tests))
        df = df.select([gene, "test", "point_of_ref"])
        return RNASeqData(
            name=self.name,
            df=df,
            df_FDR=self.df_FDR[gene],
            processed_dfs=self.processed_dfs,
            gene=gene,
            ordered_test_names=[test for test in tests if test != "ID"],
            stats=stats,
        )

    # TODO - post init to check that indexes are equal. found examples of missign and duplicates gene names

    @property
    def plot_df(self):
        return self.df.to_pandas()

    @property
    def test_names(self):
        return sorted(set(self.df["test"]))

    @property
    def degs(self):
        return set(self.processed_dfs.keys())

    @property
    def point_of_reference(self):
        """frgrg"""

        ref_df = self.df.filter(pl.col("point_of_ref") == "yes")
        point_of_ref = set(ref_df["test"])
        if len(point_of_ref) == 1:
            return point_of_ref.pop()
        elif len(point_of_ref) > 1:
            raise ValueError("Multiple comparisons have POR tag!")
        else:
            return None


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
        # df = pl.read_csv(csv, separator=delimitor).to_pandas(
        #     use_pyarrow_extension_array=True
        # )
        df = pd.read_csv(csv, sep=delimitor)
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


def read_CPM(path: Path) -> pl.DataFrame:
    fin = list(path.glob("*.csv"))
    if fin:
        fin = fin.pop()
        # df = pl.read_csv(str(fin)).to_pandas(use_pyarrow_extension_array=True)
        df = pd.read_csv(fin)  # TODO como
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
    FDR: dict[str, pd.DataFrame] = load_processed_rna_files(path)
    merged_FDRs = merge_FDR(FDR)
    CPM = pl.from_pandas(CPM)

    return RNASeqData(
        name=path.stem,
        df=CPM,
        df_FDR=merged_FDRs,
        processed_dfs=FDR,
        gene="NA",
        ordered_test_names=["NA"],
        stats={'log': False, 'brackets': True},
    )


class Schema:
    UNQIUE_ID = "ID"
    POS_IN_PROT = "Positions in Proteins"
    MOD_IN_PROT = "Modifications in Proteins"
    COMBINED_PROT_DESCRIPTION = (
        "Master Protein Descriptions"  # was just  "Description" in prot
    )
    GENE = "gene"
    CAT = "test"


Preprocessor = Callable[[pd.DataFrame], pd.DataFrame]


def compose(*functions: Preprocessor) -> Preprocessor:
    """Helper function to call all df functions sequencially"""
    return reduce(lambda func1, func2: lambda x: func2(func1(x)), functions)


def description_to_str(df: pd.DataFrame) -> pd.DataFrame:
    df[Schema.COMBINED_PROT_DESCRIPTION] = df[Schema.COMBINED_PROT_DESCRIPTION].astype(
        "str"
    )

    return df


def make_gene_col(df: pd.DataFrame) -> pd.DataFrame:
    df[Schema.GENE] = df[Schema.COMBINED_PROT_DESCRIPTION].apply(
        lambda x: [val for val in x.strip().split(" ") if "GN=" in val]
        .pop()
        .replace("GN=", "")
        if "GN=" in x
        else x.split(" ")[-1]
    )
    return df


def make_gene_col_unique(df: pd.DataFrame, col: str) -> pd.DataFrame:
    non_unique_genes = list(df[col])
    counts = Counter(non_unique_genes)
    unique_genes = []
    for i, gene in enumerate(non_unique_genes):
        if counts[gene] > 1:
            unique_genes.append(
                f"{gene}_{counts[gene]-non_unique_genes[:i].count(gene)}"
            )
        else:
            unique_genes.append(gene)

    df[col] = unique_genes
    return df


def remove_unwanted_cols(df: pd.DataFrame, misc: list[str]) -> pd.DataFrame:
    abundances: list[str] = [
        col
        for col in list(df.columns)
        if col.startswith(("Abundance: F", "Abundance Ratio Adj. P-Value:"))
    ]
    df = df[misc + abundances]

    return df


def add_categories(df: pd.DataFrame) -> pd.DataFrame:
    df = df.set_index(Schema.GENE)
    df = df.T
    index_str = "index_str"
    df[index_str] = df.index
    df[Schema.CAT] = df[index_str].apply(
        lambda x: "P_value" if "P-Value:" in x else str(x).split(",")[-1].strip()
    )
    df[Schema.CAT] = df[Schema.CAT].astype("category")
    del df[index_str]

    return df.T


def read_excel(path: Path) -> pd.DataFrame:
    # load the data
    files = [fin for fin in list(path.glob("*xlsx")) if not rubbish(fin.name)]
    if len(files) != 1:
        raise ValueError(
            f"The len of input should be one but it is {len(files)}; {files}"
        )
    return pd.read_excel(files[0])
    # pls caused too many dramas ie The current offset in the file is 3909947 bytes.
    # can't handle values like 0.001000


def split_dfs(df):
    df = df.T
    df_FDR = df.query('test == "P_value"')
    if Schema.UNQIUE_ID in df.index:  # this changes test to ID, doe sit matter?
        df_FDR.columns = df.iloc[0]
        df_FDR.drop(df_FDR.index[0])
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
        description_to_str,
        make_gene_col,
        partial(make_gene_col_unique, col=Schema.GENE),
        partial(remove_unwanted_cols, misc=[Schema.GENE]),
        add_categories,
    ]
    preprocessor = compose(*pipe)
    df = preprocessor(df)
    df, df_FDR = split_dfs(df)
    # df_FDR.to_csv(f"~/Downloads/{path.stem}_df_FDR.csv")
    return ProtData(
        name=path.stem,
        df=pl.from_pandas(df),
        df_FDR=df_FDR,
        gene="NA",
        ordered_test_names=["NA"],
        stats={'log': False, 'brackets': True},
    )


def create_unqiue_id(df: pd.DataFrame) -> pd.DataFrame:
    def choose_pos(gene: str, pos_extact: str | float, pos_range: str) -> str:
        if not isinstance(pos_extact, float):
            poses = pos_extact.strip().split(" ")[-1]
            pos: str = "_".join(bit.split("(")[0][1:] for bit in poses.split(";"))
        else:
            pos: str = pos_range.split(" ")[-1][1:-1]
        return "_".join([gene, pos])

    df[Schema.UNQIUE_ID] = df.apply(
        lambda x: choose_pos(
            x[Schema.GENE],
            x[Schema.MOD_IN_PROT],
            x[Schema.POS_IN_PROT],
        ),
        axis=1,
    )
    return df


def load_phospho_data(path: Path) -> PhosphoProtData:
    df = read_excel(path)
    pipe = [
        description_to_str,
        make_gene_col,
        create_unqiue_id,
        partial(make_gene_col_unique, col=Schema.UNQIUE_ID),
        partial(remove_unwanted_cols, misc=[Schema.GENE, Schema.UNQIUE_ID]),
        add_categories,
    ]
    preprocessor = compose(*pipe)
    df = preprocessor(df)
    df, df_FDR = split_dfs(df)
    # df_FDR.to_csv(f"~/Downloads/{path.stem}_df_FDR.csv")
    # df.to_csv(f"~/Downloads/{path.stem}_df.csv")
    return PhosphoProtData(
        name=path.stem,
        df=df,
        df_FDR=df_FDR,
        gene="NA",
        ordered_test_names=["NA"],
        stats={'log': True, 'brackets': True},
    )


class SchemaFunction:
    ARRHYTHMIA = "Arrhythmia"
    CONDITION = "Condition"
    NORM1 = "Normalisation 1"
    NORM2 = "Normalisation 2"
    KEY = "key"
    DRUG = "Drug"
    TIME = "Timepoint"
    DOSE = "Dose"
    WELL = "Well Number"
    FORCE = "Force"
    RATE = "Rate"
    TA50 = "Ta50"
    TR50 = "Tr50"
    TPEAK = "Tpeak 85"
    TA_15_30 = "Ta 15-30"
    TA_30_85 = "Ta 30-85"
    TR_85_50 = "Tr 85-50"
    TR_50_15 = "Tr 50-15"
    RR_SCAT = "RRscat"
    DATASET = "dataset"


@dataclass(slots=True, frozen=True)
class FunctionData:
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
    ec50: pd.DataFrame = field(repr=False)
    names: list[str] | None = None
    metric: str | None = None
    drug: str | None = None

    def filter(
        self, test: str, conditions: list[str], datasets: list[str], metric: str
    ):  # test is drug, here
        filt = "Drug == @test & Condition == @conditions & dataset == @datasets"
        df = self.df.query(filt).copy(deep=True)
        return FunctionData(
            names=datasets, df=df, ec50=self.ec50, metric=metric, drug=test
        )

    def _create_pivot_df(self, df: pd.DataFrame) -> pd.DataFrame:
        cols = [self.metric, SchemaFunction.TIME, SchemaFunction.WELL]
        df: pd.DataFrame = df[cols]
        df = df.reset_index()[cols]
        df = df.pivot(
            index=SchemaFunction.TIME, columns=SchemaFunction.WELL, values=self.metric
        ).T
        df = df[["Baseline", "1", "2", "3", "4", "5"]]
        df.columns = [0, 1, 2, 3, 4, 5]

        return df

    def _arrhythmias(self, df: pd.DataFrame) -> pd.DataFrame:
        return (
            df[[SchemaFunction.DRUG, SchemaFunction.DOSE, SchemaFunction.ARRHYTHMIA]]
            .groupby(SchemaFunction.DRUG)
            .value_counts()
            .to_frame()
            .reset_index()
        )

    def _get_mean_values(
        self, df: pd.DataFrame, dose: dict[str, float]
    ) -> tuple[list[str], list[str], list[str]]:
        datas = df.describe().to_dict()
        xs = [dose[str(time_point)] for time_point in list(datas.keys())]
        ys = [v.get("50%") for v in datas.values()]
        stds = [v.get("std") for v in datas.values()]

        return xs, ys, stds

    def _get_arrhythmia_counts(self, df, xs):
        arrhythmia = self._arrhythmias(df)
        arrhythmia_counts: list[str] = []
        for x in xs:
            arrhythmia_at_dose = arrhythmia.query("Dose == @x")
            dd: dict[str, int] = dict(
                zip(arrhythmia_at_dose["Arrhythmia"], arrhythmia_at_dose["count"])
            )
            total = sum(dd.values())
            non_arrhythmic: int = total - dd.get("Y", 0)
            arrhythmia_counts.append(f"{str(non_arrhythmic)}/{str(total)}")

        return arrhythmia_counts

    def possible_dataset_names(self, drug: str) -> list[str]:
        return list(set(self.df.query("Drug == @drug")[SchemaFunction.DATASET]))

    def possible_condition_names(self, drug: str) -> list[str]:
        return list(set(self.df.query("Drug == @drug")[SchemaFunction.CONDITION]))

    @property
    def plot_df(self):
        i = 0
        tmp_d = defaultdict(list)

        for name, df in self.df.groupby(SchemaFunction.DATASET):
            i += 1
            mapping = get_mappings(
                df,
                index_cols=[SchemaFunction.TIME],
                value_cols=[
                    SchemaFunction.DOSE,
                ],
            )
            condition = set(df[SchemaFunction.CONDITION])
            assert len(condition) == 1
            pivot_df = self._create_pivot_df(df)

            dose = mapping[SchemaFunction.DOSE]
            dose["0"] = 0

            xs, ys, stds = self._get_mean_values(pivot_df, dose)
            arrhythmia_counts = self._get_arrhythmia_counts(df, xs)

            tmp_d[self.drug] += xs
            tmp_d[self.metric] += ys
            tmp_d["stds"] += stds
            tmp_d["arrhythmia_counts"] += arrhythmia_counts
            tmp_d["name & condition"] += [f"{name}, {condition.pop()}"] * len(xs)
        df = pd.DataFrame.from_dict(tmp_d)
        return df

    @property
    def experiments(self):
        return [name.strip().split("_")[1] for name in self.names]

    @property
    def test_names(self):
        return sorted(set(self.df[SchemaFunction.DRUG]))

    @property
    def metrics(self):
        return [
            name for Schema, name in vars(SchemaFunction).items() if "__" not in Schema
        ][9:-1]


# Preprocessor = Callable[[pd.DataFrame], pd.DataFrame]


# def compose_functiona_data(*functions: Preprocessor) -> Preprocessor:
#     """Helper function to call all df functions sequencially"""
#     return reduce(lambda func1, func2: lambda x: func2(func1(x)), functions)


def remove_leading_0(df: pd.DataFrame) -> pd.DataFrame:
    df[SchemaFunction.WELL] = df[SchemaFunction.WELL].apply(
        lambda x: x.replace("Pos0", "Pos") if "Pos0" in x else x
    )
    return df


def time_str(df: pd.DataFrame) -> pd.DataFrame:
    df[SchemaFunction.TIME] = df[SchemaFunction.TIME].astype(str)
    return df


def remove_arrhythmia(df: pd.DataFrame) -> pd.DataFrame:
    if SchemaFunction.ARRHYTHMIA in df.columns:
        df = df.query('Arrhythmia == "N"')
    return df


def key(df: pd.DataFrame) -> pd.DataFrame:
    df[SchemaFunction.KEY] = df[SchemaFunction.TIME] + ":" + df[SchemaFunction.WELL]
    return df


def del_uneeded(df: pd.DataFrame) -> pd.DataFrame:
    del df["Plate"]
    return df


def index_key(df: pd.DataFrame) -> pd.DataFrame:
    return df.set_index("key")


def add_dataset_name(df: pd.DataFrame, name: str) -> pd.DataFrame:
    df[SchemaFunction.DATASET] = name
    return df


def put_required_cols_back(df: pd.DataFrame, original_df: pd.DataFrame) -> pd.DataFrame:
    original_df = original_df[[col for col in original_df if col not in df.columns]]
    df = df.join(original_df, how="inner")
    return df


def index_copy(df: pd.DataFrame) -> pd.DataFrame:
    df["index_copy"] = df.index
    return df


def index_copy2(df: pd.DataFrame) -> pd.DataFrame:
    df["index_copy"] = df["index_copy"].apply(lambda x: x[1])
    return df


def zeros(df: pd.DataFrame) -> pd.DataFrame:
    df = df.fillna(0.0)
    return df


def time_and_pos(df: pd.DataFrame) -> pd.DataFrame:
    df[[SchemaFunction.TIME, SchemaFunction.WELL]] = df["index_copy"].str.split(
        ":", expand=True
    )
    return df


def get_mappings(df, index_cols, value_cols):
    df = df[index_cols + value_cols]
    return df.set_index(index_cols).to_dict()


def find_arrhythmias(df):
    return (
        df[[SchemaFunction.DRUG, SchemaFunction.DOSE, SchemaFunction.ARRHYTHMIA]]
        .groupby("Drug")
        .value_counts()
        .to_frame()
        .reset_index()
    )


def metric_names():
    return [
        name for Schema, name in vars(SchemaFunction).items() if "__" not in Schema
    ][9:-2]


def normalize_group(group):
    drug_data_dict = group.to_dict()
    norm = {}
    for data in metric_names():
        norm[data] = {}
        for k, v in drug_data_dict[data].items():
            _, pos = k.split(":")
            baseline_key = "Baseline:" + pos
            norm[data][k] = v / drug_data_dict[data][baseline_key]
    return pd.DataFrame(norm)


def replace_time(df):
    df[SchemaFunction.TIME] = df.index
    df[SchemaFunction.TIME] = df[SchemaFunction.TIME].apply(lambda x: x.split(":")[0])
    return df


def normalize_group2(group, DMSO):
    group = replace_time(group)
    DMSO = replace_time(DMSO)
    drug_data_dict = group.to_dict()
    norm = {}
    for data in metric_names():
        norm[data] = {}
        time_point_1 = ""
        for k, v in drug_data_dict[data].items():
            time_point_2, _ = k.split(":")
            if time_point_1 != time_point_2:
                mean = DMSO.loc[DMSO["Timepoint"] == time_point_2, data].mean()
                time_point_1 = time_point_2
            norm[data][k] = v / mean

    return pd.DataFrame(norm)


def load_function_data(path: Path) -> FunctionData:
    files = [fin for fin in list(path.glob("*xlsx")) if not rubbish(fin.name)]
    dfs = []
    for fin in files:
        if "Experiment_" not in fin.stem:
            continue
        df = pd.read_excel(fin)

        pipe = [
            remove_leading_0,
            time_str,
            remove_arrhythmia,
            key,
            del_uneeded,
            index_key,
        ]
        preprocessor = compose(*pipe)
        df = preprocessor(df)
        normalized_df = df.groupby("Drug", group_keys=True).apply(normalize_group)
        normalized_df = normalized_df.reset_index(level=["Drug"])
        # normalized_df.to_csv("~/Downloads/normalized_dfgg.csv")
        sub_dfs = []
        for name, drug_df in df.groupby(SchemaFunction.DRUG):
            mapping = get_mappings(
                drug_df,
                index_cols=[SchemaFunction.TIME],
                value_cols=[
                    SchemaFunction.NORM2,
                    SchemaFunction.NORM1,
                ],
            )
            assert mapping[SchemaFunction.NORM1]["Baseline"] == "Baseline"
            norm2 = mapping[SchemaFunction.NORM2]["Baseline"]
            DMSO = normalized_df.query("Drug == @norm2")
            normalized_df2 = normalize_group2(
                normalized_df.query("Drug == @name"), DMSO
            )
            pipe = [
                index_copy,
                zeros,
                time_and_pos,
                partial(add_dataset_name, name=fin.stem),
                partial(put_required_cols_back, original_df=drug_df),
            ]
            preprocessor = compose(*pipe)
            normalized_df2 = preprocessor(normalized_df2)
            sub_dfs.append(normalized_df2)
        dfs.append(pd.concat(sub_dfs))
    df = pd.concat(dfs)
    # arrhythmia = pd.concat(arrhythmias)

    df_EC50 = pd.read_excel([fin for fin in files if "EC50" in fin.stem][0])
    return FunctionData(df=df, ec50=df_EC50)
