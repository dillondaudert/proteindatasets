
import pandas as pd
from pathlib import Path
datadir = str(Path(Path.home(), "data", "dssp"))

def merge():

    cpdb_df = pd.read_csv(datadir+"/cullpdb_pc30_res2.5_R1.0_d180208_chains15102.txt", delim_whitespace=True)

    # change the name of the IDs column to dssp_id
    cpdb_df["dssp_id"] = cpdb_df["IDs"].str.lower()
    cpdb_df.drop(labels=["IDs"], axis=1)

    dssp_df = pd.read_csv(datadir+"/dssp_records.csv")
    dssp_df["dssp_id"] = dssp_df["dssp_id"].str.lower()
    merged = dssp_df.merge(cpdb_df, how="inner", on="dssp_id")

    return merged