
import pandas as pd
from pathlib import Path
datadir = str(Path(Path.home(), "data", "dssp"))


cpdb_df = pd.read_csv(datadir+"/cullpdb_pc30_res2.5_R1.0_d180208_chains15102.txt", delim_whitespace=True)

# get a dataframe consisting of only the first 4 characters of each cpdb id
cpdb_ids = cpdb_df["IDs"].str[0:4].str.lower()
# only take the unique entries
cpdb_ids = cpdb_ids[~cpdb_ids.duplicated()]

# get the cpdb dataset with only unique, 4-letter IDs
cpdb_df_unique = cpdb_df[~cpdb_df["IDs"].duplicated()]

# change the name of the IDs column to dssp_id
cpdb_df_unique["dssp_id"] = cpdb_ids
cpdb_df_unique.drop(labels=["IDs"], axis=1)

merged_frames = []
# now, for each dssp file, take an inner join with the unique cpdb IDs
for i in range(12):
    dssp_df = pd.read_csv(datadir+"/dssp_%d.csv" % i)
    merged = dssp_df.merge(cpdb_df_unique, how="inner", on="dssp_id")
    print(len(merged))
    merged_frames.append(merged)

all_merged = pd.concat(merged_frames)
# Write out the sequence and the secondary structure information to a file
all_merged[["dssp_id", "seq", "ss"]].to_csv(datadir+"/cpdb_dssp_%d.csv" % len(all_merged), index=False)