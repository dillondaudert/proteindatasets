# convert sequences and ids to fasta format
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio import SeqIO

# load the data
data = pd.read_csv("/home/dillon/data/dssp/cpdb2_dssp_14726.csv")
data["seq"] = data["seq"].str.replace("!", "*", regex=False)

seqrecs = [SeqRecord(Seq(data.loc[prot]["seq"], generic_protein), id=data.loc[prot]["dssp_id"]) \
           for prot in range(data.shape[0])]

SeqIO.write(seqrecs, "cpdb2_dssp_14726.fasta", "fasta")