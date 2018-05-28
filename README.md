
# Overview

This repo consists of a walkthrough for creating a dataset of protein sequences in a similar fashion as was done by [Zhou & Troyanskaya, 2014](https://arxiv.org/pdf/1403.1347.pdf). The authors of that paper used a dataset constructed as outlined below to train convolutional generative stochastic networks for protein secondary structure prediction (PSSP). The modest goal is to introduce individuals with experience in machine/deep learning to the basics of structural bioinformatics by constructing a dataset for PSSP as well as other, related tasks.

### Packages Used

- Python 3
- Pandas
- BioPython
- TensorFlow

## Proteins - primary sequence and secondary structure

Proteins are one of the most important classes of macromolecules in biology. They are composed of linear chains of 20 (sometimes 22) proteinogenic amino acids. Each amino acid can be represented by a single letter, allowing us to express proteins as a string. For example, the sequence for hemoglobin is:
>VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR

When found in the appropriate biological context, proteins assume - through a complex interplay of physicochemical forces - a 3-dimensional conformation, referred to as the protein's *tertiary* or *native* structure. This native state determines how a particular protein behaves; that is, form determines function. Accurately predicting tertiary structure is therefore the first and most fundamental step for understanding the function of discovered proteins as well as for designing new proteins with entirely novel functions.

During the folding process, regions along the protein backbone assume local conformations called secondary structures. These secondary structures are characterized by the patterns of hydrogen bonds formed between amino acids in a particular region. One common way of categorizing such structures is given by the dictionary of protein secondary structires [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/). Since nearly all secondary structures are localized to a particular region of the backbone, they can also be written as a string of letters, where each letter corresponds to a particular secondary structure. Hemoglobin, aligned with its secondary structures, is as follows:

>       VLSPADKTNV KAAWGKVGAH AGEYGAEALE RMFLSFPTTK TYFPHFDLSH GSAQVKGHGK KVADALTNAV AHVDDMPNAL SALSDLHAHK LRVDPVNFKL LSHCLLVTLA AHLPAEFTPA VHASLDKFLA SVSTVLTSKY R
>          HHHHHHH HHHHHHHGGG HHHHHHHHHH HHHHH GGGG GG TTS  ST T HHHHHHHH HHHHHHHHHH HTTTSHHHHT HHHHHHHHHT T   THHHHH HHHHHHHHHH HH TTT  HH HHHHHHHHHH HHHHHHHTT
           

## Combining Cull PDB and DSSP

The protein data bank [(PDB)](https://www.rcsb.org/) is an online repository of hundreds of thousands of proteins whose 3d structure has been resolved experimentally. We will use the web server CullPDB: PISCES to extract a diverse set of proteins from the PDB, and the use DSSP to find their secondary structures.

### Cull PDB: PISCES

Proteins from the PDB can be queried based on criteria such as resolution, sequence identity, etc. It's possible (as of 20/03/2018) to download different lists [here](http://dunbrack.fccc.edu/Guoli/pisces_download.php).

### DSSP files

The PISCES lists provide PDB ID's, but they do not have the secondary structure information. To get that, we need to download the DSSP information from the PDB. This can be done directly [here](http://swift.cmbi.ru.nl/gv/dssp/).

By syncing the database locally, the individual \*.dssp files can be parsed by the script [parse_dssp.py](./parse_dssp.py).

### Note on CPDB chains and DSSP ids
The PISCES server checks each CHAIN of a PDB entry individually. As such, the cpdb IDs may contain all or only some of the chains of a particular PDB entry. On the other hand, the DSSP outputs a single file / entry per PDB ID, which will include all of the chains for that entry.

The `parse_dssp.py` script will split dssp files into their constituent chains, appending the chain id to the end of the dssp_id, and creating records for each chain. Some edges cases where the parser incorrectly identifies the chain id exist, but those are skipped.

### Next Steps
The rest of this notebook will assume that a list downloaded from PISCES as well as some number of .csv files containing the parsed DSSP data exist in a `data/dssp` folder.

## Joining the Data

We want to do a join on the PDB id field of the PISCES and DSSP datasets. Since these are both in either tab-separated or csv format, Pandas is an ideal candidate for doing this.

We want to concatenate the two datasets, joining on the two id's. Since the cpdb data is a subset of the dssp data, we join on the cpdb id field. See [merge_data.py](./merge_data.py) for more details.

In the code below, we drop sequences below 26 residues in length and save the results of the merge to the file `cpdb2_records.csv`.


```python
import pandas as pd
from pathlib import Path
datadir = str(Path(Path.home(), "data", "dssp"))
```


```python
from merge_data import merge
merged = merge()

merged = merged[merged.seq.str.len() > 25]
merged.to_csv(datadir+"/cpdb2_records.csv")
```

## Generate the Position-Specific Similary Matrices
Similar to CPDB, we can calculate position-specific profile similarity scores using PSI-BLAST. The process is as follows:

### CPDB2 to FASTA format
In order to make calculating the PSSMs amenable to the CPDB2 dataset, we save the `dssp_id` and `seq` fields in FASTA file format.
Some of these records contain leftover `!` gap symbols, so these are replaced with `*` to indicate a gap of indeterminate length
in the sequence, as described [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp)

See the [csv_to_fasta.py](./csv_to_fasta.py) file.

### Download and preprocess UniRef90
Download the [UniRef90](http://www.uniprot.org/downloads) dataset, and filter using [pfilt](http://bioinf.cs.ucl.ac.uk/psipred/) to remove low information content and coiled-coil regions.

#### pfilt
See the [README](http://bioinfadmin.cs.ucl.ac.uk/downloads/pfilt/).


### Download BLAST+ and create a local database
Create a BLAST database out of the filtered sequences in FASTA format using the blast command line tool, described [here](https://www.ncbi.nlm.nih.gov/books/NBK279688/). The BLAST+ software is available [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

### Run PSIBLAST to get the scores
Run PSI-BLAST on the CullPDB dataset downloaded via PISCES against the BLAST database just created, 
with an inclusion threshold of 0.001 for 3 iterations.

Transform the profile scores into the 
range [0, 1) via logistic sigmoid
```bash
psiblast -db uniref90_filt_db -query example.fasta -out output.psiblast -evalue 0.001 -num_threads 4 -num_iterations 3 -out_pssm output.pssm -out_ascii_pssm asci_test.pssm -save_pssm_after_last_round
```

**Note** that PSIBLAST will only save the PSSM for a SINGLE protein at a time. Each new profile score overwrites the previous, so the sequences need to be run independently.

# Creating a dataset without PSIBLAST

## Creating .TFRecords files

We could use the string pairs directly as inputs to a learning model. The TensorFlow tf.data API allows for reading text data and converting to feature vectors as a preprocessing step in a model. However, this places a heavy computational bottleneck at the CPU that could slow down training. Since the dataset is ~14k sequences, we can save them directly as feature vectors in a .tfrecords file that is loaded into memory at training.

## Features

Oftentimes, position-specific features are calculated using PSIBLAST or hidden markov models. To keep things simple at this stage, we can simply append feature vectors that correspond to the amino acids / secondary structures of the sequence and save those as TF records. The features are the following:


```python
aa_feats = pd.read_csv("./cpdb2_aa_features.csv", index_col=0)
aa_feats
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>A</th>
      <th>C</th>
      <th>D</th>
      <th>E</th>
      <th>F</th>
      <th>G</th>
      <th>H</th>
      <th>I</th>
      <th>K</th>
      <th>L</th>
      <th>...</th>
      <th>X</th>
      <th>!</th>
      <th>SOS</th>
      <th>EOS</th>
      <th>hydrophobicity</th>
      <th>polar</th>
      <th>hydropathy intensity</th>
      <th>hydrophilicity</th>
      <th>pH_l</th>
      <th>vdW_vol</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>A</th>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.8</td>
      <td>3.0</td>
      <td>6.01</td>
      <td>67.0</td>
    </tr>
    <tr>
      <th>C</th>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>-1.0</td>
      <td>2.5</td>
      <td>-1.0</td>
      <td>5.05</td>
      <td>86.0</td>
    </tr>
    <tr>
      <th>D</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>-3.5</td>
      <td>3.0</td>
      <td>2.85</td>
      <td>91.0</td>
    </tr>
    <tr>
      <th>E</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>-3.5</td>
      <td>3.0</td>
      <td>3.15</td>
      <td>109.0</td>
    </tr>
    <tr>
      <th>F</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>-1.0</td>
      <td>2.8</td>
      <td>-2.5</td>
      <td>5.49</td>
      <td>135.0</td>
    </tr>
    <tr>
      <th>G</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-0.4</td>
      <td>0.0</td>
      <td>6.06</td>
      <td>48.0</td>
    </tr>
    <tr>
      <th>H</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>1.0</td>
      <td>-3.2</td>
      <td>-0.5</td>
      <td>7.60</td>
      <td>118.0</td>
    </tr>
    <tr>
      <th>I</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>-1.0</td>
      <td>4.5</td>
      <td>-1.8</td>
      <td>6.05</td>
      <td>124.0</td>
    </tr>
    <tr>
      <th>K</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>-3.9</td>
      <td>3.0</td>
      <td>9.60</td>
      <td>135.0</td>
    </tr>
    <tr>
      <th>L</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>-1.0</td>
      <td>3.8</td>
      <td>-1.8</td>
      <td>6.01</td>
      <td>124.0</td>
    </tr>
    <tr>
      <th>M</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>-1.0</td>
      <td>1.9</td>
      <td>-1.3</td>
      <td>5.74</td>
      <td>124.0</td>
    </tr>
    <tr>
      <th>N</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-2.0</td>
      <td>1.0</td>
      <td>-3.5</td>
      <td>0.2</td>
      <td>5.41</td>
      <td>96.0</td>
    </tr>
    <tr>
      <th>P</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-1.0</td>
      <td>0.0</td>
      <td>1.6</td>
      <td>0.0</td>
      <td>6.30</td>
      <td>90.0</td>
    </tr>
    <tr>
      <th>Q</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-2.0</td>
      <td>1.0</td>
      <td>-3.5</td>
      <td>0.2</td>
      <td>5.65</td>
      <td>114.0</td>
    </tr>
    <tr>
      <th>R</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>-4.5</td>
      <td>3.0</td>
      <td>10.76</td>
      <td>148.0</td>
    </tr>
    <tr>
      <th>S</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-2.0</td>
      <td>0.0</td>
      <td>-0.8</td>
      <td>0.3</td>
      <td>5.68</td>
      <td>73.0</td>
    </tr>
    <tr>
      <th>T</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-2.0</td>
      <td>0.0</td>
      <td>-0.7</td>
      <td>-0.4</td>
      <td>5.60</td>
      <td>93.0</td>
    </tr>
    <tr>
      <th>V</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>-1.0</td>
      <td>4.2</td>
      <td>-1.5</td>
      <td>6.00</td>
      <td>105.0</td>
    </tr>
    <tr>
      <th>W</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>-1.0</td>
      <td>-0.9</td>
      <td>-3.4</td>
      <td>5.89</td>
      <td>163.0</td>
    </tr>
    <tr>
      <th>Y</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-2.0</td>
      <td>-1.0</td>
      <td>-1.3</td>
      <td>-2.3</td>
      <td>5.64</td>
      <td>141.0</td>
    </tr>
    <tr>
      <th>X</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>!</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>SOS</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>EOS</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>24 rows × 30 columns</p>
</div>




```python
ss_feats = pd.read_csv("./cpdb2_ss_features.csv", index_col=0)
ss_feats
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>H</th>
      <th>B</th>
      <th>E</th>
      <th>G</th>
      <th>I</th>
      <th>T</th>
      <th>S</th>
      <th>U</th>
      <th>SOS</th>
      <th>EOS</th>
    </tr>
    <tr>
      <th>labels</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>H</th>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>B</th>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>E</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>G</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>I</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>T</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>S</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>U</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>SOS</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>EOS</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
    </tr>
  </tbody>
</table>
</div>



The secondary structures are the targets, and so their "features" are simply one-hot encodings of the letters, including special SOS/EOS tokens.
The amino acid features include a one-hot encoding the 20 proteinogenic amino acids as well as six physicochemical properties: hydrophobicity, polarity, hydropathy intensity, hydrophilicity, isoelectric point, and van der Waals volume. 

See the [make_tfrecords.py](./make_tfrecords.py) script for how the protein sequences are processed into tf records.

These files are read into TensorFlow models using the tf.data API.

### Note on cpdb2 sequences with the character 'b', and 'j', and 'o', and 'u', and 'z'
According to the description of DSSP, lowercase characters indicate a SS-bridge Cysteine. These come in pairs; only some of them show up as bad characters, because when capitalized, some are valid amino acid codes.

The strategy for these is thus to replace them with 'C' for cysteine.
