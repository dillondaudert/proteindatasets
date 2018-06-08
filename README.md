
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

### Note on CPDB sequences with the character 'b', and 'j', and 'o', and 'u', and 'z'
According to the description of DSSP, lowercase characters indicate a SS-bridge Cysteine. These come in pairs; only some of them show up as bad characters, because when capitalized, some are valid amino acid codes.

The strategy for these is thus to replace them with 'C' for cysteine.

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
# from merge_data import merge
# merged = merge()

# merged = merged[merged.seq.str.len() > 25]
# merged.to_csv(datadir+"/cpdb2_records.csv")
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

# Physicochemical properties

Along with position-specific similarity matrix scores, we also use a number of physicochemical properties to describe amino acids:


```python
aa_feats = pd.read_csv("aa_feats.csv", index_col=0)
aa_feats = aa_feats.drop(labels=["X"], axis=0)
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
      <th>hydrophobicity</th>
      <th>polar</th>
      <th>hydropathy intensity</th>
      <th>hydrophilicity</th>
      <th>pH_l</th>
      <th>vdW_vol</th>
      <th>pK1</th>
      <th>pK2</th>
      <th>steric</th>
      <th>polarizability</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>A</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.8</td>
      <td>3.0</td>
      <td>6.01</td>
      <td>67.0</td>
      <td>2.35</td>
      <td>9.87</td>
      <td>1.28</td>
      <td>0.05</td>
    </tr>
    <tr>
      <th>C</th>
      <td>1.0</td>
      <td>-1.0</td>
      <td>2.5</td>
      <td>-1.0</td>
      <td>5.05</td>
      <td>86.0</td>
      <td>1.92</td>
      <td>10.70</td>
      <td>1.77</td>
      <td>0.13</td>
    </tr>
    <tr>
      <th>D</th>
      <td>2.0</td>
      <td>1.0</td>
      <td>-3.5</td>
      <td>3.0</td>
      <td>2.85</td>
      <td>91.0</td>
      <td>1.99</td>
      <td>9.90</td>
      <td>1.60</td>
      <td>0.11</td>
    </tr>
    <tr>
      <th>E</th>
      <td>2.0</td>
      <td>1.0</td>
      <td>-3.5</td>
      <td>3.0</td>
      <td>3.15</td>
      <td>109.0</td>
      <td>2.10</td>
      <td>9.47</td>
      <td>0.00</td>
      <td>0.15</td>
    </tr>
    <tr>
      <th>F</th>
      <td>1.0</td>
      <td>-1.0</td>
      <td>2.8</td>
      <td>-2.5</td>
      <td>5.49</td>
      <td>135.0</td>
      <td>2.20</td>
      <td>9.31</td>
      <td>2.94</td>
      <td>0.29</td>
    </tr>
    <tr>
      <th>G</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>-0.4</td>
      <td>0.0</td>
      <td>6.06</td>
      <td>48.0</td>
      <td>2.35</td>
      <td>9.78</td>
      <td>0.00</td>
      <td>0.00</td>
    </tr>
    <tr>
      <th>H</th>
      <td>-1.0</td>
      <td>1.0</td>
      <td>-3.2</td>
      <td>-0.5</td>
      <td>7.60</td>
      <td>118.0</td>
      <td>1.80</td>
      <td>9.33</td>
      <td>2.99</td>
      <td>0.23</td>
    </tr>
    <tr>
      <th>I</th>
      <td>1.0</td>
      <td>-1.0</td>
      <td>4.5</td>
      <td>-1.8</td>
      <td>6.05</td>
      <td>124.0</td>
      <td>2.32</td>
      <td>9.76</td>
      <td>4.19</td>
      <td>0.19</td>
    </tr>
    <tr>
      <th>K</th>
      <td>2.0</td>
      <td>1.0</td>
      <td>-3.9</td>
      <td>3.0</td>
      <td>9.60</td>
      <td>135.0</td>
      <td>2.16</td>
      <td>9.06</td>
      <td>1.89</td>
      <td>0.22</td>
    </tr>
    <tr>
      <th>L</th>
      <td>1.0</td>
      <td>-1.0</td>
      <td>3.8</td>
      <td>-1.8</td>
      <td>6.01</td>
      <td>124.0</td>
      <td>2.33</td>
      <td>9.74</td>
      <td>2.59</td>
      <td>0.19</td>
    </tr>
    <tr>
      <th>M</th>
      <td>1.0</td>
      <td>-1.0</td>
      <td>1.9</td>
      <td>-1.3</td>
      <td>5.74</td>
      <td>124.0</td>
      <td>2.13</td>
      <td>9.28</td>
      <td>2.35</td>
      <td>0.22</td>
    </tr>
    <tr>
      <th>N</th>
      <td>-2.0</td>
      <td>1.0</td>
      <td>-3.5</td>
      <td>0.2</td>
      <td>5.41</td>
      <td>96.0</td>
      <td>2.14</td>
      <td>8.72</td>
      <td>1.60</td>
      <td>0.13</td>
    </tr>
    <tr>
      <th>P</th>
      <td>-1.0</td>
      <td>0.0</td>
      <td>1.6</td>
      <td>0.0</td>
      <td>6.30</td>
      <td>90.0</td>
      <td>1.95</td>
      <td>10.64</td>
      <td>2.67</td>
      <td>0.00</td>
    </tr>
    <tr>
      <th>Q</th>
      <td>-2.0</td>
      <td>1.0</td>
      <td>-3.5</td>
      <td>0.2</td>
      <td>5.65</td>
      <td>114.0</td>
      <td>2.17</td>
      <td>9.13</td>
      <td>1.56</td>
      <td>0.18</td>
    </tr>
    <tr>
      <th>R</th>
      <td>2.0</td>
      <td>1.0</td>
      <td>-4.5</td>
      <td>3.0</td>
      <td>10.76</td>
      <td>148.0</td>
      <td>1.82</td>
      <td>8.99</td>
      <td>2.34</td>
      <td>0.29</td>
    </tr>
    <tr>
      <th>S</th>
      <td>-2.0</td>
      <td>0.0</td>
      <td>-0.8</td>
      <td>0.3</td>
      <td>5.68</td>
      <td>73.0</td>
      <td>2.19</td>
      <td>9.21</td>
      <td>1.31</td>
      <td>0.06</td>
    </tr>
    <tr>
      <th>T</th>
      <td>-2.0</td>
      <td>0.0</td>
      <td>-0.7</td>
      <td>-0.4</td>
      <td>5.60</td>
      <td>93.0</td>
      <td>2.09</td>
      <td>9.10</td>
      <td>3.03</td>
      <td>0.11</td>
    </tr>
    <tr>
      <th>V</th>
      <td>1.0</td>
      <td>-1.0</td>
      <td>4.2</td>
      <td>-1.5</td>
      <td>6.00</td>
      <td>105.0</td>
      <td>2.39</td>
      <td>9.74</td>
      <td>3.67</td>
      <td>0.14</td>
    </tr>
    <tr>
      <th>W</th>
      <td>1.0</td>
      <td>-1.0</td>
      <td>-0.9</td>
      <td>-3.4</td>
      <td>5.89</td>
      <td>163.0</td>
      <td>2.46</td>
      <td>9.41</td>
      <td>3.21</td>
      <td>0.41</td>
    </tr>
    <tr>
      <th>Y</th>
      <td>-2.0</td>
      <td>-1.0</td>
      <td>-1.3</td>
      <td>-2.3</td>
      <td>5.64</td>
      <td>141.0</td>
      <td>2.20</td>
      <td>9.21</td>
      <td>2.94</td>
      <td>0.30</td>
    </tr>
  </tbody>
</table>
</div>



### PCA on features
There is some redundancy in these features, which we can see by examining the correlation coefficients.


```python
aa_feats.corr()
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
      <th>hydrophobicity</th>
      <th>polar</th>
      <th>hydropathy intensity</th>
      <th>hydrophilicity</th>
      <th>pH_l</th>
      <th>vdW_vol</th>
      <th>pK1</th>
      <th>pK2</th>
      <th>steric</th>
      <th>polarizability</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>hydrophobicity</th>
      <td>1.000000</td>
      <td>-0.071685</td>
      <td>0.116361</td>
      <td>0.256511</td>
      <td>0.098194</td>
      <td>0.304599</td>
      <td>0.018649</td>
      <td>0.237080</td>
      <td>-0.023572</td>
      <td>0.265580</td>
    </tr>
    <tr>
      <th>polar</th>
      <td>-0.071685</td>
      <td>1.000000</td>
      <td>-0.856103</td>
      <td>0.793689</td>
      <td>0.175193</td>
      <td>-0.167100</td>
      <td>-0.515392</td>
      <td>-0.352830</td>
      <td>-0.514872</td>
      <td>-0.220739</td>
    </tr>
    <tr>
      <th>hydropathy intensity</th>
      <td>0.116361</td>
      <td>-0.856103</td>
      <td>1.000000</td>
      <td>-0.579074</td>
      <td>-0.185182</td>
      <td>-0.116596</td>
      <td>0.459863</td>
      <td>0.544784</td>
      <td>0.452738</td>
      <td>-0.141111</td>
    </tr>
    <tr>
      <th>hydrophilicity</th>
      <td>0.256511</td>
      <td>0.793689</td>
      <td>-0.579074</td>
      <td>1.000000</td>
      <td>0.144248</td>
      <td>-0.319574</td>
      <td>-0.385994</td>
      <td>-0.071089</td>
      <td>-0.613583</td>
      <td>-0.391803</td>
    </tr>
    <tr>
      <th>pH_l</th>
      <td>0.098194</td>
      <td>0.175193</td>
      <td>-0.185182</td>
      <td>0.144248</td>
      <td>1.000000</td>
      <td>0.364222</td>
      <td>-0.214373</td>
      <td>-0.301306</td>
      <td>0.253380</td>
      <td>0.286892</td>
    </tr>
    <tr>
      <th>vdW_vol</th>
      <td>0.304599</td>
      <td>-0.167100</td>
      <td>-0.116596</td>
      <td>-0.319574</td>
      <td>0.364222</td>
      <td>1.000000</td>
      <td>-0.009948</td>
      <td>-0.400525</td>
      <td>0.567074</td>
      <td>0.939633</td>
    </tr>
    <tr>
      <th>pK1</th>
      <td>0.018649</td>
      <td>-0.515392</td>
      <td>0.459863</td>
      <td>-0.385994</td>
      <td>-0.214373</td>
      <td>-0.009948</td>
      <td>1.000000</td>
      <td>-0.057576</td>
      <td>0.081364</td>
      <td>0.057177</td>
    </tr>
    <tr>
      <th>pK2</th>
      <td>0.237080</td>
      <td>-0.352830</td>
      <td>0.544784</td>
      <td>-0.071089</td>
      <td>-0.301306</td>
      <td>-0.400525</td>
      <td>-0.057576</td>
      <td>1.000000</td>
      <td>-0.005957</td>
      <td>-0.450413</td>
    </tr>
    <tr>
      <th>steric</th>
      <td>-0.023572</td>
      <td>-0.514872</td>
      <td>0.452738</td>
      <td>-0.613583</td>
      <td>0.253380</td>
      <td>0.567074</td>
      <td>0.081364</td>
      <td>-0.005957</td>
      <td>1.000000</td>
      <td>0.484432</td>
    </tr>
    <tr>
      <th>polarizability</th>
      <td>0.265580</td>
      <td>-0.220739</td>
      <td>-0.141111</td>
      <td>-0.391803</td>
      <td>0.286892</td>
      <td>0.939633</td>
      <td>0.057177</td>
      <td>-0.450413</td>
      <td>0.484432</td>
      <td>1.000000</td>
    </tr>
  </tbody>
</table>
</div>




```python
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
```


```python
pca = PCA(n_components=3)
pca.fit(aa_feats.values)
X = pca.transform(aa_feats.values)
```


```python
fig = plt.figure(1, figsize=(8, 6))
plt.clf()
ax = Axes3D(fig)
ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=X[:, 0], cmap=plt.cm.viridis)
for i, txt in enumerate(list(aa_feats.index)):
    ax.text(X[i, 0], X[i, 1], X[i, 2], "%s" % (txt), size=15, color="k")
plt.show()
```


![png](output_33_0.png)


We remove any features that has an absolute correlation coefficient with any other feature greater than 0.6, leaving us with:


```python
low_corr_feats = ["hydrophobicity", "hydropathy intensity", "pH_l", "pK1", "pK2", "steric", "polarizability"]
aa_feats_2 = aa_feats[low_corr_feats]
```


```python
aa_feats_2
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
      <th>hydrophobicity</th>
      <th>hydropathy intensity</th>
      <th>pH_l</th>
      <th>pK1</th>
      <th>pK2</th>
      <th>steric</th>
      <th>polarizability</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>A</th>
      <td>0.0</td>
      <td>1.8</td>
      <td>6.01</td>
      <td>2.35</td>
      <td>9.87</td>
      <td>1.28</td>
      <td>0.05</td>
    </tr>
    <tr>
      <th>C</th>
      <td>1.0</td>
      <td>2.5</td>
      <td>5.05</td>
      <td>1.92</td>
      <td>10.70</td>
      <td>1.77</td>
      <td>0.13</td>
    </tr>
    <tr>
      <th>D</th>
      <td>2.0</td>
      <td>-3.5</td>
      <td>2.85</td>
      <td>1.99</td>
      <td>9.90</td>
      <td>1.60</td>
      <td>0.11</td>
    </tr>
    <tr>
      <th>E</th>
      <td>2.0</td>
      <td>-3.5</td>
      <td>3.15</td>
      <td>2.10</td>
      <td>9.47</td>
      <td>0.00</td>
      <td>0.15</td>
    </tr>
    <tr>
      <th>F</th>
      <td>1.0</td>
      <td>2.8</td>
      <td>5.49</td>
      <td>2.20</td>
      <td>9.31</td>
      <td>2.94</td>
      <td>0.29</td>
    </tr>
    <tr>
      <th>G</th>
      <td>0.0</td>
      <td>-0.4</td>
      <td>6.06</td>
      <td>2.35</td>
      <td>9.78</td>
      <td>0.00</td>
      <td>0.00</td>
    </tr>
    <tr>
      <th>H</th>
      <td>-1.0</td>
      <td>-3.2</td>
      <td>7.60</td>
      <td>1.80</td>
      <td>9.33</td>
      <td>2.99</td>
      <td>0.23</td>
    </tr>
    <tr>
      <th>I</th>
      <td>1.0</td>
      <td>4.5</td>
      <td>6.05</td>
      <td>2.32</td>
      <td>9.76</td>
      <td>4.19</td>
      <td>0.19</td>
    </tr>
    <tr>
      <th>K</th>
      <td>2.0</td>
      <td>-3.9</td>
      <td>9.60</td>
      <td>2.16</td>
      <td>9.06</td>
      <td>1.89</td>
      <td>0.22</td>
    </tr>
    <tr>
      <th>L</th>
      <td>1.0</td>
      <td>3.8</td>
      <td>6.01</td>
      <td>2.33</td>
      <td>9.74</td>
      <td>2.59</td>
      <td>0.19</td>
    </tr>
    <tr>
      <th>M</th>
      <td>1.0</td>
      <td>1.9</td>
      <td>5.74</td>
      <td>2.13</td>
      <td>9.28</td>
      <td>2.35</td>
      <td>0.22</td>
    </tr>
    <tr>
      <th>N</th>
      <td>-2.0</td>
      <td>-3.5</td>
      <td>5.41</td>
      <td>2.14</td>
      <td>8.72</td>
      <td>1.60</td>
      <td>0.13</td>
    </tr>
    <tr>
      <th>P</th>
      <td>-1.0</td>
      <td>1.6</td>
      <td>6.30</td>
      <td>1.95</td>
      <td>10.64</td>
      <td>2.67</td>
      <td>0.00</td>
    </tr>
    <tr>
      <th>Q</th>
      <td>-2.0</td>
      <td>-3.5</td>
      <td>5.65</td>
      <td>2.17</td>
      <td>9.13</td>
      <td>1.56</td>
      <td>0.18</td>
    </tr>
    <tr>
      <th>R</th>
      <td>2.0</td>
      <td>-4.5</td>
      <td>10.76</td>
      <td>1.82</td>
      <td>8.99</td>
      <td>2.34</td>
      <td>0.29</td>
    </tr>
    <tr>
      <th>S</th>
      <td>-2.0</td>
      <td>-0.8</td>
      <td>5.68</td>
      <td>2.19</td>
      <td>9.21</td>
      <td>1.31</td>
      <td>0.06</td>
    </tr>
    <tr>
      <th>T</th>
      <td>-2.0</td>
      <td>-0.7</td>
      <td>5.60</td>
      <td>2.09</td>
      <td>9.10</td>
      <td>3.03</td>
      <td>0.11</td>
    </tr>
    <tr>
      <th>V</th>
      <td>1.0</td>
      <td>4.2</td>
      <td>6.00</td>
      <td>2.39</td>
      <td>9.74</td>
      <td>3.67</td>
      <td>0.14</td>
    </tr>
    <tr>
      <th>W</th>
      <td>1.0</td>
      <td>-0.9</td>
      <td>5.89</td>
      <td>2.46</td>
      <td>9.41</td>
      <td>3.21</td>
      <td>0.41</td>
    </tr>
    <tr>
      <th>Y</th>
      <td>-2.0</td>
      <td>-1.3</td>
      <td>5.64</td>
      <td>2.20</td>
      <td>9.21</td>
      <td>2.94</td>
      <td>0.30</td>
    </tr>
  </tbody>
</table>
</div>



To make the features more amenable to training, we normalize the features except for hydrophobicity.


```python
for feat in ["hydropathy intensity", "pH_l", "pK1", "pK2", "steric", "polarizability"]:
    aa_feats_2[feat] = (aa_feats_2[feat] - aa_feats_2[feat].mean()) / aa_feats_2[feat].std()
aa_feats_2
```

    /home/dillon/.conda/envs/tf1_7/lib/python3.6/site-packages/ipykernel_launcher.py:2: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      





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
      <th>hydrophobicity</th>
      <th>hydropathy intensity</th>
      <th>pH_l</th>
      <th>pK1</th>
      <th>pK2</th>
      <th>steric</th>
      <th>polarizability</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>A</th>
      <td>0.0</td>
      <td>0.707641</td>
      <td>-0.009696</td>
      <td>1.058258</td>
      <td>0.691227</td>
      <td>-0.838612</td>
      <td>-1.144703</td>
    </tr>
    <tr>
      <th>C</th>
      <td>1.0</td>
      <td>0.940199</td>
      <td>-0.557225</td>
      <td>-1.251645</td>
      <td>2.318798</td>
      <td>-0.390254</td>
      <td>-0.381568</td>
    </tr>
    <tr>
      <th>D</th>
      <td>2.0</td>
      <td>-1.053156</td>
      <td>-1.811980</td>
      <td>-0.875614</td>
      <td>0.750055</td>
      <td>-0.545807</td>
      <td>-0.572351</td>
    </tr>
    <tr>
      <th>E</th>
      <td>2.0</td>
      <td>-1.053156</td>
      <td>-1.640877</td>
      <td>-0.284709</td>
      <td>-0.093144</td>
      <td>-2.009831</td>
      <td>-0.190784</td>
    </tr>
    <tr>
      <th>F</th>
      <td>1.0</td>
      <td>1.039867</td>
      <td>-0.306274</td>
      <td>0.252478</td>
      <td>-0.406893</td>
      <td>0.680314</td>
      <td>1.144703</td>
    </tr>
    <tr>
      <th>G</th>
      <td>0.0</td>
      <td>-0.023256</td>
      <td>0.018821</td>
      <td>1.058258</td>
      <td>0.514744</td>
      <td>-2.009831</td>
      <td>-1.621663</td>
    </tr>
    <tr>
      <th>H</th>
      <td>-1.0</td>
      <td>-0.953488</td>
      <td>0.897150</td>
      <td>-1.896269</td>
      <td>-0.367674</td>
      <td>0.726065</td>
      <td>0.572351</td>
    </tr>
    <tr>
      <th>I</th>
      <td>1.0</td>
      <td>1.604651</td>
      <td>0.013118</td>
      <td>0.897102</td>
      <td>0.475525</td>
      <td>1.824083</td>
      <td>0.190784</td>
    </tr>
    <tr>
      <th>K</th>
      <td>2.0</td>
      <td>-1.186046</td>
      <td>2.037835</td>
      <td>0.037603</td>
      <td>-0.897125</td>
      <td>-0.280452</td>
      <td>0.476960</td>
    </tr>
    <tr>
      <th>L</th>
      <td>1.0</td>
      <td>1.372093</td>
      <td>-0.009696</td>
      <td>0.950821</td>
      <td>0.436307</td>
      <td>0.360059</td>
      <td>0.190784</td>
    </tr>
    <tr>
      <th>M</th>
      <td>1.0</td>
      <td>0.740864</td>
      <td>-0.163688</td>
      <td>-0.123553</td>
      <td>-0.465720</td>
      <td>0.140455</td>
      <td>0.476960</td>
    </tr>
    <tr>
      <th>N</th>
      <td>-2.0</td>
      <td>-1.053156</td>
      <td>-0.351902</td>
      <td>-0.069834</td>
      <td>-1.563840</td>
      <td>-0.545807</td>
      <td>-0.381568</td>
    </tr>
    <tr>
      <th>P</th>
      <td>-1.0</td>
      <td>0.641196</td>
      <td>0.155704</td>
      <td>-1.090489</td>
      <td>2.201142</td>
      <td>0.433260</td>
      <td>-1.621663</td>
    </tr>
    <tr>
      <th>Q</th>
      <td>-2.0</td>
      <td>-1.053156</td>
      <td>-0.215019</td>
      <td>0.091322</td>
      <td>-0.759860</td>
      <td>-0.582407</td>
      <td>0.095392</td>
    </tr>
    <tr>
      <th>R</th>
      <td>2.0</td>
      <td>-1.385382</td>
      <td>2.699433</td>
      <td>-1.788832</td>
      <td>-1.034390</td>
      <td>0.131305</td>
      <td>1.144703</td>
    </tr>
    <tr>
      <th>S</th>
      <td>-2.0</td>
      <td>-0.156146</td>
      <td>-0.197909</td>
      <td>0.198759</td>
      <td>-0.602985</td>
      <td>-0.811161</td>
      <td>-1.049311</td>
    </tr>
    <tr>
      <th>T</th>
      <td>-2.0</td>
      <td>-0.122924</td>
      <td>-0.243536</td>
      <td>-0.338428</td>
      <td>-0.818688</td>
      <td>0.762665</td>
      <td>-0.572351</td>
    </tr>
    <tr>
      <th>V</th>
      <td>1.0</td>
      <td>1.504983</td>
      <td>-0.015399</td>
      <td>1.273133</td>
      <td>0.436307</td>
      <td>1.348275</td>
      <td>-0.286176</td>
    </tr>
    <tr>
      <th>W</th>
      <td>1.0</td>
      <td>-0.189369</td>
      <td>-0.078137</td>
      <td>1.649163</td>
      <td>-0.210800</td>
      <td>0.927368</td>
      <td>2.289406</td>
    </tr>
    <tr>
      <th>Y</th>
      <td>-2.0</td>
      <td>-0.322259</td>
      <td>-0.220723</td>
      <td>0.252478</td>
      <td>-0.602985</td>
      <td>0.680314</td>
      <td>1.240095</td>
    </tr>
  </tbody>
</table>
</div>




```python
aa_feats_2.mean(axis=0)
```




    hydrophobicity          1.500000e-01
    hydropathy intensity   -5.551115e-18
    pH_l                    3.733125e-16
    pK1                     2.375877e-15
    pK2                    -3.314016e-15
    steric                  1.110223e-16
    polarizability         -1.110223e-16
    dtype: float64




```python
aa_feats_2.std(axis=0)
```




    hydrophobicity          1.531253
    hydropathy intensity    1.000000
    pH_l                    1.000000
    pK1                     1.000000
    pK2                     1.000000
    steric                  1.000000
    polarizability          1.000000
    dtype: float64




```python
# re-add the "unknown" amino acid character, X
aa_feats_2.loc["X"] = [0 for _ in range(aa_feats_2.shape[1])]
```

    /home/dillon/.conda/envs/tf1_7/lib/python3.6/site-packages/ipykernel_launcher.py:2: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      



```python
aa_feats_2
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
      <th>hydrophobicity</th>
      <th>hydropathy intensity</th>
      <th>pH_l</th>
      <th>pK1</th>
      <th>pK2</th>
      <th>steric</th>
      <th>polarizability</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>A</th>
      <td>0.0</td>
      <td>0.707641</td>
      <td>-0.009696</td>
      <td>1.058258</td>
      <td>0.691227</td>
      <td>-0.838612</td>
      <td>-1.144703</td>
    </tr>
    <tr>
      <th>C</th>
      <td>1.0</td>
      <td>0.940199</td>
      <td>-0.557225</td>
      <td>-1.251645</td>
      <td>2.318798</td>
      <td>-0.390254</td>
      <td>-0.381568</td>
    </tr>
    <tr>
      <th>D</th>
      <td>2.0</td>
      <td>-1.053156</td>
      <td>-1.811980</td>
      <td>-0.875614</td>
      <td>0.750055</td>
      <td>-0.545807</td>
      <td>-0.572351</td>
    </tr>
    <tr>
      <th>E</th>
      <td>2.0</td>
      <td>-1.053156</td>
      <td>-1.640877</td>
      <td>-0.284709</td>
      <td>-0.093144</td>
      <td>-2.009831</td>
      <td>-0.190784</td>
    </tr>
    <tr>
      <th>F</th>
      <td>1.0</td>
      <td>1.039867</td>
      <td>-0.306274</td>
      <td>0.252478</td>
      <td>-0.406893</td>
      <td>0.680314</td>
      <td>1.144703</td>
    </tr>
    <tr>
      <th>G</th>
      <td>0.0</td>
      <td>-0.023256</td>
      <td>0.018821</td>
      <td>1.058258</td>
      <td>0.514744</td>
      <td>-2.009831</td>
      <td>-1.621663</td>
    </tr>
    <tr>
      <th>H</th>
      <td>-1.0</td>
      <td>-0.953488</td>
      <td>0.897150</td>
      <td>-1.896269</td>
      <td>-0.367674</td>
      <td>0.726065</td>
      <td>0.572351</td>
    </tr>
    <tr>
      <th>I</th>
      <td>1.0</td>
      <td>1.604651</td>
      <td>0.013118</td>
      <td>0.897102</td>
      <td>0.475525</td>
      <td>1.824083</td>
      <td>0.190784</td>
    </tr>
    <tr>
      <th>K</th>
      <td>2.0</td>
      <td>-1.186046</td>
      <td>2.037835</td>
      <td>0.037603</td>
      <td>-0.897125</td>
      <td>-0.280452</td>
      <td>0.476960</td>
    </tr>
    <tr>
      <th>L</th>
      <td>1.0</td>
      <td>1.372093</td>
      <td>-0.009696</td>
      <td>0.950821</td>
      <td>0.436307</td>
      <td>0.360059</td>
      <td>0.190784</td>
    </tr>
    <tr>
      <th>M</th>
      <td>1.0</td>
      <td>0.740864</td>
      <td>-0.163688</td>
      <td>-0.123553</td>
      <td>-0.465720</td>
      <td>0.140455</td>
      <td>0.476960</td>
    </tr>
    <tr>
      <th>N</th>
      <td>-2.0</td>
      <td>-1.053156</td>
      <td>-0.351902</td>
      <td>-0.069834</td>
      <td>-1.563840</td>
      <td>-0.545807</td>
      <td>-0.381568</td>
    </tr>
    <tr>
      <th>P</th>
      <td>-1.0</td>
      <td>0.641196</td>
      <td>0.155704</td>
      <td>-1.090489</td>
      <td>2.201142</td>
      <td>0.433260</td>
      <td>-1.621663</td>
    </tr>
    <tr>
      <th>Q</th>
      <td>-2.0</td>
      <td>-1.053156</td>
      <td>-0.215019</td>
      <td>0.091322</td>
      <td>-0.759860</td>
      <td>-0.582407</td>
      <td>0.095392</td>
    </tr>
    <tr>
      <th>R</th>
      <td>2.0</td>
      <td>-1.385382</td>
      <td>2.699433</td>
      <td>-1.788832</td>
      <td>-1.034390</td>
      <td>0.131305</td>
      <td>1.144703</td>
    </tr>
    <tr>
      <th>S</th>
      <td>-2.0</td>
      <td>-0.156146</td>
      <td>-0.197909</td>
      <td>0.198759</td>
      <td>-0.602985</td>
      <td>-0.811161</td>
      <td>-1.049311</td>
    </tr>
    <tr>
      <th>T</th>
      <td>-2.0</td>
      <td>-0.122924</td>
      <td>-0.243536</td>
      <td>-0.338428</td>
      <td>-0.818688</td>
      <td>0.762665</td>
      <td>-0.572351</td>
    </tr>
    <tr>
      <th>V</th>
      <td>1.0</td>
      <td>1.504983</td>
      <td>-0.015399</td>
      <td>1.273133</td>
      <td>0.436307</td>
      <td>1.348275</td>
      <td>-0.286176</td>
    </tr>
    <tr>
      <th>W</th>
      <td>1.0</td>
      <td>-0.189369</td>
      <td>-0.078137</td>
      <td>1.649163</td>
      <td>-0.210800</td>
      <td>0.927368</td>
      <td>2.289406</td>
    </tr>
    <tr>
      <th>Y</th>
      <td>-2.0</td>
      <td>-0.322259</td>
      <td>-0.220723</td>
      <td>0.252478</td>
      <td>-0.602985</td>
      <td>0.680314</td>
      <td>1.240095</td>
    </tr>
    <tr>
      <th>X</th>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
  </tbody>
</table>
</div>




```python
aa_feats_2.to_csv("aa_feats_final.csv")
```
