# Protein sequence and structure datasets

This repo contains scripts for creating various protein sequence and structure
datasets, as well as some guides for how to use them. 

## Contents

### proteinfeatures
Protein amino acid features.

### cpdb
Working with the cullPDB dataset created in [Zhou & Troyanskaya, 
2014](https://arxiv.org/abs/1403.1347).

### cpdb2
Creating a new protein sequence-structure dataset following the methods used for
the cullPDB dataset, referred to as cpdb2.

### psiblast
Scripts for calling NCBI+ psiblast on large fasta files from BioPython and
handling the results using multiprocessing.
