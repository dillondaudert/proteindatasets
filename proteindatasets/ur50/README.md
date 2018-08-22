# Creating a dataset from uniref50
These scripts create a tfrecords dataset from the uniref50 fasta file using
multiprocessing. To make this less RAM intensive, the fasta file is split up
until several CSVs using `fasta_to_csv.jl`.

More description to come.
