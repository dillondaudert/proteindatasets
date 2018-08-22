# Using BioPython and PSIBLAST

`do_psiblast` calls PSIBLAST on multiple sequences in a single fasta file and
save the results in individual directories. These directories are traversed
by `process_psiblast` to aggregate sequences into a single fasta file with
some filters applied.

These are just examples and aren't readily generalizable.
