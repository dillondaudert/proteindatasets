# This script takes a fasta file as input and runs PSIBLAST
# on each one individually. The results are saved in their
# own files, with the names equal to their fasta IDs

import os
from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbipsiblastCommandline

def do_psiblast(dirpath, rec):
    """
    Run a PSIBLAST query on the given sequence.
    """

    # save the query to a fasta file
    query_basepath = Path(dirpath, rec.id)
    SeqIO.write(rec, str(query_basepath)+".fasta", "fasta")

    # build the query
    query=str(query_basepath)+".fasta"
    db="cpdb2_db"
    evalue=0.001
    outfmt=5
    out=str(query_basepath)+"_blast.xml"
    num_threads=8
    num_iterations=3
    # out_pssm=str(query_basepath)+"_blast.pssm"
    # out_ascii_pssm=str(query_basepath)+"_blast.ascii_pssm"
    # save_pssm_after_last_round=True
    psib_cline = NcbipsiblastCommandline(query=query, db=db, evalue=evalue, outfmt=outfmt, out=out, num_threads=num_threads, num_iterations=num_iterations)
    #print(psib_cline)
    stdout, stderr = psib_cline()

    return stdout, stderr


def create_directory(rec):
    dirname = rec.id+"_pb"
    dirpath = Path(Path.cwd(), dirname)

    # check that directory doesn't exist
    count = 0
    while dirpath.exists():
        count += 1
        print("dir %s exists, incrementing name" % (dirname))
        dirname = rec.id+"_pb_"+str(count)
        dirpath = Path(Path.cwd(), dirname)

    # create directory
    #print("creating dir %s" % str(dirpath))
    os.makedirs(str(dirpath))

    return dirpath

def parse_results():
    # an example of how to parse a single file
    file = open("/path/to/file")
    record = NCBIXML.read(file)
    # record.alignments
    # record.query_length
    # record.num_hits
    # etc




if __name__ == "__main__":
    # open the file

    records = SeqIO.parse("/home/dillonbbailly/main/uniref50/uniref50_filt.fasta", "fasta")

    for i, rec in enumerate(records):
        if i % 1000 == 999:
            print("Handled %d records." % i)
        # check if dir / files already exist
        # create directory for this record
        dirpath = create_directory(rec)
        # execute psiblast
        do_psiblast(dirpath, rec)
