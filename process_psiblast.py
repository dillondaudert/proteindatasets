# This script takes a fasta file as input and runs PSIBLAST
# on each one individually. The results are saved in their
# own files, with the names equal to their fasta IDs

import os
from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIXML

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

def calc_identities(record):
    # calculate the % identity of a record with all of its
    # alignments; hsps
    identities = []
    for align in record.alignments:
        ident = sum((hsp.identities for hsp in align.hsps))
        identities.append(ident/record.query_length)
    return identities

def count_identities(raw_dir):
    perc30 = 0
    perc50 = 0
    perc80 = 0
    sum_perc = 0
    num_align = 0
    for i, rec in enumerate(read_records(raw_dir)):
        if i % 10000 == 0:
            print("%-9d: 30%%: %-5d, 50%%: %-5d, 80%%: %-5d" % (i, perc30, perc50, perc80))
            if num_align > 0:
                print("\tAverage %% identity: %.4f" % (sum_perc/num_align))
        identities = calc_identities(rec)
        if len(identities) > 0:
            sum_perc += sum(identities)
            num_align += len(rec.alignments)
            if max(identities) > .3:
                perc30 += 1
            if max(identities) > .5:
                perc50 += 1
            if max(identities) > .8:
                perc80 += 1
    return perc30, perc50, perc80, sum_perc, num_align



def read_records(raw_dir: Path):
    # traverse the directories of psiblast results
    for d in raw_dir.iterdir():
        for f in d.iterdir():
            if f.suffix == ".xml":
                try:
                    recs = list(NCBIXML.parse(open(f)))
                except Exception as e:
                    print("Exception %s encountered for file %s" % (str(e), str(f)))
                    continue

                for rec in recs:
                    yield rec



if __name__ == "__main__":

    raw_dir = Path(Path.home(), "main", "cull_uniref50", "raw")
    perc30, perc50, perc80, sum_perc, num_align = count_identities(raw_dir)
