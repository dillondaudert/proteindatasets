# This script takes a fasta file as input and runs PSIBLAST
# on each one individually. The results are saved in their
# own files, with the names equal to their fasta IDs

import os
from multiprocessing import Process, Queue
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
    num_threads=6
    num_iterations=3
    # out_pssm=str(query_basepath)+"_blast.pssm"
    # out_ascii_pssm=str(query_basepath)+"_blast.ascii_pssm"
    # save_pssm_after_last_round=True
    try:
        psib_cline = NcbipsiblastCommandline(query=query, db=db, evalue=evalue, outfmt=outfmt, out=out, num_threads=num_threads, num_iterations=num_iterations)
    #print(psib_cline)
        stdout, stderr = psib_cline()
    except:
        print("Failed to run PSIBLAST on record %s" % rec.id)
        return -1

    return 1


def create_directory(rec):
    dirname = rec.id+"_pb"
    dirpath = Path(Path.cwd(), dirname)

    # check that directory doesn't exist
    count = 0
    if dirpath.exists():
        return None
#        count += 1
#        print("dir %s exists, incrementing name" % (dirname))
#        dirname = rec.id+"_pb_"+str(count)
#        dirpath = Path(Path.cwd(), dirname)

    # create directory
    #print("creating dir %s" % str(dirpath))
    os.makedirs(str(dirpath))

    return dirpath

def worker(rec_queue, done_queue):
    #
    count = 0
    skipped = 0
    while True:
        rec = rec_queue.get()
        if rec is None:
            done_queue.put((count, skipped))
            return

        count += 1

        # check if dir / files already exist
        # create directory for this record
        dirpath = create_directory(rec)
        if dirpath is None:
            skipped += 1
            continue
        # execute psiblast
        if do_psiblast(dirpath, rec) == -1:
            skipped += 1



if __name__ == "__main__":

    num_workers = 3

    records = SeqIO.parse("/home/dillonbbailly/main/uniref50/uniref50_filt.fasta", "fasta")

    rec_queue = Queue(1000)
    done_queue = Queue()

    workers = []
    for i in range(num_workers):
        p = Process(target=worker, args=(rec_queue, done_queue))
        workers.append(p)
        p.start()

    for i, rec in enumerate(records):
        if len(rec.seq) < 25 or len(rec.seq) > 2000:
            continue
        rec_queue.put(rec)
        if i % 5000 == 0:
            print("Handled %d records" % (i))

    for i in range(num_workers):
        rec_queue.put(None)

    total = 0
    skipped = 0
    for i in range(num_workers):
        (count, skip) = done_queue.get()
        total += count
        skipped += skip

    print("DONE: %d sequences processed, of which %d were skipped." % (total, skipped))

    for p in workers:
        p.join()





