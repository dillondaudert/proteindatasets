
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from multiprocessing import Process, Queue

num_readers = 6

def del_dirs(cleanup_queue: Queue):
    """
    Check the directories passed to the queue. If it is a subdirectory
    of raw/, then delete it recursively.
    Only send directories that are safe to delete here.
    """

    while True:
        directory = cleanup_queue.get()
        # check if done; this only happens is combine_records is done
        if directory is None:
            break

        # delete
        print("Would delete %s" % str(directory))

def combine_records(save_queue: Queue, cleanup_queue: Queue):
    """
    Process psiblast results, aggregate, and save to a fasta file.
    Once written, send directories to cleanup.
    """

    readers_done = 0
    file_counter = 1
    seqrecs = []
    dirs = []
    while readers_done < num_readers:
        directory = save_queue.get()
        # check if done
        if directory is None:
            readers_done += 1
            print("Reader %d done!" % readers_done)
            continue

        dirs.append(directory)

        # read the sequence from the fasta file in this directory
        # create a SeqRecord TODO
        seqrec = None

        seqrecs.append(seqrec)

        if len(seqrecs) == 100000:
            # write out to file
            try:
                #SeqIO.write(seqrecs, "cull_uniref50_%d.fasta"%(file_counter), "fasta")
                print("Would write %d seqs to file" % len(seqrecs))
                file_counter += 1
                seqrecs.clear()
            except Exception as e:
                print("Exception: %s" % str(e))
                print("Can't write file, but directories not deleted.")
                quit()

            while len(dirs) > 0:
                cleanup_queue.put(dirs.pop())

            assert len(seqrecs) == 0
            assert len(dirs) == 0

    try:
        #SeqIO.write(seqrecs, "cull_uniref50_%d.fasta"%(file_counter), "fasta")
        print("Would write %d seqs to file" % len(seqrecs))
        file_counter += 1
        seqrecs.clear()
    except Exception as e:
        print("Exception: %s" % str(e))
        print("Can't write file, but directories not deleted.")
        quit()

    while len(dirs) > 0:
        cleanup_queue.put(dirs.pop())

    cleanup_queue.put(None)

def reader(dir_queue: Queue, cleanup_queue: Queue, save_queue: Queue):
    """
    Calculate sequence identities from blast results in directories
    passed to a queue. Forward the low % ones to be saved, and the
    rest forward to be cleaned up.
    """

    while True:
        d = dir_queue.get()

        if d is None:
            # signal to writer this reader is done
            save_queue.put(None)
            break

        # find xml
        for f in d.iterdir():
            if f.suffix == ".xml":
                try:
                    recs = list(NCBIXML.parse(open(f)))
                except Exception as e:
                    print("Exception %s encountered for file %s" % (str(e), str(f)))


                    # TODO
                    # calc identity
                    # decide where it goes


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
