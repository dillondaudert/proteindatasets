from Bio import SeqIO
from multiprocessing import Process, Queue
from pathlib import Path
from glob import glob
import argparse as ap
import numpy as np
import tensorflow as tf
from proteinfeatures.features import prot_to_vector

HOME = Path.home()
DATADIR = str(Path(HOME, "data", "cUR50"))

def _bytes_feature(value):
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))

def _int64_feature(value):
    return tf.train.Feature(int64_list=tf.train.Int64List(value=[value]))

def _floats_feature(value):
    return tf.train.Feature(float_list=tf.train.FloatList(value=value))

def write_file(filename, data):
    """
    Write a dataframe of records to a file.
    """
    num_written = 0
    num_skipped = 0
    with tf.python_io.TFRecordWriter(filename) as writer:
        for rec in data:
            # convert the strings to
            try:
                seq_id = bytes(rec.id, "utf-8")
                seq_len = len(str(rec.seq))
                seq = bytes(str(rec.seq), "utf-8")
                seq_phyche = prot_to_vector(str(rec.seq)).reshape(-1)

                tf_example = tf.train.Example(features=tf.train.Features(feature={
                    "id": _bytes_feature(seq_id),
                    "seq_len": _int64_feature(seq_len),
                    "seq": _bytes_feature(seq),
                    "seq_phyche": _floats_feature(seq_phyche)}))

                writer.write(tf_example.SerializeToString())
                num_written += 1

            except Exception as e:
                print("Exception encountered while processing id %s" % rec.id)
                print(e)
                num_skipped += 1

    return num_written, num_skipped

def worker(wid, worker_queue, done_queue):
    """
    A worker will read tuples of records from the `worker_queue`, calculate
    features, and write those records to a file.
    It also tracks the total number of records written and skipped.
    """

    files_written = 0
    total_written = 0
    total_skipped = 0
    while True:
        (filename, data) = worker_queue.get()
        # check if done
        if filename is None:
            done_queue.put((total_written, total_skipped))
            return

        if total_written % 10000 == 0:
            print(data[0].id, len(str(data[0].seq)), str(data[0].seq))

        written, skipped = write_file(filename, data)
        files_written += 1
        total_written += written
        total_skipped += skipped
        if files_written % 5 == 0:
            print("Worker %d:: %d total files, last file: %s \
                   \n\t- records written: %d, records skipped: %d\n" % (wid, files_written, filename, total_written, total_skipped))



def cUR50_to_tfrecords():
    """
    Convert a pandas dataframe of protein sequences into a TFRecord
    format.
    """

    num_workers = 5
    worker_queue = Queue(maxsize=10000)
    done_queue = Queue()

    print("Spawning %d workers." % (num_workers))
    workers = []
    for i in range(num_workers):
        p = Process(target=worker, args=(i, worker_queue, done_queue))
        workers.append(p)
        p.start()


    filename = DATADIR+"/cur50.fasta"
    # the global count of output files
    outfile_count = 0
    outfile_prefix = DATADIR+"/tfrecords"

    # find the last output file that was written
    sorted_files = sorted(glob(outfile_prefix+"/cur50_*.tfrecords"))
    if len(sorted_files) > 0:
        last_file = Path(sorted_files[-1]).stem
        start_count = int(last_file.split("_")[-1]) - 2*num_workers
        start_count = start_count if start_count > 0 else 0
    else:
        start_count = 0
    print("Starting at outfile #%d\n" % (start_count))

    filesize = 1000

    print("Processing %s\n" % filename)

    # read records from fasta file and yield `filesize` number at a time
    def recs_producer():
        producer = SeqIO.parse(filename, "fasta")
        while True:
            records = []
            try:
                for _ in range(filesize):
                    rec = next(producer)
                    records.append(rec)
                yield tuple(records)
            except StopIteration:
                break

        if len(records) != 0:
            yield tuple(records)


    for filerecs in recs_producer():
        # NOTE: if the file already exists, skip it
        if outfile_count < start_count:
            outfile_count += 1
            continue

        if len(filerecs) != filesize:
            print("# Recs in Tuple is not %d; instead got %d" % (filesize, len(filerecs)))

        outfile = outfile_prefix+"/cur50_%05d.tfrecords" % (outfile_count)
        worker_queue.put((outfile, filerecs))
        outfile_count += 1

    # pass stop signal to workers
    for _ in range(num_workers):
        worker_queue.put((None, None))

    total_written = 0
    total_skipped = 0
    for _ in range(num_workers):
        (records_written, records_skipped) = done_queue.get()
        total_written += records_written
        total_skipped += records_skipped

    print("%d records written, %d records skipped" % (total_written, total_skipped))

    print("Joining workers")
    for p in workers:
        p.join()


if __name__ == "__main__":
    cUR50_to_tfrecords()
