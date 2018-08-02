# convert a pandas DataFrame of amino acid sequence, secondary structure
# sequence pairs into TF records
from multiprocessing import Process, Queue
from pathlib import Path
from glob import glob
import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.model_selection import KFold
from features import prot_to_vector

HOME = str(Path.home())

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
        for index in range(data.shape[0]):
            sample = data.iloc[index]
            if len(sample.seq) > 1050 or len(sample.seq) < 20:
                num_skipped += 1
                continue
            # convert the strings to
            try:
                seq_id = bytes(sample.id, "utf-8")
                seq_len = len(sample.seq)
                seq = bytes(sample.seq, "utf-8")
                seq_phyche = prot_to_vector(sample.seq).reshape(-1)

                tf_example = tf.train.Example(features=tf.train.Features(feature={
                    "uniref_id": _bytes_feature(seq_id),
                    "seq_len": _int64_feature(seq_len),
                    "seq": _bytes_feature(seq),
                    "seq_phyche": _floats_feature(seq_phyche)}))

                writer.write(tf_example.SerializeToString())
                num_written += 1

            except Exception as e:
                print("Exception encountered while processing index %d" % index)
                print(e)
                print(sample.id)
                num_skipped += 1

    return num_written, num_skipped

def worker(wid, worker_queue, done_queue):
    """
    A worker processes views of the dataframe of records and writes them to
    files.
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
    worker_queue = Queue(maxsize=10)
    done_queue = Queue()

    print("Spawning %d workers." % (num_workers))
    workers = []
    for i in range(num_workers):
        p = Process(target=worker, args=(i, worker_queue, done_queue))
        workers.append(p)
        p.start()


    files = ["uniref50_%d.csv" % (i) for i in range(1, 8)]
    # the global count of output files
    outfile_count = 0
    outfile_prefix = HOME+"/data/uniref50/tfrecords/"

    # find the last output file that was written
    last_file = Path(sorted(glob(outfile_prefix+"ur50_*.tfrecords"))[-1]).stem
    start_count = int(last_file.split("_")[-1]) - num_workers
    start_count = start_count if start_count > 0 else 0
    print("Starting at outfile #%d\n" % (start_count))

    filesize = 1000
    # for each subset of uniref50, containing a few million proteins
    for f in files:
        print("Processing %s\n" % f)
        data = pd.read_csv(HOME+"/data/uniref50/"+f)
        num_seqs = data.shape[0]
        # split into tfrecord files, each with filesize (1000) proteins
        num_outfiles = num_seqs // filesize if num_seqs % filesize == 0 else (num_seqs // filesize) + 1

        # pass views of the dataframe into the queue
        for i in range(num_outfiles):
            # NOTE: if the file already exists, skip it
            if outfile_count < start_count:
                outfile_count += 1
                continue

            outfile = outfile_prefix+"ur50_%05d.tfrecords" % (outfile_count)
            start_index = i*filesize
            end_index = (i+1)*filesize if (i+1)*filesize < num_seqs else num_seqs
            worker_queue.put((outfile, data.iloc[start_index:end_index]))
            outfile_count += 1

        print("Final index for %s: %d, written to %s" % (f, end_index, outfile))

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

