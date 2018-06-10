# convert a pandas DataFrame of amino acid sequence, secondary structure
# sequence pairs into TF records
from pathlib import Path
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

def cUR50_to_tfrecords():
    """
    Convert a pandas dataframe of protein sequences into a TFRecord
    format.
    Create a training and validation set 90/10 split.
    """

    data = pd.read_csv(HOME+"/data/cullpdb2_psiblast/cullUR50_20pc_189K.csv")

    num_seqs = data.shape[0]
    perm = np.random.permutation(num_seqs)
    train_inds = perm[0:int(num_seqs*.9)]
    valid_inds = perm[int(num_seqs*.9):]

    train_file = HOME+"/data/cullpdb2_psiblast/cUR50_train.tfrecords"
    valid_file = HOME+"/data/cullpdb2_psiblast/cUR50_valid.tfrecords"

    print("Writing ", train_file)
    train_writer = tf.python_io.TFRecordWriter(train_file)
    for index in train_inds:
        # convert the strings to
        try:
            sample = data.iloc[index]
            seq_id = bytes(sample.uniref_id, "utf-8")
            seq_len = len(sample.seq)
            seq = bytes(sample.seq, "utf-8")
            seq_phyche = prot_to_vector(sample.seq).reshape(-1)

            tf_example = tf.train.Example(features=tf.train.Features(feature={
                "uniref_id": _bytes_feature(seq_id),
                "seq_len": _int64_feature(seq_len),
                "seq": _bytes_feature(seq),
                "seq_phyche": _floats_feature(seq_phyche)}))

            train_writer.write(tf_example.SerializeToString())

        except Exception as e:
            print("Exception encountered while processing index %d" % index)
            print(e)
            print(sample.uniref_id)
    train_writer.close()

    print("Writing ", valid_file)
    valid_writer = tf.python_io.TFRecordWriter(valid_file)
    for index in valid_inds:
        # convert the strings to
        try:
            sample = data.iloc[index]
            seq_id = bytes(sample.uniref_id, "utf-8")
            seq_len = len(sample.seq)
            seq = bytes(sample.seq, "utf-8")
            seq_phyche = prot_to_vector(sample.seq)

            tf_example = tf.train.Example(features=tf.train.Features(feature={
                "uniref_id": _bytes_feature(seq_id),
                "seq_len": _int64_feature(seq_len),
                "seq": _bytes_feature(seq),
                "seq_phyche": _floats_feature(seq_phyche)}))

            valid_writer.write(tf_example.SerializeToString())

        except Exception as e:
            print("Exception encountered while processing index %d" % index)
            print(e)
            print(sample.uniref_id)
    valid_writer.close()


if __name__ == "__main__":
    cUR50_to_tfrecords()

