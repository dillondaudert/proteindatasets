#!/usr/bin/env python3
import os, sys
import argparse as ap
import numpy as np, tensorflow as tf
from features import prot_to_vector

def _int64_feature(value):
    return tf.train.Feature(int64_list=tf.train.Int64List(value=[value]))

def _bytes_feature(value):
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))

def _floats_feature(value):
    return tf.train.Feature(float_list=tf.train.FloatList(value=value))

aa_list = ["A", "C", "E", "D", "G", "F", "I", "H", "K", "M", "L", "N", "Q",
           "P", "S", "R", "T", "W", "V", "Y", "X", "NoSeq"]
ss_list = ["L", "B", "E", "G", "I", "H", "S", "T", "NoSeq"]


def cpdb_to_tfrecord(datadir: str):
    """
    Convert the numpy array format for cpdb files to TFRecord format
    Save training and validation set with 256 samples
    Args:
        datadir: the directory where the data is located. Saves the tfrecords here.
    """
    # convert onehot to string for both AA and SS
    # NOTE: the ordering of these comes from the dataset readme
    aa_list = ["A", "C", "E", "D", "G", "F", "I", "H", "K", "M", "L", "N", "Q",
               "P", "S", "R", "T", "W", "V", "Y", "X", "NoSeq"]
    ss_list = ["L", "B", "E", "G", "I", "H", "S", "T", "NoSeq"]

    datadir = os.path.abspath(datadir)
    data = np.load(os.path.join(datadir, "cpdb+profile_6133_filtered.npy.gz")).reshape(-1, 700, 57)
    num_samples = data.shape[0]

    # shuffle data
    data = np.random.permutation(data)

    # NOTE: convert one-hot sequence to string
    seqs_inds = np.argmax(data[:, :, 0:22], axis=2)
    ss_inds = np.argmax(data[:, :, 22:31], axis=2)
    seqs = []
    ss = []
    for i in range(seqs_inds.shape[0]):
        # convert the indices to letters, ignoring noseqs
        seq = "".join([aa_list[seqs_inds[i,j]] for j in range(seqs_inds.shape[1]) if aa_list[seqs_inds[i,j]] != "NoSeq"])
        ss_labels = "".join([ss_list[ss_inds[i,j]] for j in range(ss_inds.shape[1]) if ss_list[ss_inds[i,j]] != "NoSeq"])
        try:
            assert len(seq) == len(ss_labels)
        except:
            print(seq)
            print(ss_labels)
            raise
        seqs.append(seq)
        ss.append(ss_labels)

    seqs_pssm = data[:, :, 35:56]

    # Get the indices for training, validation set
    train_examples = range(0, num_samples-256)
    valid_examples = range(num_samples-256, num_samples)
    print("train range: ", train_examples)
    print("valid range: ", valid_examples)

    train_file = "TESTcpdb_train.tfrecords"
    valid_file = "TESTcpdb_valid.tfrecords"

    print("Writing ", train_file)
    train_writer = tf.python_io.TFRecordWriter(train_file)

    for index in train_examples:
        example = tf.train.Example(features=tf.train.Features(feature={
            'dssp_id': _bytes_feature(bytes(str(index), "utf-8")),
            'seq_len': _int64_feature(len(seqs[index])),
            'seq': _bytes_feature(bytes(seqs[index], "utf-8")),
            'seq_phyche': _floats_feature(prot_to_vector(seqs[index]).reshape(-1)),
            'seq_pssm': _floats_feature(seqs_pssm[index, 0:len(seqs[index]), :].reshape(-1)),
            'ss': _bytes_feature(bytes(ss[index], "utf-8"))}))
        train_writer.write(example.SerializeToString())
    train_writer.close()

    print("Writing ", valid_file)
    valid_writer = tf.python_io.TFRecordWriter(valid_file)
    for index in valid_examples:
        example = tf.train.Example(features=tf.train.Features(feature={
            'dssp_id': _bytes_feature(bytes(str(index), "utf-8")),
            'seq_len': _int64_feature(len(seqs[index])),
            'seq': _bytes_feature(bytes(seqs[index], "utf-8")),
            'seq_phyche': _floats_feature(prot_to_vector(seqs[index]).reshape(-1)),
            'seq_pssm': _floats_feature(seqs_pssm[index, 0:len(seqs[index]), :].reshape(-1)),
            'ss': _bytes_feature(bytes(ss[index], "utf-8"))}))
        valid_writer.write(example.SerializeToString())
    valid_writer.close()

# TODO: Update this function
def cpdb_513_to_tfrecord(datadir: str):
    """
    Convert the numpy array format for cpdb_513 to a TFRecord file.
    """

    datadir = os.path.abspath(datadir)
    data = np.load(os.path.join(datadir, "cb513+profile_split1.npy.gz")).reshape(-1, 700, 57)

    # NOTE: convert one-hot sequence to string
    seqs_inds = np.argmax(data[:, :, 0:22], axis=2)
    ss_inds = np.argmax(data[:, :, 22:31], axis=2)
    seqs = []
    ss = []
    for i in range(seqs_inds.shape[0]):
        # convert the indices to letters, ignoring noseqs
        seq = "".join([aa_list[seqs_inds[i,j]] for j in range(seqs_inds.shape[1]) if aa_list[seqs_inds[i,j]] != "NoSeq"])
        ss_labels = "".join([ss_list[ss_inds[i,j]] for j in range(ss_inds.shape[1]) if ss_list[ss_inds[i,j]] != "NoSeq"])
        try:
            assert len(seq) == len(ss_labels)
        except:
            print(seq)
            print(ss_labels)
            raise
        seqs.append(seq)
        ss.append(ss_labels)

    seqs_pssm = data[:, :, 35:56]


    test_file = "cpdb513_test.tfrecords"

    print("Writing ", test_file)
    test_writer = tf.python_io.TFRecordWriter(test_file)

    for index in range(data.shape[0]):
        example = tf.train.Example(features=tf.train.Features(feature={
            'dssp_id': _bytes_feature(bytes(str(index), "utf-8")),
            'seq_len': _int64_feature(len(seqs[index])),
            'seq': _bytes_feature(bytes(seqs[index], "utf-8")),
            'seq_phyche': _floats_feature(prot_to_vector(seqs[index]).reshape(-1)),
            'seq_pssm': _floats_feature(seqs_pssm[index, 0:len(seqs[index]), :].reshape(-1)),
            'ss': _bytes_feature(bytes(ss[index], "utf-8"))}))
        test_writer.write(example.SerializeToString())
    test_writer.close()

parser = ap.ArgumentParser(description="Convert the CPDB dataset from numpy arrays to TF records.")
parser.add_argument("-d", "--datadir", type=str, required=True,
                    help="The directory where the data will be read from and written to.")
args = parser.parse_args()

if not os.path.isdir(args.datadir):
    print("Invalid directory %s, quitting." % (args.datadir))

print("Processing data.")
cpdb_to_tfrecord(args.datadir)
#cpdb_513_to_tfrecord(args.datadir)
