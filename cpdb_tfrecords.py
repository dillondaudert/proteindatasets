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

# Count the protein sequence lengths for all samples
def get_length(seq_labels):
    assert seq_labels.shape == (700, 9)
    noseq = np.array([[0., 0., 0., 0., 0., 0., 0., 0., 1.]])
    return np.logical_not(np.all(np.equal(seq_labels, noseq), axis=1)).sum()

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
    seqs_onehot = data[:, :, 0:22]
    ss_onehot = data[:, :, 22:31]
    seqs_inds = np.argmax(seqs_onehot, axis=2)
    ss_inds = np.argmax(ss_onehot, axis=2)
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
    # TODO: calculate phyche features

    # Get the indices for training, validation set
    train_examples = range(0, num_samples-256)
    valid_examples = range(num_samples-256, num_samples)
    print("train range: ", train_examples)
    print("valid range: ", valid_examples)

    train_file = "cpdb_train.tfrecords"
    valid_file = "cpdb_valid.tfrecords"

    print("Writing ", train_file)
    train_writer = tf.python_io.TFRecordWriter(train_file)

    for index in train_examples:
        example = tf.train.Example(features=tf.train.Features(feature={
            'dssp_id': _bytes_feature(bytes("No ID", "utf-8")),
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
            'dssp_id': _bytes_feature(bytes("No ID", "utf-8")),
            'seq_len': _int64_feature(len(seqs[index])),
            'seq': _bytes_feature(bytes(seqs[index], "utf-8")),
            'seq_phyche': _floats_feature(prot_to_vector(seqs[index]).reshape(-1)),
            'seq_pssm': _floats_feature(seqs_pssm[index, 0:len(seqs[index]), :].reshape(-1)),
            'ss': _bytes_feature(bytes(ss[index], "utf-8"))}))
        valid_writer.write(example.SerializeToString())
    valid_writer.close()

def cpdb_513_to_tfrecord(datadir: str):
    """
    Convert the numpy array format for cpdb_513 to a TFRecord file.
    """

    datadir = os.path.abspath(datadir)
    data = np.load(os.path.join(datadir, "cb513+profile_split1.npy.gz")).reshape(-1, 700, 57)
    # get indices for train/valid sets
    num_samples = data.shape[0]

    seqs = np.concatenate([data[:, :, 0:22].copy(), data[:, :, 35:56].copy()], axis=2).reshape(num_samples, -1)
    labels = data[:, :, 22:31].copy().reshape(num_samples, 700, -1)

    num_features = 43
    num_labels = 9

    seq_lengths = [get_length(labels[l, :, :]) for l in range(num_samples)]

    # Flatten labels
    labels = labels.reshape(num_samples, -1)

    filename = os.path.join(datadir, "cpdb_513.tfrecords")
    print("Writing ", filename)
    writer = tf.python_io.TFRecordWriter(filename)

    for index in range(num_samples):
        example = tf.train.Example(features=tf.train.Features(feature={
            'seq_len': _int64_feature(seq_lengths[index]),
            'seq_data': _floats_feature(seqs[index, 0:num_features*seq_lengths[index]]),
            'label_data': _floats_feature(labels[index, 0:num_labels*seq_lengths[index]])}))
        writer.write(example.SerializeToString())
    writer.close()

parser = ap.ArgumentParser(description="Convert the CPDB dataset from numpy arrays to TF records.")
parser.add_argument("-d", "--datadir", type=str, required=True,
                    help="The directory where the data will be read from and written to.")
args = parser.parse_args()

if not os.path.isdir(args.datadir):
    print("Invalid directory %s, quitting." % (args.datadir))

print("Processing data.")
cpdb_to_tfrecord(args.datadir)
#cpdb_513_to_tfrecord(args.datadir)
