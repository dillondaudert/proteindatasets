"""Parsers for tf records files."""

from pathlib import Path
import tensorflow as tf, numpy as np

def cpdb_parser(record):
    """
    Parse a CPDB tfrecord Record into a tuple of tensors.
    """

    keys_to_features = {
        "dssp_id": tf.VarLenFeature(tf.string),
        "seq_len": tf.FixedLenFeature([], tf.int64),
        "seq": tf.VarLenFeature(tf.string),
        "seq_phyche": tf.VarLenFeature(tf.float32),
        "seq_pssm": tf.VarLenFeature(tf.float32),
        "ss": tf.VarLenFeature(tf.string),
        }

    parsed = tf.parse_single_example(record, keys_to_features)

    dssp_id = parsed["dssp_id"]
    dssp_id = tf.cast(dssp_id, tf.string)
    seq_len = parsed["seq_len"]
    seq_len = tf.cast(seq_len, tf.int32)
    seq = parsed["seq"]
    seq = tf.cast(seq, tf.string)
    seq_phyche = tf.sparse_tensor_to_dense(parsed["seq_phyche"])
    seq_pssm = tf.sparse_tensor_to_dense(parsed["seq_pssm"])
    ss = parsed["ss"]
    ss = tf.cast(ss, tf.string)

    return dssp_id, seq_len, seq, seq_phyche, seq_pssm, ss

def cpdb_dataset(tfrecords):
    """
    Open a tfrecords file in the cpdb format, parse, and
    return a tf.data.Dataset object
    """

    dataset = tf.data.TFRecordDataset(tfrecords)
    dataset = dataset.map(lambda x: cpdb_parser(x))
    return dataset
