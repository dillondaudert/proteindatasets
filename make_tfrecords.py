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

def cpdb2_to_tfrecord():
    """
    Convert a pandas dataframe of protein, structure string pairs to TFRecord
    format.
    Creates 5 pairs of training, validation data using 5-fold cross validation.
    """

    data = pd.read_csv(HOME+"/data/dssp/cpdb_dssp_14335.csv")

    folds = 5

    kf = KFold(folds)


    curr_fold = 0

    for train_inds, valid_inds in kf.split(np.arange(data.shape[0])):
        curr_fold += 1
        print("Creating fold %d: Train %d, Valid %d\n" % (curr_fold, train_inds.size, valid_inds.size))

        train_file = HOME+"/data/cpdb2/cpdb2_14335_train_%d.tfrecords" % curr_fold
        valid_file = HOME+"/data/cpdb2/cpdb2_14335_valid_%d.tfrecords" % curr_fold

        print("Writing ", train_file)
        train_writer = tf.python_io.TFRecordWriter(train_file)
        for index in train_inds:
            # convert the strings to
            try:
                sample = data.iloc[index]
                seq_id = bytes(sample.dssp_id, "utf-8")
                seq_data = prot_to_vector(sample.seq, kind="aa").reshape(-1)
                label_data = prot_to_vector(sample.ss, kind="ss").reshape(-1)
                seq_len = len(sample.seq)

                tf_example = tf.train.Example(features=tf.train.Features(feature={
                    "dssp_id": _bytes_feature(seq_id),
                    "seq_len": _int64_feature(seq_len),
                    "seq_data": _floats_feature(seq_data),
                    "label_data": _floats_feature(label_data)}))

                train_writer.write(tf_example.SerializeToString())

            except Exception as e:
                print("Exception encountered while processing index %d" % index)
                print(e)
                print(sample.dssp_id)
#                train_writer.close()
#                quit()


        train_writer.close()

        #
        print("Writing ", valid_file)
        valid_writer = tf.python_io.TFRecordWriter(valid_file)
        for index in valid_inds:
            # convert the strings to
            try:
                sample = data.iloc[index]
                seq_id = bytes(sample.dssp_id, "utf-8")
                seq_data = prot_to_vector(sample.seq, kind="aa").reshape(-1)
                label_data = prot_to_vector(sample.ss, kind="ss").reshape(-1)
                seq_len = len(sample.seq)

                tf_example = tf.train.Example(features=tf.train.Features(feature={
                    "dssp_id": _bytes_feature(seq_id),
                    "seq_len": _int64_feature(seq_len),
                    "seq_data": _floats_feature(seq_data),
                    "label_data": _floats_feature(label_data)}))

                valid_writer.write(tf_example.SerializeToString())

            except Exception as e:
                print("Exception encountered while processing index %d" % index)
                print(e)
                print(sample.dssp_id)
#                valid_writer.close()
#                quit()


        valid_writer.close()


    print("Finished.")

if __name__ == "__main__":
    cpdb2_to_tfrecord()

