#!/usr/bin/env python3
import os, sys
import argparse as ap
from cpdb_tfrecords import cpdb_to_tfrecord, cpdb_513_to_tfrecord

parser = ap.ArgumentParser(description="Convert the CPDB dataset from numpy arrays to TF records.")
parser.add_argument("-d", "--datadir", type=str, required=True,
                    help="The directory where the data will be read from and written to.")
args = parser.parse_args()

if not os.path.isdir(args.datadir):
    print("Invalid directory %s, quitting." % (args.datadir))

print("Processing data.")
cpdb_to_tfrecord(args.datadir)
cpdb_513_to_tfrecord(args.datadir)
