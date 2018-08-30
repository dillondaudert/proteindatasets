#!/usr/bin/env python3
import os, sys
import argparse as ap
from .cpdb_tfrecords import cpdb_to_tfrecord, cpdb_513_to_tfrecord

parser = ap.ArgumentParser(description="Convert the CUR50 dataset from fasta to tfrecord format.")
parser.add_argument("infile", type=str,
                    help="The file where FASTA records will be read from")
parser.add_argument("outdir", type=str,
                    help="The directory where tfrecords files will be written.")
parser.add_argument("-np", "--num_procs", type=int, default=1,
                    help="The number of worker processes to use.")
parser.add_argument("-sz", "--file_size", type=int, default=1000,
                    help="The number of records to include in each tfrecord file.")
args = parser.parse_args()

if not os.path.isdir(args.datadir):
    print("Invalid directory %s, quitting." % (args.datadir))

print("Processing data.")
cpdb_to_tfrecord(args.datadir)
cpdb_513_to_tfrecord(args.datadir)
