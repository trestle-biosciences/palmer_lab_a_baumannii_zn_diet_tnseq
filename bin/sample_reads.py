#!/usr/bin/env python

# sample command:  sample_reads.py A4_S9_L003_R1_001.fastq.gz

import gzip
import argparse
import logging
import sys
from pathlib import Path
from Bio import SeqIO
from datetime import datetime
import random

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    description="Peek at the fastq reads"
    )

parser.add_argument(
    "fq_in",
    metavar="FILE_IN",
    type=Path,
    help="fastq, fq, fq.gz, fastq.gz",
)

args = parser.parse_args()
file_in = args.fq_in
logger.info(f"Processing file {file_in} starting at {datetime.now()}")

# if str(file_in).lower()[-2:] == 'gz':
#     with gzip.open(file_in, "rt") as handle:
#         for record in SeqIO.parse(handle, "fastq"):
#             if random.uniform(0,1) <= 0.0001:
#                 print(record.seq)
# else:
#     for record in SeqIO.parse(file_in, "fastq"):
#         if random.uniform(0,1) <= 0.0001:
#             print(record.seq)

### PORTION WITH ADAPTER
count_without = 0
count_with = 0

if str(file_in).lower()[-2:] == 'gz':
    with gzip.open(file_in, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if str(record.seq).find('ACAGGTTG') != -1:
                count_with +=1
            else:
                count_without +=1
                
else:
    for record in SeqIO.parse(file_in, "fastq"):
        if str(record.seq).find('ACAGGTTG') != -1:
            count_with +=1
        else:
            count_without +=1

print(count_with)
print(count_without)
print(count_with/(count_with + count_without))

# logger.critical(f"The given sample sheet does not appear to contain a header.")
#         sys.exit(1)