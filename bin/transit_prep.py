#!/usr/bin/env python

import pandas as pd 
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(
description="Generate everything needed to run Transit Zinb for each design and each contig"
)

parser.add_argument(
    "--tests-to-run",
    metavar="TESTS_TO_RUN",
    type=Path,
    help="CSV file with columns design_id and ref_condition minimally",
)
parser.add_argument(
    "--design-file",
    metavar="DESIGN_FILE",
    type=Path,
    help="CSV file with columns design_id, sample, and condition",
)
parser.add_argument(
    "--fastas-and-gffs",
    metavar="FASTAS_AND_GFFS_FILE",
    type=Path,
    help="CSV file with coluns name, fasta, and gff",
)

parser.add_argument(
    "--outfile",
    metavar="OUT_FILE",
    type=Path,
    help="CSV file with all tests to run and associated wig files",
)


args = parser.parse_args()

design_file = Path(args.design_file)
tests_to_run = Path(args.tests_to_run)
fastas_and_gffs = Path(args.fastas_and_gffs)
outfile = Path(args.outfile)

design_file_df = pd.read_csv(design_file)
tests_to_run_df = pd.read_csv(tests_to_run)
fastas_and_gffs_df = pd.read_csv(fastas_and_gffs)

# get all contig names
if "name" not in fastas_and_gffs_df.columns:
    AssertionError("column name must appear in the file")
else: 
    all_contigs = list(fastas_and_gffs_df.name)

# ensure design file has all designs in tests_to_run
tests_to_run_designs = set(tests_to_run_df.design_id)
design_file_designs = set(design_file_df.design_id)

if not tests_to_run_designs.issubset(design_file_designs):
    AssertionError("Not all designs are included in the designs.csv file")

contig_designs = []
# get all tests to run
for contig in all_contigs:
    for index, row in tests_to_run_df.iterrows():
        design = row['design_id']
        ref_condition = row['ref_condition']
        try:
            covariates = row['covariates'].lower()
        except KeyError:
            print("No covariates found")
            covariates = None
        else:
            print("Something unexpected went wrong with covariates")
            covariates = None 
        if covariates is not None:
            designs_sub_columns = ['sample', 'condition','Filename', 'replicate'] + covariates.split(",")
        else:
            designs_sub_columns = ['sample', 'condition','Filename', 'replicate']
        if not set(designs_sub_columns).issubset(design_file_df.columns):
            AssertionError("Not all required columns are in the designs file. Check that any designated covariates have a column.")        
        designs_sub = design_file_df.loc[design_file_df.design_id==design,:].copy()
        designs_sub['Filename'] = designs_sub['sample'].apply(lambda x: f"{contig}_{x}.wig")
        wig_files = '-'.join(designs_sub['Filename'])
        design_file = f"{contig}_{design}_design.txt"
        designs_sub = designs_sub.loc[:,designs_sub_columns]
        designs_sub.columns = ['ID', 'Condition', 'Filename', 'Replicate']
        designs_sub.to_csv(f"{contig}_{design}_design.txt", sep="\t", index=False)
        contig_designs.append([contig, design, design_file, wig_files, ref_condition, covariates])

all_designs_df = pd.DataFrame.from_records(contig_designs, columns=['contig', 'design_id', 'design_file', 'wigs', 'ref_condition', 'covariates'])

all_designs_df.to_csv(outfile, index=False)