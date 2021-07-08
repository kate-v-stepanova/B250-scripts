#!/usr/bin/env python

import sys
import os
import glob
import pandas as pd

project_id = sys.argv[1]

sample_file = None
if len(sys.argv) >= 3:
    sample_file = sys.argv[2]
    if sample_file in ['all_samples', 'all']:
        sample_file = None

mem='30G'
if len(sys.argv) >= 4:
    mem = sys.argv[3]


BASE_DIR = os.getenv('BASE_DIR')
DIRICORE_PATH = os.path.join(BASE_DIR, "software/diricore")

project_dir = os.path.join(BASE_DIR, project_id)
input_dir = os.path.join(project_dir, "analysis/output/rrna/blat_results")
fastq_dir = os.path.join(project_dir, "analysis/input/fastq")
clean_dir = os.path.join(project_dir, "analysis/output/clean")

os.makedirs(clean_dir, exist_ok=True)
   
def process_file(infile):
    filename = os.path.basename(infile)
    sample = filename.replace('.psl', '')
    fastq_file = os.path.join(fastq_dir, "{}.fastq.gz".format(sample))
    df = pd.read_csv(infile, sep="\t", header=None, usecols=[9, 13, 15, 16], skiprows=5)
    df.columns = ['counts_sequence', 'gene', 'start', 'end']
   
    df[['sequence', 'counts']] = df.get('counts_sequence').str.split('_', expand=True)
    sequences = df['sequence'].unique()
    
    fastq_df = pd.read_csv(fastq_file, sep='\n', header=None)
    idx0 = range(0, len(fastq_df), 4)
    idx1 = range(1, len(fastq_df), 4)
    idx2 = range(2, len(fastq_df), 4)
    idx3 = range(3, len(fastq_df), 4)

    col0 = fastq_df.iloc[idx0]
    col1 = fastq_df.iloc[idx1]
    col2 = fastq_df.iloc[idx2]
    col3 = fastq_df.iloc[idx3]

    df1 = pd.DataFrame()
    df1['col0'] = col0[0].reset_index()[0]
    df1['col1'] = col1[0].reset_index()[0]
    df1['col2'] = col2[0].reset_index()[0]
    df1['col3'] = col3[0].reset_index()[0]

    fastq_df = df1

    # remove rrna reads
    fastq_df = fastq_df.loc[~fastq_df['col1'].isin(sequences)]

    clean_file = os.path.join(clean_dir, '{}.fastq.gz'.format(sample))
    print('Writing file: {}'.format(clean_file))
    fastq_df.to_csv(clean_file, sep="\n", header=False, index=False, compression='gzip')    



if sample_file is None:
    blat_results = sorted(glob.glob(os.path.join(input_dir, "*.psl")))
    for infile in blat_results:
        print('bsub -q long -R "rusage[mem={}]" python {} {} {}'.format(mem, sys.argv[0], project_id, infile))
else:
    process_file(sample_file)
 
