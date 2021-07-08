#!/usr/bin/env python

import sys
import os
import glob
import pandas as pd

project_id = sys.argv[1]

BASE_DIR = os.getenv('BASE_DIR')
DIRICORE_PATH = "$BASE_DIR/software/diricore"

project_dir = os.path.join(BASE_DIR, project_id)
input_dir = os.path.join(project_dir, "analysis/output/rrna/blat_results")
gene_outdir = os.path.join(project_dir, "analysis/output/rrna/reads_per_gene")
pos_outdir = os.path.join(project_dir, 'analysis/output/rrna/reads_per_position')
align_dir = os.path.join(project_dir, 'analysis/output/alignment_stats')
os.makedirs(gene_outdir, exist_ok=True)
os.makedirs(pos_outdir, exist_ok=True)
os.makedirs(align_dir, exist_ok=True)

blat_results = glob.glob(os.path.join(input_dir, "*.psl"))
stats_df = pd.DataFrame(columns=['sample', 'rrna_reads'])
for infile in blat_results:
    filename = os.path.basename(infile)
    print("Processing ", infile)
    gene_outfile = os.path.join(gene_outdir, filename.replace('.psl', '.tsv'))
    pos_outfile = os.path.join(pos_outdir, filename.replace('.psl', '.tsv'))
    # a line looks like: 345202_TTCGCGCGGGTCGGGGGGCGGGGCGGACTGT 18 8 25 31 100.0 chr6_cox_hap2 - 4700634 4700651 18 None None
    # we only need our id (which is counts+sequence), then gene, start and end positions
    df = pd.read_csv(infile, sep="\t", header=None, usecols=[9, 13, 15, 16], skiprows=5)
    df.columns = ['counts_sequence', 'gene', 'start', 'end']
    if df.empty:
         print("No data loaded! Skipping sample")
         continue
    # in hg19 rRNA genes are named in a wierd way: hg19_wgEncodeGencodeBasicV19_ENST00000607521.1__RNA28S5
    # so, let's get rid of the left part
    df['gene'] = df['gene'].str.split('__').str[-1] # in case genes are named in a normal way, the name will not be modified
    df = df.drop_duplicates('counts_sequence') # in case of multi-alignment, let's keep only the first alignment (the first gene which was in fasta)
    df[['sequence', 'counts']] = df.get('counts_sequence').str.split('_', expand=True)
    df['counts'] = df['counts'].astype(int)
    df1 = df.groupby('gene').agg({'counts': 'sum'}).reset_index()
    df1 = df1[['counts', 'gene']]
    df1.columns = ['gene_counts', 'gene_name']
    print('Writing file: {}'.format(gene_outfile))
    df1.to_csv(gene_outfile, sep="\t", index=False)

    sample = filename.replace('.psl', '')
    stats_df = stats_df.append({'sample': sample, 'rrna_reads': int(df1['gene_counts'].sum())}, ignore_index=True)
    df['length'] = df['sequence'].str.len()
    df2 = df.groupby(['gene', 'start', 'length']).agg({'counts': 'sum'}).reset_index()
    df2['reads_info'] = df2['counts'].astype(str) + ' reads of length ' + df2['length'].astype(str) + ', '
    df2 = df2.groupby(['gene', 'start']).agg({'reads_info': 'sum', 'counts': 'sum'}).reset_index()
    df2['reads_info'] = df2['reads_info'].str[:-2] # remove the last comma
    df2 = df2[['start', 'counts', 'gene', 'reads_info']]
    print('Writing file: {}'.format(pos_outfile))
    df2.to_csv(pos_outfile, sep="\t", index=False)

stats_df['rrna_reads'] = stats_df['rrna_reads'].astype(int)
stats_outfile = os.path.join(project_dir, 'analysis/output/alignment_stats/rrna_stats.txt')
print('Writing stats file: {}'.format(stats_outfile))
stats_df.to_csv(stats_outfile, sep="\t", index=False)
