#!/usr/bin/env python
import pandas as pd
import glob
import os
import sys
import matplotlib.pyplot as plt
import itertools
import numpy as np

project_id = sys.argv[1]
bam_type = sys.argv[2]
aa = sys.argv[3] # e.g. pro, asp, leu ...

x = 1
if len(sys.argv) >= 5:
    x = int(sys.argv[4])

min_reads = 50 #
if len(sys.argv) >= 6:
    min_reads = int(sys.argv[5])

genes = []
if len(sys.argv) >= 7:
    if sys.argv[6] in ['all', 'all_genes']:
        genes = []
    elif os.path.isfile(sys.argv[6]):
        subset = os.path.basename(sys.argv[6]).replace('.tsv', '')
        genes = pd.read_csv(sys.argv[6], sep="\t", usecols=[0], header=None)[0].tolist()
    else:
        subset = sys.argv[6]
        genes = [sys.argv[6]]

x_axis = None
if len(sys.argv) >= 8:
    x_axis = sys.argv[7]
    if x_axis not in ['log2', 'log10']:
        print('Warning: Unknown x_axis: {}. Can be [None, log2, log10]. Will use default value: None')
        x_axis = None

codons_file = '/icgc/dkfzlsdf/analysis/OE0532/static/hg19/codons.txt'
codons_df = pd.read_csv(codons_file, sep="\t", header=None, usecols=[0,1])
codons_df = codons_df.loc[codons_df[1].str.lower() == aa]
codons = codons_df[0].tolist()
codons = codons + codons
codons = list(itertools.permutations(codons, x))
codons = [''.join(c) for c in codons]
codons = list(set(codons))

fasta_file = '/icgc/dkfzlsdf/analysis/OE0532/static/hg19/transcripts_cds.fa'
fasta_df = pd.read_csv(fasta_file, sep=">", header=None)
trans = fasta_df[1].dropna().tolist()
seq = fasta_df[0].dropna().tolist()
fasta_df = pd.DataFrame(columns=['trans', 'seq'])
fasta_df['trans'] = trans
fasta_df['seq'] = seq

codon_file = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/codons/{}x/{}_codon_positions.tsv".format(x, aa)
c_df = pd.read_csv(codon_file, sep="\t", header=None, names=['trans', 'codon', 'pos'])
c_df = c_df.sort_values(['trans', 'pos'])
# keep only first codon
c_df = c_df.drop_duplicates('trans')

trans_file = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/transcriptLength.txt"
trans_df = pd.read_csv(trans_file, sep="\t")

# select genes of interest
if len(genes) != 0:
    trans_df = trans_df.loc[trans_df['gene_name'].isin(genes)]

trans_df = trans_df[['transcript', 'cds_len', '5utr_len']]
trans_df.columns = ['trans', 'cds_len', '5utr_len']
trans = trans_df['trans'].tolist()

indir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/tsv".format(project_id, bam_type)
if len(genes) >= 1:
    plot_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}//analysis/output/figures/ext_diricore/{}_cumulative/{}/{}".format(project_id, aa, bam_type, subset) 
else:
    plot_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}//analysis/output/figures/ext_diricore/{}_cumulative/{}".format(project_id, aa, bam_type)

os.makedirs(plot_dir, exist_ok=True)

for f in glob.glob("{}/*.tsv".format(indir)):
    sample = os.path.basename(f).replace('.tsv', '')
    df = pd.read_csv(f, sep="\t", header=None, names=['trans', 'start', 'seq'])
    df['start'] = df['start'].astype(int)
    # convert 1-based coordinates into 0-based coordinates
    df['start'] = df['start'] - 1
    df = df.loc[df['trans'].isin(trans)]
    # convert coordinates. 0 is the beginning of cds. 
    df = pd.merge(df, trans_df, on='trans', how='inner')
    df['start'] = df['start'] - df['5utr_len']
    # drop reads which aligned before the pos of the double codon
    df = pd.merge(df, c_df, on='trans')
    df = df.loc[df['start'] >= df['pos']]
    # drop reads from 3' UTR
    df = df.loc[df['start'] < df['cds_len']]
    # normalize by cds length
    df['start'] = df['start'] - df['pos']
    df['start'] = df['start'] / (df['cds_len'] - df['pos'])
    # check how many reads we have in each transcript
    df['counts'] = 1
    counts = df.groupby('trans')['counts'].sum().reset_index()
    df = pd.merge(counts, df.drop('counts', axis=1), on='trans')
    df = df.loc[df['counts'] >= min_reads] 
    # 
    df['start'] = df['start'].round(2)
    df = df[['trans', 'counts', 'start']]
    df = df.sort_values(['trans', 'start'])
    df['counts'] = 1 / df['counts'] # normalize to total counts per transcript
    # group by trans & start and count reads in each group
    df = df.groupby(['trans', 'start'])['counts'].sum().reset_index()
    trans_num = len(df['trans'].unique())
    df = df.groupby('start')['counts'].sum().reset_index()
    df['counts'] = df['counts'] / trans_num
    # get cumulative counts
    df['counts'] = df['counts'].cumsum()
    if x_axis == 'log2':
        df = df.loc[df['start'] != 0]
        df['start'] = np.log2(df['start'])
    # df['cum_sum'] = df.groupby(['trans'])['counts'].cumsum()
#    plot_df = df.groupby('start')['cum_sum'].mean().reset_index()
#    plot_df['cum_sum'] = plot_df['cum_sum'].round(2)
#    plot_df.columns = ['transcript', sample]
    df.columns = ['transcript', sample]
    if len(df) == 0:
        print('No data to plot for sample {}'.format(sample))
        continue
    plot_file = "{}/{}.pdf".format(plot_dir, sample)
    fig = df.plot(kind='line', x='transcript', y=sample)
    fig.set_title('Sample: {}. Amino Acid: {}. Min reads: {}'.format(sample, aa, min_reads))
    fig.set_ylabel('Percent reads')
    if x_axis == 'log2':
        fig.set_xlabel('log2(transcript)')
    plt.savefig(plot_file)
    print('Saved figure: {}'.format(plot_file))
