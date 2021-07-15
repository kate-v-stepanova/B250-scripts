import pandas as pd
import sys
import os
import glob
import matplotlib.pyplot as plt

BASE_DIR = os.getenv('BASE_DIR')

project_id = sys.argv[1]
bam_type = sys.argv[2]
genome = 'hg19'

if len(sys.argv) >= 4:
    genome = sys.argv[3]

genes = []
subset = 'all_genes'
if len(sys.argv) >= 5:
    gene_file = sys.argv[4] # can contain: gene names, gene ids, transcript ids
    if gene_file != 'all_genes':
        subset = os.path.basename(gene_file).replace('.txt', '').replace('.tsv', '')
        genes = pd.read_csv(gene_file, sep="\t", header=None)[0].tolist()

exclude = False
if len(sys.argv) >= 6:
    exclude = sys.argv[5] == 'exclude'


indir = "{}/{}/analysis/output/ext_diricore/{}/tsv/".format(BASE_DIR, project_id, bam_type)
plot_dir = "{}/{}/analysis/output/figures/size_dist/{}/{}".format(BASE_DIR, project_id, bam_type, subset)
if exclude:
    plot_dir = "{}/{}/analysis/output/figures/size_dist/{}/exclude_{}".format(BASE_DIR, project_id, bam_type, subset)

gene_names = "{}/static/{}/gene_names.txt".format(BASE_DIR, genome)
genes_df = pd.read_csv(gene_names, sep="|", header=None)
genes_df.columns = ['trans_id', 'gene_id', 'gene']

# if gene_file contains gene names
by_name = genes_df.loc[genes_df['gene'].isin(genes)]
if len(by_name) == len(genes) and len(genes) != 0:
   if exclude:
        genes_df = genes_df.loc[~genes_df['gene'].isin(genes)]
   else:
       genes_df = by_name

# if gene_file contains transcript ids
by_trans = genes_df.loc[genes_df['trans_id'].isin(genes)]
if len(by_trans) == len(genes) and len(genes) != 0:
    if exclude:
        genes_df = genes_df.loc[~genes_df['trans_id'].isin(genes)]
    else:
        genes_df = by_trans

# if gene file contrains gene ids
by_gene_id = genes_df.loc[genes_df['gene_id'].isin(genes)]
if len(by_gene_id) == len(genes) and len(genes) != 0:
    if exclude:
        genes_df = genes_df.loc[~genes_df['gene_id'].isin(genes)]
    else:
        genes_df = by_gene_id

# create output dirs
os.makedirs(plot_dir, exist_ok=True)

for f in glob.glob(os.path.join(indir, '*.tsv')):
    sample = os.path.basename(f).replace('.tsv', '')
    df = pd.read_csv(f, sep="\t", header=None, usecols=[0,2])
    df.columns = ['trans_id', 'seq']
    if len(genes_df) != 0:
        df = pd.merge(df, genes_df)
    df['length'] = df['seq'].str.len()
    df = df[['gene', 'length']]
    df = df.groupby('length').count().reset_index()
    df.columns = ['length', 'reads']
    # row counts
    outfile = os.path.join(plot_dir, os.path.basename(f))
    print('Writing file: {}'.format(outfile))
    df.to_csv(outfile, sep="\t", index=False)
    plot_file = outfile.replace('.tsv', '.pdf')
    print('Writing file: {}'.format(plot_file))
    title = '{} - {} - {}\nSize distribution {} (total reads: {})'.format(sample, project_id, bam_type, subset, df['reads'].sum())
    if exclude:
        title = '{} - {} - {}\nSize distribution. Exclude {}. Total reads: {}'.format(sample, project_id, bam_type, subset, df['reads'].sum())
    ax = df.plot(kind='bar', x='length', y='reads', legend=False, title=title)
    fig = ax.get_figure()
    fig.savefig(plot_file)
    plt.close(fig)

    # percentage
    df1 = pd.DataFrame()
    df1['length'] = df['length'].to_list()
    df1['percent'] = df['reads'].to_list()
    df1['percent'] = df1['percent'] / df['reads'].sum() * 100

    outfile1 = os.path.join(plot_dir, '{}_percent.tsv'.format(sample))
    print('Writing file: {}'.format(outfile1))
    df1.to_csv(outfile1, sep="\t", index=False)

    plot_file1 = os.path.join(plot_dir, '{}_percent.pdf'.format(sample))
    print('Writing file: {}'.format(plot_file1))
    ax = df1.plot(kind='bar', x='length', y='percent', legend=False, title=title)
    ax.set_ylabel("reads percentage")
    fig = ax.get_figure()
    fig.savefig(plot_file1)
    plt.close(fig)
    


