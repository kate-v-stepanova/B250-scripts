import pandas as pd
import sys
import glob
import os

BASE_DIR = os.getenv('BASE_DIR')

project_id = sys.argv[1]
bam_type = sys.argv[2]
norm = 'rpkm'
if len(sys.argv) >= 4:
    norm = sys.argv[3]
    if norm == 'raw':
        norm = 'counts'

min_reads = 1
if len(sys.argv) >= 5:
    min_reads = int(sys.argv[4])

genes = [] # can be a list of: transcript IDs, gene IDs, gene names.
subset = ''
if len(sys.argv) >= 6:
    gene_file = sys.argv[5]
    genes = pd.read_csv(gene_file, sep="\t", usecols=[0])
    genes.columns = ['gene']
    genes = genes['gene'].tolist()
    subset = os.path.basename(gene_file).replace('.tsv', '').replace('.txt', '') + '_'

groupby_key = 'gene_id'
if len(sys.argv) >= 7:
    groupby_key = sys.argv[6]

if groupby_key in ['transcript', 'trans', 'trans_id', 'transcript_id']:
   groupby_key = 'transcript'
if groupby_key in ['gene', 'gene_name']:
   groupby_key = 'gene_name'
else:
   groupby_key = 'gene_id'    


indir = os.path.join(BASE_DIR, project_id, "analysis/output/alignments/reads_per_gene/normalized//{}/".format(bam_type))

full_df = None
samples = []

full_df = None
infiles = sorted(glob.glob(os.path.join(indir, '*.tsv')))
if len(infiles) == 0:
    print('No files found in {}'.format(indir))
    exit()

raw_counts = None
for f in infiles:
    sample = os.path.basename(f).replace('.tsv', '')
    samples.append(sample)
    df = pd.read_csv(f, sep="\t")
    if norm not in df.columns:
        continue
    if len(genes) > 0:
        df1 = df.loc[df['transcript'].isin(genes)]
        if len(df1) == 0:
            df1 = df.loc[df['gene_name'].isin(genes)]
        if len(df1) == 0:
            df1 = df.loc[df['gene_id'].isin(genes)]
        df = df1
    import pdb; pdb.set_trace()
    if raw_counts is None:
        raw_counts = df
    else:
        raw_counts = pd.merge(raw_counts, df)
#    df = df.loc[df['counts'] >= min_reads]
    df = df[[groupby_key, norm]]
    df = df.groupby(groupby_key).agg({norm: 'sum'}).reset_index()
    df.columns = ['Gene', sample]
    df['Gene'] = df['Gene'].str.split('.').str[0]
    if full_df is None:
        full_df = df
    else:
        full_df = pd.merge(full_df, df, how='outer', on='Gene')

#full_df = full_df.fillna(0)
full_df = full_df.dropna()

outfile = os.path.join(indir, '{}{}_{}min_reads.tsv'.format(subset, norm, min_reads))

print('Writing file: {}'.format(outfile))


full_df.to_csv(outfile, sep="\t", index=False)
