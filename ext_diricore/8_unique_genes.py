import pandas as pd
import os
import sys
import glob

project_id = sys.argv[1]
bam_type = sys.argv[2]

indir = '/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/tsv/'.format(project_id, bam_type)
outdir = '/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/unique_genes/{}'.format(project_id, bam_type)

os.makedirs(outdir, exist_ok=True)
genome = 'hg19'
if len(sys.argv) >= 4:
    genome = sys.argv[3]

threshold = 100

sample = None
if len(sys.argv) >= 5:
    sample = sys.argv[4]

if sample is None:
    for f in sorted(glob.glob('{}/*.tsv'.format(indir))):
        sample = os.path.basename(f).replace('.tsv', '')
        print('python {} {}'.format(' '.join(sys.argv), sample))
    exit(0)

f = "{}/{}.tsv".format(indir, sample)
df = pd.read_csv(f, sep="\t", header=None, names=['trans', 'pos', 'seq'])

trans_file = '/icgc/dkfzlsdf/analysis/OE0532/static/{}/gene_names.txt'.format(genome)
t_df = pd.read_csv(trans_file, sep="|", header=None, names=['trans', 'gene'], usecols=[0,2])

df = pd.merge(df, t_df, on='trans')
total = len(df)
df['counts'] = 1
df = df.groupby('gene')['counts'].sum().reset_index()
df = df.loc[df['counts'] >= threshold]
df['norm'] = df['counts'] / total

outfile = '{}/{}.tsv'.format(outdir, sample)
print('Writing file: {}'.format(outfile))
df.to_csv(outfile, sep="\t", index=False)
