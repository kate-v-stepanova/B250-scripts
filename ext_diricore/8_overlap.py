import pandas as pd
import os
import sys
import glob

project_id = sys.argv[1]
bam_type = sys.argv[2]

s1 = sys.argv[3]
s2 = sys.argv[4]

indir = '/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/unique_genes/{}'.format(project_id, bam_type)
outdir = indir

f1 = "{}/{}.tsv".format(indir, s1)
f2 = "{}/{}.tsv".format(indir, s2)

df1 = pd.read_csv(f1, sep="\t")
df2 = pd.read_csv(f2, sep="\t") 

df1.columns = ['gene', s1, '{}_norm'.format(s1)]
df2.columns = ['gene', s2, '{}_norm'.format(s2)]

df = pd.merge(df1, df2, on='gene', how='inner')

outfile = '{}/{}__vs__{}.tsv'.format(outdir, s1, s2)
print('Writing file: {}'.format(outfile))
df.to_csv(outfile, sep="\t", index=False)

df1 = df1.loc[~df1['gene'].isin(df['gene'])]
outfile1 = '{}/{}_excl_overlap.tsv'.format(outdir, s1)
print('Writing file: {}'.format(outfile1))
df1.to_csv(outfile1, sep="\t", header=False, index=False)
