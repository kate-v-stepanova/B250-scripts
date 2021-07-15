import pandas as pd
import sys
import os

trans = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/transcriptLength.txt"
trans_df = pd.read_csv(trans, sep="\t")

hg19_fasta = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/transcripts.fa"
fasta_df = pd.read_csv(hg19_fasta, sep=">", header=None)
seq = fasta_df[0].dropna()
trans_id = fasta_df[1].dropna()
trans_id = trans_id.str.split('|').str[0]
fasta_df = pd.DataFrame(columns=['transcript', 'seq'])
fasta_df['seq'] = seq.tolist()
fasta_df['transcript'] = trans_id.tolist()
df = pd.merge(fasta_df, trans_df, on='transcript')
# cds = seq[5utr_len:]. axis=1 means rows
df['cds'] = df.apply(lambda row: row['seq'][row['5utr_len']:], axis=1)
# drop transcripts which don't start with ATG -> only 114 out of 78K transcripts. Probably annotation is wrong, or something about they are non-coding.
df = df.loc[df['cds'].str.startswith('ATG')]
df['transcript'] = '>' +  df['transcript']
outfile = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/transcripts_cds.fa"
print('Writing file: {}'.format(outfile))
df[['transcript', 'cds']].to_csv(outfile, sep="\n", header=False, index=False)
