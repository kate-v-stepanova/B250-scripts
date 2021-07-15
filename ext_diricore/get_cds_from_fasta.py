import pandas as pd
import glob
import os
import sys

project_id = sys.argv[1]
genome = sys.argv[2]
minreads = sys.argv[3]
bam_type = sys.argv[4]

BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = os.path.join(BASE_DIR, project_id)
indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/tsv".format(bam_type))
transcripts = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/transcriptLength.txt".format(genome)
trans_df = pd.read_csv(transcripts, sep="\t")
cds_file = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/transcript_cds.fa"


trans_fasta = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/transcripts_single_header.fa".format(genome)
trans_seq = pd.read_csv(trans_fasta, sep=">", header=None)
seq = trans_seq[0].dropna()
trans = trans_seq[1].dropna()
trans_seq = pd.DataFrame()
trans_seq['transcript'] = trans.tolist()
trans_seq['seq'] = seq.tolist()
trans_df = pd.merge(trans_df, trans_seq, on="transcript", how="inner")

def trim_utrs(row):
    # trim 5' UTR
    cds = row['seq'][row['5utr_len']:]
    # trim 3' UTR
    cds = row['seq'][:row['cds_len']]
    return cds

trans_df['cds'] = trans_df.apply(trim_utrs, axis=1)
cds_df = trans_df[['transcript', 'cds']]
cds_df['transcript'] = ">" + cds_df['transcript']
cds_df.to_csv(cds_file, sep="\n", header=False, index=False)

for f in glob.glob(os.path.join(indir, '*.tsv')):
    df = pd.read_csv(f, sep="\t", header=None, names=['transcript', 'start', 'seq'])
    import pdb; pdb.set_trace()
    
