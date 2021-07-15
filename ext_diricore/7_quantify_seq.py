import pandas as pd
import sys
import os
import glob
import matplotlib.pyplot as plt

# example: 16210
# quantify reads containing certain sequences

project_id = sys.argv[1]
bam_type = sys.argv[2]
seq_file = sys.argv[3]
seq_df = pd.read_csv(seq_file, sep="\t", header=None, names=['subset', 'seq'])
subset_name = os.path.basename(seq_file).replace('.txt', '').replace('.tsv', '')

base_dir = "/icgc/dkfzlsdf/analysis/OE0532/"
project_dir = os.path.join(base_dir, project_id)
indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/tsv".format(bam_type))
outdir = os.path.join(project_dir, "analysis/output/quantify_by_seq/{}".format(bam_type))
if not os.path.isdir(outdir):
    os.makedirs(outdir, exist_ok=True)

for seq in seq_df['seq'].tolist():
    subset = seq_df.loc[seq_df['seq'] == seq, 'subset'].iloc[0]
    for f in glob.glob(os.path.join(indir, '*.tsv')):
        sample = os.path.basename(f).replace('.tsv', '')
        if sample not in seq_df.columns:
            seq_df[sample] = 0
        df = pd.read_csv(f, sep="\t", header=None, names=['trans', 'start', 'seq'])
        df = df.loc[df['seq'].str.contains(seq)]
        seq_df.loc[seq_df['seq'] == seq, sample] = len(df)
        
outfile = os.path.join(outdir, "{}.tsv".format(subset_name))
print("Writing file: {}".format(outfile))
 
seq_df.to_csv(outfile, sep="\t", index=False)
