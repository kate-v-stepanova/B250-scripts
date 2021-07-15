import pandas as pd
import os
import sys
import glob

project_id = sys.argv[1]

bam_type = "all_unique"
if len(sys.argv) >= 3:
    bam_type = sys.argv[2]

genome = 'hg19'
if len(sys.argv) >= 4:
    genome = sys.argv[3]

# TODO:
genes = []
if len(sys.argv) >= 5:
    gene_file = sys.argv[4]
    gene_df = pd.read_csv(gene_file, sep="\t", usecols=[0])
    genes = gene_df[0].tolist()

BASE_DIR = os.getenv('BASE_DIR')

bam_dir = "{}/{}/analysis/output/ext_diricore/{}/tsv/".format(BASE_DIR, project_id, bam_type)
output_dir = "{}/{}/analysis/output/alignments/reads_per_gene/normalized/{}".format(BASE_DIR, project_id, bam_type)

trans_len = "{}/static/{}/transcriptLength.txt".format(BASE_DIR, genome)
trans_df = pd.read_csv(trans_len, sep="\t")
trans_df = trans_df[['transcript', 'gene_name', 'gene_id', 'length']]
os.makedirs(output_dir, exist_ok=True)
for f in glob.glob(os.path.join(bam_dir, '*.tsv')):
    samplename = os.path.basename(f).replace('.tsv', '')
    df = pd.read_csv(f, sep="\t", header=None, names=['transcript', 'start', 'seq'])
    df = df.drop(['start', 'seq'], axis=1)
    df['counts'] = 1
    df = df.groupby('transcript').sum().reset_index()
    # df = df.loc[df['counts'] >= min_reads]
    df = pd.merge(df, trans_df, on='transcript', how='left')
    df = df.dropna()

    # tpm
    df['rpk'] = df['counts'] / df['length']
    per_million_factor = df['rpk'].sum() / 1000000
    df['tpm'] = df['rpk'] / per_million_factor 

    # rpkm
    per_million_factor = df['counts'].sum() / 1000000
    df['rpm'] = df['counts'] / per_million_factor
    df['rpkm'] = df['rpm'] / df['length']

    # cpm
    df['cpm'] = df['counts'] / (df['counts'].sum() * 1000000)

    norm_df = df[['transcript', 'gene_name', 'gene_id', 'counts', 'cpm', 'tpm', 'rpkm']]
    outfile = os.path.join(output_dir, "{}.tsv".format(samplename))
    print('Writing file: {}'.format(outfile))
    norm_df.to_csv(outfile, sep="\t", index=False)
