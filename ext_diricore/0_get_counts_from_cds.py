import sys
import os
import glob
import pandas as pd

BASE_DIR = os.getenv('BASE_DIR')

project_id = sys.argv[1]
bam_type = sys.argv[2]
genome = 'hg19'
if len(sys.argv) >= 4:
    genome = sys.argv[3]


indir = "{}/{}/analysis/output/ext_diricore/{}/tsv".format(BASE_DIR, project_id, bam_type) # {}/software/ext_diricore/1_get_seq_from_bam.sh 18927 all_unique
#indir = "{}/{}/analysis/output/alignments/reads_per_gene/sam/{}".format(BASE_DIR, project_id, bam_type)
outdir = "{}/{}/analysis/output/alignments/reads_per_gene/tsv/{}_cds".format(BASE_DIR, project_id, bam_type)
trans_file = "{}/static/{}/cds_info.tsv".format(BASE_DIR, genome)
trans_df = pd.read_csv(trans_file, sep="\t")
print(indir)
os.makedirs(outdir, exist_ok=True)
for f in glob.glob("{}/*.tsv".format(indir)):
    outfile = os.path.join(outdir, os.path.basename(f))
    if os.path.isfile(outfile):
        print('File exists. Skipping: {}'.format(outfile))
        continue
    df = pd.read_csv(f, sep="\t", header=None, names=['transcript', 'start', 'seq'], usecols=[0,1,2])
    df['tx_id'] = df['transcript'].str.split('.').str[0]
    df = pd.merge(df, trans_df, on='tx_id', how='inner')
    df = df.loc[(df['start'] >= df['cds_start']) & (df['start'] <= df['cds_end'])]
    # df = df.loc[df['read_strand'] == df['strand']] 
    df['counts'] = 1
    df1 = df.groupby('gene_name').counts.sum().reset_index()
    df = pd.merge(df1, df.drop('counts', axis=1), on='gene_name', how='left')
    df = df.drop_duplicates('gene_name')
    df = df[['gene_name', 'counts']]
    df.columns = ['gene', 'counts']
    print('Writing file: {}'.format(outfile))
    df.to_csv(outfile, sep="\t", index=False)

