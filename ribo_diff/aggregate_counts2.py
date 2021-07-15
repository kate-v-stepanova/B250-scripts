import pandas as pd
import os
import sys
import random
import glob

BASE_DIR = os.getenv('BASE_DIR')

rna_id = sys.argv[1]
rp_id = sys.argv[2]

bam_type = sys.argv[3]

sample1 = sys.argv[4]
sample2 = sys.argv[5]

# whatever is in the 5th argument, we include non-coding genes too
coding_only = True
noncoding = ['non_coding', 'non-coding', 'noncoding']
if len(sys.argv) >= 7:
    if sys.argv[6] in noncoding :
        coding_only = False
prefix = "coding" if coding_only else "all"


min_reads = 100
if len(sys.argv) >= 8:
    min_reads = int(sys.argv[7])

outfile_rna = "{}/{}/analysis/output/ribo_diff/reads/{}/{}_vs_{}_{}.tsv".format(BASE_DIR, rna_id, bam_type, sample1, sample2, prefix)
outfile_rp = "{}/{}/analysis/output/ribo_diff/reads/{}/{}_vs_{}_{}.tsv".format(BASE_DIR, rp_id, bam_type, sample1, sample2, prefix)

# create output dirs
os.makedirs("{}/{}/analysis/output/ribo_diff/reads/{}".format(BASE_DIR, rna_id, bam_type), exist_ok=True)
os.makedirs("{}/{}/analysis/output/ribo_diff/reads/{}".format(BASE_DIR, rp_id, bam_type), exist_ok=True)

indir_rna = "{}/{}/analysis/output/ext_diricore/{}/tsv".format(BASE_DIR, rna_id, bam_type)
indir_rp = "{}/{}/analysis/output/ext_diricore/{}/tsv".format(BASE_DIR, rp_id, bam_type)

gene_names = "{}/static/hg19/transcriptLength.txt".format(BASE_DIR)
coding_genes_file = "{}/static/hg19/coding_genes.txt".format(BASE_DIR)

sno_rnas_file ="{}/static/hg19/snoRNAs/snoRNAs.txt"
sno_rnas = pd.read_csv(sno_rnas_file, sep="|", header=None, names=['gene_id', 'gene_name'])


names_df = pd.read_csv(gene_names, sep="\t", header=None, names=["trans", "gene_id", "gene_name"], usecols=[0,1,2], skiprows=1)
names_df = names_df.drop_duplicates()

with open(coding_genes_file, 'r') as fff:
    coding_genes = fff.read().splitlines()

samples = [sample1, sample2]

print("Samples: {}".format(",".join(samples)))

# aggregate
df = None
for sample in samples:
    rna_files = sorted(glob.glob(os.path.join(indir_rna, '{}*.tsv'.format(sample))))
    rp_files = sorted(glob.glob(os.path.join(indir_rp, '{}*.tsv'.format(sample))))
    if len(rna_files) != len(rp_files):
        print('Files do not match!')
        print('RP directory: {}'.format(indir_rp))
        print('RNA directory: {}'.format(indir_rna))
        sys.exit(1)
    for rna_file in rna_files:
        samplename = os.path.basename(rna_file).replace('.tsv', '')
        rna_df = pd.read_csv(rna_file, skiprows=4, names=["trans", "RNA_{}".format(samplename)], usecols=[0,1], sep="\t")
        rna_df['RNA_{}'.format(samplename)] = 1
        rna_df = rna_df.groupby('trans')['RNA_{}'.format(samplename)].sum().reset_index()
        rna_df = pd.merge(rna_df, names_df, how='inner')
        rna_df = rna_df.drop('trans', axis=1)
        rna_df = rna_df.groupby(['gene_id', 'gene_name'])['RNA_{}'.format(samplename)].sum().reset_index()
        if df is None:
            df = rna_df
        else:
            df = pd.merge(df, rna_df, on=['gene_id', 'gene_name'], how='outer')
    for rp_file in rp_files:
        samplename = os.path.basename(rp_file).replace('.tsv', '')
        rp_df = pd.read_csv(rp_file, skiprows=4, names=["trans", "RP_{}".format(sample)], usecols=[0,1], sep="\t")
        rp_df['RP_{}'.format(samplename)] = 1
        rp_df = rp_df.groupby('trans')['RP_{}'.format(samplename)].sum().reset_index()
        rp_df = pd.merge(rp_df, names_df, how='inner')
        rp_df = rp_df.drop('trans', axis=1)
        rp_df = rp_df.groupby(['gene_id', 'gene_name'])['RP_{}'.format(samplename)].sum().reset_index()
        if df is None:
            df = rp_df
        else:
            df = pd.merge(df, rp_df, on=['gene_id', 'gene_name'], how='outer')

df = df.fillna(0)

# filter 
# remove if all zeros
cols = list(df.columns)
cols.remove("gene_id")
#df = df.loc[(df[cols]>=20).any(1)]
cols = list(filter(lambda x: x.startswith("RP"), cols))
df = df.loc[(df[cols]>=min_reads).any(1)]

#cols = list(df.columns)
#cols.remove("gene_id")
#for col in cols:
#   new_col1 = "R1_" + str(col)
#   new_col2 = "R2_" + str(col)
#   df[new_col1] = df[col].apply(split)
#   df[new_col2] = df[col].sub(df[new_col1], axis=0)
#   df = df.drop(columns=[col])

#if coding_only:
#    df = df.loc[df['gene_id'].isin(coding_genes)]

# fill na with gene_id
df['gene_name'] = df['gene_name'].fillna(df['gene_id'])
df = df.loc[~df['gene_name'].str.startswith('MT-')]
df = df.loc[~df['gene_name'].str.startswith('HIST')]
df = df.loc[~df['gene_name'].isin(sno_rnas['gene_name'].tolist())]

df['Entry'] = df['gene_name']
df = df.drop(columns=["gene_name", "gene_id"])  

# reordering
cols = list(df.columns)
cols.remove("Entry")
df[cols] = df[cols].astype(int)
cols = ["Entry"] + cols
df = df[cols]

print("Writing file: {}".format(outfile_rna))
print("Writing file: {}".format(outfile_rp))
df.to_csv(outfile_rna, sep="\t", header=True, index=False)
df.to_csv(outfile_rp, sep="\t", index=False)


