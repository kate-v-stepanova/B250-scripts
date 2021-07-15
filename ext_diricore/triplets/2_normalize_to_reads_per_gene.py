import pandas as pd
import glob
import os
import sys
import matplotlib.pyplot as plt

project_id = sys.argv[1]
genome = sys.argv[2]
minreads = int(sys.argv[3])
bam_type = sys.argv[4]

x = 1 # number of repetitions of codons. E.g. x=1 means we look for "ATG" and x=2 means we look for "ATGATG" and so on. 
# x can be theoretically any number, but for now it's intended to look either for duplets or triplets (x=2 and x=3)
if len(sys.argv) >= 6:
    x = int(sys.argv[5])

site = 12 # p_site can be: 15: a_site, 18: e_site
if len(sys.argv) >= 7:
    site = int(sys.argv[6])

window = 0
if len(sys.argv) >= 8:
    window = int(sys.argv[7])


BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = os.path.join(BASE_DIR, project_id)
indir1 = os.path.join(project_dir, "analysis/output/ext_diricore/{}/tsv".format(bam_type))
indir2 = os.path.join(project_dir, "analysis/output/ext_diricore/{}/codons/{}x.window{}".format(bam_type, x, window))
outdir = indir1

gene_names = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/gene_names.txt".format(genome)
genes_df = pd.read_csv(gene_names, sep="|", header=None, names=['transcript', 'gene_id', 'gene'])

for f in sorted(glob.glob(os.path.join(indir1, '*.tsv'))):
    sample = os.path.basename(f).replace('.tsv', '')
    f2 = os.path.join(indir2, "{}_{}.tsv".format(site, sample))
    full_df = pd.DataFrame()
    df = pd.read_csv(f, sep="\t", header=None, names=['transcript', 'start', 'seq'])
    df['counts_per_gene'] = 1

    df = df[['transcript', 'counts_per_gene']]
    df = df.groupby('transcript')['counts_per_gene'].sum().reset_index()
    df2 = pd.read_csv(f2, sep="\t") 
#    df2.columns = ['transcript', 'codon',  'counts',   'aa','norm_counts']
    import pdb; pdb.set_trace()
    df3 = pd.merge(df, df2, on="transcript", how='right')
    df3['norm_counts2'] = df3['counts'] / df3['counts_per_gene'] 
    print("Writing data file: {}".format(f2))
    df3.to_csv(f2, sep="\t", index=False)

