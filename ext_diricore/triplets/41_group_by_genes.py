import pandas as pd
import glob
import sys
import os

project_id = sys.argv[1]
bam_type = sys.argv[2]
x = sys.argv[3]
site = sys.argv[4]
window = sys.argv[5]

indir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/codons/{}x.window{}".format(project_id, bam_type, x, window)
gene_df = pd.read_csv('/icgc/dkfzlsdf/analysis/OE0532/static/hg19/gene_names.txt', sep='|', header=None, names=['transcript', 'gene_id', 'gene_name'])

for f in glob.glob(os.path.join(indir, '{}*.tsv'.format(site))):
    df = pd.read_csv(f, sep="\t")
    import pdb; pdb.set_trace()
    df = pd.merge(df, gene_df, on='transcript', how='left')
    df = df[['gene_name','aa', 'codon',  'counts', 'norm_counts2']]
    df.columns = ['gene', 'aa', 'codon', 'counts', 'norm_counts']
    df = df.sort_values(by=['aa', 'gene', 'codon', 'norm_counts'])
    df = df.groupby(['aa', 'gene']).agg({'counts':'sum','norm_counts':'sum'}).reset_index()
    outfile = "gene_name.{}".format(os.path.basename(f))
    outfile = os.path.join(indir, outfile)
    print("Writing file {}".format(outfile))
    df.to_csv(outfile, sep="\t", index=False)

