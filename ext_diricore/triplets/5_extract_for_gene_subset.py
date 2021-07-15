import pandas as pd
import glob
import os
import sys
import matplotlib.pyplot as plt

project_id = sys.argv[1]
bam_type = sys.argv[2]
x = int(sys.argv[3])
site = int(sys.argv[4])

window = 0
if len(sys.argv) >= 6:
    window = int(sys.argv[5])

genes = []
if len(sys.argv) >= 7:
    gene = sys.argv[6]
    if sys.argv[6] == "all":
        genes = []
        subset_name = sys.argv[7]
    elif os.path.isfile(gene):
         print('Reading file: {}'.format(gene))
         subset_name = sys.argv[7]
         subset_name = os.path.basename(gene).replace('.tsv', '').replace('.txt', '')
         genes = pd.read_csv(gene, sep="\t", header=None)
         genes = genes[0].tolist()
    else:
         subset_name = sys.argv[7]
         genes = [gene]

if len(sys.argv) >= 9:
   aa = sys.argv[8].split(',')


BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = os.path.join(BASE_DIR, project_id)
indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/codons/{}x.window{}".format(bam_type, x, window))
outdir = os.path.join(indir, subset_name)

infile = os.path.join(indir, "all_samples.{}.tsv".format(site))

df = pd.read_csv(infile, sep="\t")
if len(genes) != 0:
    df = df.loc[df['gene'].isin(genes)]
df = df.loc[df['aa'].isin(aa)]
import pdb; pdb.set_trace()
outfile = os.path.join(outdir, "all_samples.{}.tsv".format(site))
print('Writing {}'.format(outfile))
df.to_csv(outfile, sep="\t", index=False)

