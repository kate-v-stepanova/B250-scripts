#!/usr/bin/env python
import pandas as pd
import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt

plt.autoscale()

project_id = sys.argv[1]
bam_type = sys.argv[2]
aa = sys.argv[3]
fc = float(sys.argv[4])

contrasts = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/input/metadata/rpf_density_contrasts.tsv".format(project_id)
contrasts_df = pd.read_csv(contrasts, sep="\t", header=None, usecols=[0,1], names=['sample', 'control'])

indir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites_dotplot/{}".format(project_id, bam_type, aa)
p_outdir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites_dotplot_outliers/psite/{}".format(project_id, bam_type, aa)
a_outdir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites_dotplot_outliers/asite/{}".format(project_id, bam_type, aa)
e_outdir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites_dotplot_outliers/esite/{}".format(project_id, bam_type, aa)

os.makedirs(p_outdir, exist_ok=True)
os.makedirs(e_outdir, exist_ok=True)
os.makedirs(a_outdir, exist_ok=True)

p_df = pd.read_csv('{}/{}_psite_dotplot.tsv'.format(indir, aa), sep="\t")
a_df = pd.read_csv('{}/{}_asite_dotplot.tsv'.format(indir, aa), sep="\t")
e_df = pd.read_csv('{}/{}_esite_dotplot.tsv'.format(indir, aa), sep="\t")

for i, row in contrasts_df.iterrows():
    sample = row['sample']
    control = row['control']
    contrast = '{}__vs__{}'.format(sample, control)
    # psites
    df = p_df[['gene', 'rpkm_{}'.format(sample), 'rpkm_{}'.format(control)]]
    df = df.dropna()
    df['log2(fc)'] = np.log2(df['rpkm_' + sample] / df['rpkm_' + control])
    df = df.loc[df['log2(fc)'].abs() >= fc]
    p_outfile = "{}/{}.tsv".format(p_outdir, contrast)
    print('Writing file: {}'.format(p_outfile))
    df.to_csv(p_outfile, sep="\t", index=False)
    # asites
    df = a_df[['gene', 'rpkm_{}'.format(sample), 'rpkm_{}'.format(control)]]
    df = df.dropna()
    df[contrast] = np.log2(df['rpkm_' + sample] / df['rpkm_' + control])
    df = df.loc[df[contrast] >= fc]
    a_outfile = "{}/{}.tsv".format(a_outdir, contrast)
    print('Writing file: {}'.format(a_outfile))
    df.to_csv(a_outfile, sep="\t", index=False)
    # esites
    df = e_df[['gene', 'rpkm_{}'.format(sample), 'rpkm_{}'.format(control)]]
    df = df.dropna()
    df[contrast] = np.log2(df['rpkm_' + sample] / df['rpkm_' + control])
    df = df.loc[df[contrast].abs() >= fc]
    e_outfile = "{}/{}.tsv".format(e_outdir, contrast)
    print('Writing file: {}'.format(e_outfile))
    df.to_csv(e_outfile, sep="\t", index=False)
    
