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

contrasts = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/input/metadata/rpf_density_contrasts.tsv".format(project_id)
contrasts_df = pd.read_csv(contrasts, sep="\t", header=None, usecols=[0,1], names=['sample', 'control'])

dotplot_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites_dotplot".format(project_id, bam_type)
outdir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites_dotplot/{}".format(project_id, bam_type, aa)
p_plot_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/figures/ext_diricore/{}/psites_dotplot/{}".format(project_id, bam_type, aa)
a_plot_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/figures/ext_diricore/{}/asites_dotplot/{}".format(project_id, bam_type, aa)
e_plot_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/figures/ext_diricore/{}/esites_dotplot/{}".format(project_id, bam_type, aa)
os.makedirs(p_plot_dir, exist_ok=True)
os.makedirs(a_plot_dir, exist_ok=True)
os.makedirs(e_plot_dir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)

# psites
p_df = None
for f in sorted(glob.glob("{}/psite_*.tsv".format(dotplot_dir))):
    sample = os.path.basename(f).replace('.tsv', '').replace('psite_', '')
    df = pd.read_csv(f, sep="\t")
    df = df.loc[df['Aa'].str.lower() == aa]
    df = df.drop('cpm_counts', axis=1)
    df.columns = ['gene', 'Aa', 'codon', sample, 'tpm_{}'.format(sample), 'rpkm_{}'.format(sample)]
    if p_df is None:
        p_df = df
    else:
        p_df = pd.merge(p_df, df, on=['gene', 'Aa', 'codon'], how='outer')

# asites
a_df = None
for f in sorted(glob.glob("{}/asite_*.tsv".format(dotplot_dir))):
    sample = os.path.basename(f).replace('.tsv', '').replace('asite_', '')
    df = pd.read_csv(f, sep="\t")
    df = df.loc[df['Aa'].str.lower() == aa]
    df = df.drop('cpm_counts', axis=1)
    df.columns = ['gene', 'Aa', 'codon', sample, 'tpm_{}'.format(sample), 'rpkm_{}'.format(sample)]
    if a_df is None:
        a_df = df
    else:
        a_df = pd.merge(a_df, df, on=['gene', 'Aa', 'codon'], how='outer')

# esites
e_df = None
for f in sorted(glob.glob("{}/esite_*.tsv".format(dotplot_dir))):
    sample = os.path.basename(f).replace('.tsv', '').replace('esite_', '')
    df = pd.read_csv(f, sep="\t")
    df = df.loc[df['Aa'].str.lower() == aa]
    df = df.drop('cpm_counts', axis=1)
    df.columns = ['gene', 'Aa', 'codon', sample, 'tpm_{}'.format(sample), 'rpkm_{}'.format(sample)]
    if e_df is None:
        e_df = df
    else:
        e_df = pd.merge(e_df, df, on=['gene', 'Aa', 'codon'], how='outer')

p_outfile = "{}/{}_psite_dotplot.tsv".format(outdir, aa)
a_outfile = "{}/{}_asite_dotplot.tsv".format(outdir, aa)
e_outfile = "{}/{}_esite_dotplot.tsv".format(outdir, aa)

print('Writing file: {}'.format(p_outfile))
p_df.to_csv(p_outfile, sep="\t", index=None)

print('Writing file: {}'.format(a_outfile))
a_df.to_csv(a_outfile, sep="\t", index=None)

print('Writing file: {}'.format(e_outfile))
e_df.to_csv(e_outfile, sep="\t", index=None)


# get the log2
for c in p_df.columns:
    if 'tpm' not in c and 'rpkm' not in c:
        continue
    p_df[c] = np.log2(p_df[c])
    e_df[c] = np.log2(e_df[c])
    a_df[c] = np.log2(a_df[c])

# axis limits  
p_tpm_min = p_df[[c for c in p_df.columns if 'tpm' in c]].dropna().min().min()
p_tpm_max = p_df[[c for c in p_df.columns if 'tpm' in c]].dropna().max().max()
a_tpm_min = a_df[[c for c in a_df.columns if 'tpm' in c]].dropna().min().min()
a_tpm_max = a_df[[c for c in a_df.columns if 'tpm' in c]].dropna().max().max()
e_tpm_min = e_df[[c for c in e_df.columns if 'tpm' in c]].dropna().min().min()
e_tpm_max = e_df[[c for c in e_df.columns if 'tpm' in c]].dropna().max().max()

# axis limits
p_rpkm_min = p_df[[c for c in p_df.columns if 'rpkm' in c]].dropna().min().min()
p_rpkm_max = p_df[[c for c in p_df.columns if 'rpkm' in c]].dropna().max().max()
a_rpkm_min = a_df[[c for c in a_df.columns if 'rpkm' in c]].dropna().min().min()
a_rpkm_max = a_df[[c for c in a_df.columns if 'rpkm' in c]].dropna().max().max()
e_rpkm_min = e_df[[c for c in e_df.columns if 'rpkm' in c]].dropna().min().min()
e_rpkm_max = e_df[[c for c in e_df.columns if 'rpkm' in c]].dropna().max().max()

for i, row in contrasts_df.iterrows():
    sample = row['sample']
    control = row['control']
    contrast = '{}__vs__{}'.format(sample, control)

    # P-site tpm
    tpm_plot = '{}/tpm__{}.pdf'.format(p_plot_dir, contrast)
    p_df1 = p_df[['tpm_{}'.format(sample), 'tpm_{}'.format(control)]]
    p_df1 = p_df1.dropna()
    p_df1.columns = [sample, control]
    fig1 = p_df1.plot(kind='scatter', x=sample, y=control, style='o', s=1)
    fig1.set_title('P-site signal in {}. TPM'.format(aa))
    fig1.set_xlabel('log2({})'.format(sample))
    fig1.set_ylabel('log2({})'.format(control))
    fig1.set_ylim(p_tpm_min, p_tpm_max)
    fig1.set_xlim(p_tpm_min, p_tpm_max)
    print('Writing file: {}'.format(tpm_plot))
    plt.savefig(tpm_plot)
    plt.clf()
    plt.close()

    # P-site rpkm
    rpkm_plot = '{}/rpkm__{}.pdf'.format(p_plot_dir, contrast)
    p_df2 = p_df[['rpkm_{}'.format(sample), 'rpkm_{}'.format(control)]]
    p_df2 = p_df2.dropna()
    fig2 = p_df2.plot(kind='scatter', x='rpkm_{}'.format(sample), y='rpkm_{}'.format(control), style='o', s=1)
    fig2.set_title('P-site signal in {}. RPKM'.format(aa))
    fig2.set_xlabel('log2({})'.format(sample))
    fig2.set_ylabel('log2({})'.format(control))
    fig2.set_ylim(p_rpkm_min, p_rpkm_max)
    fig2.set_xlim(p_rpkm_min, p_rpkm_max)
    print('Writing file: {}'.format(rpkm_plot))
    plt.savefig(rpkm_plot)
    plt.clf()
    plt.close()

    # A-site tpm
    tpm_plot = '{}/tpm__{}.pdf'.format(a_plot_dir, contrast)
    a_df1 = a_df[['tpm_{}'.format(sample), 'tpm_{}'.format(control)]]
    a_df1 = a_df1.dropna()
    fig1 = a_df1.plot(kind='scatter', x='tpm_{}'.format(sample), y='tpm_{}'.format(control), style='o', s=1)
    fig1.set_title('A-site signal in {}. TPM'.format(aa))
    fig1.set_xlabel('log2({})'.format(sample))
    fig1.set_ylabel('log2({})'.format(control))
    fig1.set_ylim(a_tpm_min, a_tpm_max)
    fig1.set_xlim(a_tpm_min, a_tpm_max)
 
    print('Writing file: {}'.format(tpm_plot))
    plt.savefig(tpm_plot)
    plt.clf()
    plt.close()

    # A-site rpkm
    rpkm_plot = '{}/rpkm__{}.pdf'.format(a_plot_dir, contrast)
    a_df2 = a_df[['rpkm_{}'.format(sample), 'rpkm_{}'.format(control)]]
    a_df2 = a_df2.dropna()
    fig2 = a_df2.plot(kind='scatter', x='rpkm_{}'.format(sample), y='rpkm_{}'.format(control), style='o', s=1)
    fig2.set_title('A-site signal in {}. RPKM'.format(aa))
    fig2.set_xlabel('log2({})'.format(sample))
    fig2.set_ylabel('log2({})'.format(control))
    fig2.set_ylim(a_rpkm_min, a_rpkm_max)
    fig2.set_xlim(a_rpkm_min, a_rpkm_max)
    print('Writing file: {}'.format(rpkm_plot))
    plt.savefig(rpkm_plot)
    plt.clf()
    plt.close()

    # E-site tpm
    tpm_plot = '{}/tpm__{}.pdf'.format(e_plot_dir, contrast)
    e_df1 = e_df[['tpm_{}'.format(sample), 'tpm_{}'.format(control)]]
    e_df1 = e_df1.dropna()
    fig1 = e_df1.plot(kind='scatter', x='tpm_{}'.format(sample), y='tpm_{}'.format(control), style='o', s=1)
    fig1.set_title('E-site signal in {}. TPM'.format(aa))
    fig1.set_xlabel('log2({})'.format(sample))
    fig1.set_ylabel('log2({})'.format(control))
    fig1.set_ylim(e_tpm_min, e_tpm_max)
    fig1.set_xlim(e_tpm_min, e_tpm_max)
    print('Writing file: {}'.format(tpm_plot))
    plt.savefig(tpm_plot)
    plt.clf()
    plt.close()



    # E-site rpkm
    rpkm_plot = '{}/rpkm__{}.pdf'.format(e_plot_dir, contrast)
    e_df2 = e_df[['rpkm_{}'.format(sample), 'rpkm_{}'.format(control)]]
    e_df2 = e_df2.dropna()
    fig2 = e_df2.plot(kind='scatter', x='rpkm_{}'.format(sample), y='rpkm_{}'.format(control), style='o', s=1)
    fig2.set_title('E-site signal in {}. RPKM'.format(aa))
    fig2.set_xlabel('log2({})'.format(sample))
    fig2.set_ylabel('log2({})'.format(control))
    fig2.set_ylim(e_rpkm_min, e_rpkm_max)
    fig2.set_xlim(e_rpkm_min, e_rpkm_max)
    print('Writing file: {}'.format(rpkm_plot))
    plt.savefig(rpkm_plot)
    plt.clf()
    plt.close()


