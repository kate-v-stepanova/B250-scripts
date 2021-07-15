#!/usr/bin/env python
import pandas as pd
import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt


project_id = sys.argv[1]
bam_type = sys.argv[2]
norm_type = 'tpm'
if len(sys.argv) >= 4:
    norm_type = sys.argv[3]

indir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites".format(project_id, bam_type)
psite_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}//analysis/output/figures/ext_diricore/{}/psite_contrasts/{}".format(project_id, bam_type, norm_type)
asite_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}//analysis/output/figures/ext_diricore/{}/asite_contrasts/{}".format(project_id, bam_type, norm_type)
esite_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}//analysis/output/figures/ext_diricore/{}/esite_contrasts/{}".format(project_id, bam_type, norm_type)
os.makedirs(psite_dir, exist_ok=True)
os.makedirs(asite_dir, exist_ok=True)
os.makedirs(esite_dir, exist_ok=True)

psites_file = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites/psites.tsv".format(project_id, bam_type)
asites_file = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites/asites.tsv".format(project_id, bam_type)
esites_file = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites/esites.tsv".format(project_id, bam_type)
p_df = pd.read_csv(psites_file, sep="\t")
a_df = pd.read_csv(asites_file, sep="\t")
e_df = pd.read_csv(esites_file, sep="\t")

# drop stop codons
p_df = p_df.loc[p_df['aa'] != 'Stp']
a_df = a_df.loc[a_df['aa'] != 'Stp']
e_df = e_df.loc[e_df['aa'] != 'Stp']

contrast_file = '/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/input/metadata/rpf_density_contrasts.tsv'.format(project_id)
contrast_df = pd.read_csv(contrast_file, sep="\t", header=None, names=['sample', 'control', 'color'])
p_outfile = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites/p_barplot.tsv".format(project_id, bam_type)
e_outfile = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites/e_barplot.tsv".format(project_id, bam_type)
a_outfile = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites/a_barplot.tsv".format(project_id, bam_type)
p_barplot = None
e_barplot = None
a_barplot = None
for i, row in contrast_df.iterrows():
    control = row['control']
    sample = row['sample']
    contrast = '{}__vs__{}'.format(sample, control)
    # psite plots
    p_df1 = p_df[['codon', 'aa', '{}_'.format(norm_type) + sample, '{}_'.format(norm_type) + control]]
    p_df1.columns = ['codon', 'aa', sample, control]
    p_df1[contrast] = (p_df1[sample] - p_df1[control]) / p_df1[control]
    if p_barplot is None:
        p_barplot = p_df1[['codon', 'aa', contrast]]
    else:
        p_barplot = pd.merge(p_barplot, p_df1[['codon', 'aa', contrast]])
    ax = p_df1.plot.bar(x='codon', y=contrast, figsize=(20, 8))
    for i, p in enumerate(ax.patches):
        aa = p_df1.iloc[i]['aa']
        h = p.get_height()
        if h > 0:
            ax.text(p.get_x(), p.get_height() + 0.001, aa, va='bottom', rotation='vertical', fontsize=10, color='black')
        else:
            ax.text(p.get_x(), p.get_height() - 0.001, aa, rotation='vertical', va='top', fontsize=10, color='black')
    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon and 3) {}".format(norm_type.upper()))
    ax.set_ylabel('({} - {}) / {} *'.format(sample, control, control))
    fig = ax.get_figure()
    p_codon_file = '{}/{}_per_codon.pdf'.format(psite_dir, contrast)
    print('Writing file: {}'.format(p_codon_file))
    fig.savefig(p_codon_file)
    plt.cla()
    plt.close(fig)
    p_df1 = p_df1.groupby('aa').sum().reset_index()
    ax = p_df1.plot.bar(x='aa', y=contrast)
    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon and 3) {}".format(norm_type.upper()), fontsize=6)
    ax.set_ylabel('({} - {}) / {} *'.format(sample, control, control))
    fig = ax.get_figure()
    p_aa_file =  '{}/{}_per_aa.pdf'.format(psite_dir, contrast)
    print('Writing file: {}'.format(p_aa_file))
    fig.savefig(p_aa_file)
    plt.cla()
    plt.close(fig)
    # asite plots
    a_df1 = a_df[['codon', 'aa', '{}_'.format(norm_type) + sample, '{}_'.format(norm_type) + control]]
    a_df1.columns = ['codon', 'aa', sample, control]
    a_df1[contrast] = (a_df1[sample] - a_df1[control]) / a_df1[control]
#    a_df1 = a_df1.sort_values(by=sample)
    ax = a_df1.plot.bar(x='codon', y=contrast, figsize=(20, 8))
    for i, p in enumerate(ax.patches):
        aa = a_df1.iloc[i]['aa']
        h = p.get_height()
        if h > 0:
            ax.text(p.get_x(), p.get_height() + 0.001, aa, rotation='vertical', va='bottom', fontsize=10, color='black')
        else:
            ax.text(p.get_x(), p.get_height() - 0.001, aa, rotation='vertical', va='top', fontsize=10, color='black')
    if a_barplot is None:
        a_barplot = a_df1[['codon', 'aa', contrast]]
    else:
        a_barplot = pd.merge(a_barplot, a_df1[['codon', 'aa', contrast]])
 
    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon and 3) {}".format(norm_type.upper()))
    ax.set_ylabel('({} - {}) / {} *'.format(sample, control, control))

    fig = ax.get_figure()
    a_codon_file = '{}/{}_per_codon.pdf'.format(asite_dir, contrast)
    print('Writing file: {}'.format(a_codon_file))
    fig.savefig(a_codon_file)
    plt.cla()
    plt.close(fig)
 
    a_df1 = a_df1.groupby(['aa']).sum().reset_index()
    ax = a_df1.plot.bar(x='aa', y=contrast)
    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon and 3) {}".format(norm_type.upper()), fontsize=6)
    ax.set_ylabel('({} - {}) / {} *'.format(sample, control, control))
    fig = ax.get_figure()
    a_aa_file = '{}/{}_per_aa.pdf'.format(asite_dir, contrast)
    print('Writing file: {}'.format(a_aa_file))
    fig.savefig(a_aa_file)
    plt.cla()
    plt.close(fig)
 
    # esite plots
    e_df1 = e_df[['codon', 'aa', '{}_'.format(norm_type) + sample, '{}_'.format(norm_type) + control]]
    e_df1.columns = ['codon', 'aa', sample, control]
    e_df1[contrast] = (e_df1[sample] - e_df1[control]) / e_df1[control]
    ax = e_df1.plot.bar(x='codon', y=contrast, figsize=(20, 8))
    for i, p in enumerate(ax.patches):
        aa = e_df1.iloc[i]['aa']
        h = p.get_height()
        if h > 0:
            ax.text(p.get_x(), p.get_height() + 0.001, aa, rotation='vertical', va='bottom', fontsize=10, color='black')
        else:
            ax.text(p.get_x(), p.get_height() - 0.001, aa, rotation='vertical', va='top', fontsize=10, color='black')
    if e_barplot is None:
        e_barplot = e_df1[['codon', 'aa', contrast]]
    else:
        e_barplot = pd.merge(e_barplot, e_df1[['codon', 'aa', contrast]])
 
    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon and 3) {}".format(norm_type.upper()))
    ax.set_ylabel('({} - {}) / {} *'.format(sample, control, control))
    fig = ax.get_figure()
    e_codon_file = '{}/{}_per_codon.pdf'.format(esite_dir, contrast)
    print('Wiriting file: {}'.format(e_codon_file))
    fig.savefig(e_codon_file)
    plt.cla()
    plt.close(fig)
 
    e_df1 = e_df1.groupby('aa').sum().reset_index()
    ax = e_df1.plot.bar(x='aa', y=contrast)
    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon and 3) {}".format(norm_type.upper()), fontsize=6)
    ax.set_ylabel('({} - {}) / {} *'.format(sample, control, control))
 
    fig = ax.get_figure()
    e_aa_file = '{}/{}_per_aa.pdf'.format(esite_dir, contrast)
    print('Writing file: {}'.format(e_aa_file))
    fig.savefig(e_aa_file)
    plt.cla()
    plt.close(fig) 

print('Writing p-site results: {}'.format(p_outfile))
p_barplot.to_csv(p_outfile, sep="\t", index=False)
print('Writing a-site results: {}'.format(a_outfile))
a_barplot.to_csv(a_outfile, sep="\t", index=False)
print('Writing e-site results: {}'.format(e_outfile))
e_barplot.to_csv(e_outfile, sep="\t", index=False)

