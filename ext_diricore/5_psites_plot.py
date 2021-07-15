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
psite_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}//analysis/output/figures/ext_diricore/{}/psite_codons/{}".format(project_id, bam_type, norm_type)
asite_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}//analysis/output/figures/ext_diricore/{}/asite_codons/{}".format(project_id, bam_type, norm_type)
esite_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}//analysis/output/figures/ext_diricore/{}/esite_codons/{}".format(project_id, bam_type, norm_type)
os.makedirs(psite_dir, exist_ok=True)
os.makedirs(asite_dir, exist_ok=True)
os.makedirs(esite_dir, exist_ok=True)

psites_file = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites/psites.tsv".format(project_id, bam_type)
asites_file = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites/asites.tsv".format(project_id, bam_type)
esites_file = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites/esites.tsv".format(project_id, bam_type)
p_df = pd.read_csv(psites_file, sep="\t")
a_df = pd.read_csv(asites_file, sep="\t")
e_df = pd.read_csv(esites_file, sep="\t")
# drop Stop codons
p_df = p_df.loc[p_df['aa'] != 'Stp']
a_df = a_df.loc[a_df['aa'] != 'Stp']
e_df = e_df.loc[e_df['aa'] != 'Stp']

cols = [c for c in p_df.columns if norm_type in c]
for c in cols:
    sample = c.replace('{}_'.format(norm_type), '')
    # psite plots
    p_df1 = p_df[['codon', 'aa', c]]
    p_df1.columns = ['codon', 'aa', sample]
    p_df1 = p_df1.sort_values(by=sample)
    ax = p_df1.plot.bar(x='codon', y=sample, figsize=(20, 8))
    for i, p in enumerate(ax.patches):
        aa = p_df1.iloc[i]['aa']
        ax.text(p.get_x(), p.get_height() + 0.001, aa, va='bottom', rotation='vertical', fontsize=10, color='black')

    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon")
    ax.set_ylabel(norm_type.upper())

    fig = ax.get_figure()
    p_codon_file = '{}/{}_per_codon.pdf'.format(psite_dir, sample)
    print('Writing file: {}'.format(p_codon_file))
    fig.savefig(p_codon_file)
    plt.cla()
    plt.close(fig)
    p_df1 = p_df1.groupby('aa').sum().reset_index()
    p_df1 = p_df1.sort_values(by=sample)
    ax = p_df1.plot.bar(x='aa', y=sample)
    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon", fontsize=6)
    ax.set_ylabel(norm_type.upper())

    fig = ax.get_figure()
    p_aa_file =  '{}/{}_per_aa.pdf'.format(psite_dir, sample)
    print('Writing file: {}'.format(p_aa_file))
    fig.savefig(p_aa_file)
    plt.cla()
    plt.close(fig)
    # asite plots
    a_df1 = a_df[['codon', 'aa', c]]
    a_df1.columns = ['codon', 'aa', sample]
    a_df1 = a_df1.sort_values(by=sample)
    ax = a_df1.plot.bar(x='codon', y=sample, figsize=(20, 8))
    for i, p in enumerate(ax.patches):
        aa = a_df1.iloc[i]['aa']
        ax.text(p.get_x(), p.get_height(), aa, rotation='vertical', va='bottom', fontsize=10, color='black')
    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon")
    ax.set_ylabel(norm_type.upper())

    fig = ax.get_figure()
    a_codon_file = '{}/{}_per_codon.pdf'.format(asite_dir, sample)
    print('Writing file: {}'.format(a_codon_file))
    fig.savefig(a_codon_file)
    plt.cla()
    plt.close(fig)
 
    a_df1 = a_df1.groupby(['aa']).sum().reset_index()
    a_df1 = a_df1.sort_values(by=sample)
    ax = a_df1.plot.bar(x='aa', y=sample)
    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon", fontsize=6)
    ax.set_ylabel(norm_type.upper())
    fig = ax.get_figure()
    a_aa_file = '{}/{}_per_aa.pdf'.format(asite_dir, sample)
    print('Writing file: {}'.format(a_aa_file))
    fig.savefig(a_aa_file)
    plt.cla()
    plt.close(fig)
 
    # esite plots
    e_df1 = a_df[['codon', 'aa', c]]
    e_df1.columns = ['codon', 'aa', sample]
    e_df1 = e_df1.sort_values(by=sample)
    ax = e_df1.plot.bar(x='codon', y=sample, figsize=(20, 8))
    for i, p in enumerate(ax.patches):
        aa = e_df1.iloc[i]['aa']
        ax.text(p.get_x(), p.get_height(), aa, rotation='vertical', va='bottom', fontsize=10, color='black')
    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon")
    ax.set_ylabel(norm_type.upper())
 
    fig = ax.get_figure()
    e_codon_file = '{}/{}_per_codon.pdf'.format(esite_dir, sample)
    print('Wiriting file: {}'.format(e_codon_file))
    fig.savefig(e_codon_file)
    plt.cla()
    plt.close(fig)
 
    e_df1 = e_df1.groupby('aa').sum().reset_index()
    e_df1 = e_df1.sort_values(by=sample)
    ax = e_df1.plot.bar(x='aa', y=sample)
    ax.set_xlabel("* normalized 1) to the number of reads per codon per transcript and 2) to the number of transcripts per codon", fontsize=6)
    ax.set_ylabel(norm_type.upper())
 
    fig = ax.get_figure()
    e_aa_file = '{}/{}_per_aa.pdf'.format(esite_dir, sample)
    print('Writing file: {}'.format(e_aa_file))
    fig.savefig(e_aa_file)
    plt.cla()
    plt.close(fig) 
