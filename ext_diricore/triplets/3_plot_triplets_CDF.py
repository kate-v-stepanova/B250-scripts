import pandas as pd
import glob
import os
import sys
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np

project_id = sys.argv[1]
bam_type = sys.argv[2]
x = int(sys.argv[3])
site = int(sys.argv[4])

aa = sys.argv[5]

window = 0
if len(sys.argv) >= 7:
    window = int(sys.argv[6])

subset_name = ''
if len(sys.argv) >= 8:
    subset_name = '_' + sys.argv[7]

contrasts = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/input/metadata/rpf_density_contrasts{}.tsv".format(project_id, subset_name)
contrasts_df = pd.read_csv(contrasts, sep="\t", header=None, names=['sample', 'ctrl'], usecols=[0,1])

BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = os.path.join(BASE_DIR, project_id)
indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/codons_per_pos/{}x.window{}".format(bam_type, x, window))
outdir = os.path.join(indir, subset_name, "CDF")
plot_dir = os.path.join(project_dir, "analysis/output/figures/ext_diricore/{}/{}x.window{}/{}".format(bam_type, x, window, subset_name))

# create output dirs
os.makedirs(plot_dir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)

# get 2 groups of transcripts - with >= 1 DD sites per transcript and with 1 site per transcript
codon_file = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/codons/{}x/{}_codon_positions.tsv".format(x, aa)
codon_df = pd.read_csv(codon_file, sep="\t", header=None, names=['transcript', 'codon', 'pos'])
codon_df = codon_df[['transcript', 'codon']].groupby('transcript').count().reset_index()
group1 = codon_df.loc[codon_df['codon'] > 1]['transcript'].tolist()
group2 = codon_df.loc[codon_df['codon'] == 1]['transcript'].tolist()

# norm counts
norm_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/alignments/reads_per_gene/normalized/{}".format(project_id, bam_type)
if not os.path.isdir(norm_dir) or os.listdir(norm_dir) == []:
    print('Normalized counts not found: {}'.format(norm_dir))
    print('Get norm counts: python /icgc/dkfzlsdf/analysis/OE0532/software/scripts/normalize_counts.py {} {}'.format(project_id, bam_type))
    exit(1)

leg_labels = []
fig, ax = plt.subplots()
for i, row in contrasts_df.iterrows():
    sample = row['sample']
    ctrl = row['ctrl']
    contrast = '{}__vs__{}'.format(sample, ctrl)

    sample_file = '{}/{}_{}.tsv'.format(indir, site, sample)
    sample_df = pd.read_csv(sample_file, sep="\t")

    ctrl_file = '{}/{}_{}.tsv'.format(indir, site, ctrl)
    ctrl_df = pd.read_csv(ctrl_file, sep="\t")

    sample_df = sample_df.loc[sample_df['aa'] == aa]
    ctrl_df = ctrl_df.loc[ctrl_df['aa'] == aa]

    sample_df = sample_df[['transcript', 'codon_pos', 'counts']]
    ctrl_df = ctrl_df[['transcript', 'codon_pos', 'counts']]

    norm_sample = os.path.join(norm_dir, '{}.tsv'.format(sample))
    norm_ctrl = os.path.join(norm_dir, '{}.tsv'.format(ctrl))

    norm_s_df = pd.read_csv(norm_sample, sep="\t")
    norm_c_df = pd.read_csv(norm_ctrl, sep="\t")

    sample_df = pd.merge(sample_df, norm_s_df[['transcript', 'rpkm']], on='transcript', how='inner')
    ctrl_df = pd.merge(ctrl_df, norm_c_df[['transcript', 'rpkm']], on='transcript', how='inner')

    sample_df['rpkm'] = sample_df['counts'] * sample_df['rpkm']
    ctrl_df['rpkm'] = ctrl_df['counts'] * ctrl_df['rpkm']

    sample_df.columns = ['transcript', 'codon_pos', '{}_raw'.format(sample), '{}_norm'.format(sample)]
    ctrl_df.columns = ['transcript', 'codon_pos', '{}_raw'.format(ctrl), '{}_norm'.format(ctrl)]

    df = pd.merge(sample_df, ctrl_df, on=['transcript', 'codon_pos'])

    # threshold: >= 3 reads in at least one condition
    df = df.loc[(df['{}_raw'.format(sample)] >= 3) &(df['{}_raw'.format(ctrl)] >=3)]

    # calculate ratio
    df[contrast] = df[sample + '_norm'] / df[ctrl + '_norm']
    df[contrast] = np.log2(df[contrast])
    df = df.sort_values(by=contrast)

    df1 = df.loc[df['transcript'].isin(group1)] # >1 DD site per transcript
    df2 = df.loc[df['transcript'].isin(group2)] # =1 DD site per transcript

    df1['cdf'] = norm.cdf(df1[contrast])
    df2['cdf'] = norm.cdf(df2[contrast])
    df1_file = os.path.join(outdir, ">1DDsite_{}.tsv".format(contrast))
    df2_file = os.path.join(outdir, "1DDsite_{}.tsv".format(contrast))
    print('Writing file: {}'.format(df1_file))
    print('Writing file: {}'.format(df2_file))
    df1.to_csv(df1_file, sep="\t", index=False)
    df2.to_csv(df2_file, sep="\t", index=False)
    import pdb; pdb.set_trace()
    df1.plot(ax=ax, kind='line', x=contrast, y='cdf', label="{} (>1 DD site)".format(contrast))
    df2.plot(ax=ax, kind='line', x=contrast, y='cdf', label="{} (=1 DD site)".format(contrast))

plot_file = os.path.join(plot_dir, "{}{}.{}x.site{}.{}.pdf".format(aa, subset_name, x, site, bam_type))
print('Writing file: {}'.format(plot_file))
ax.legend()
fig = ax.get_figure()
fig.savefig(plot_file)
