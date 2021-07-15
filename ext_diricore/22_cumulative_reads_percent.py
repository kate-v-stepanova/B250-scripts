import pandas as pd
import sys
import os
import glob
import matplotlib.pyplot as plt

project_id = sys.argv[1]
bam_type = sys.argv[2]
gene = sys.argv[3] # can be a gene name or a file with the list of genes

base_dir = "/icgc/dkfzlsdf/analysis/OE0532/"
project_dir = os.path.join(base_dir, project_id)
indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/tsv".format(bam_type))
#indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/psl".format(bam_type))
outdir = os.path.join(project_dir, "analysis/output/figures/cumulative_reads")

if not os.path.isdir(outdir):
    os.makedirs(outdir)

if os.path.isfile(gene):
    print('Reading file: {}'.format(gene))
    df = pd.read_csv(gene, sep="\t", header=None)
    genes = df[0].tolist()
    subset_name = os.path.basename(gene).replace('.txt', '').replace('.tsv', '')

else:
    genes = [gene]
    if len(sys.argv) >= 5:
        subset_name = "{}_{}".format(sys.argv[4], gene)
    else:
        subset_name = gene


no_title = False
x_axis = 'perc'
if len(sys.argv) >= 5:
   no_title = sys.argv[4] == 'no_title'
   x_axis = 'trans_len' if sys.argv[4] in ['trans', 'len', 'trans_len', 'trans_length'] else 'perc'
   
palette = plt.get_cmap('Set3') # max 12 samples

trans_len = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/transcriptLength.txt"
trans_df = pd.read_csv(trans_len, sep="\t")
#subset_name = os.path.basename(gene).replace('.txt', '').replace('.tsv', '')
#subset = pd.read_csv(subset_name, sep="\t", names=['transcript'], usecols=[0])
#trans_df = pd.merge(trans_df, subset, on="transcript", how="inner")
trans_df = trans_df.loc[trans_df['gene_name'].isin(genes)]

plt.rcParams["figure.figsize"] = (10,4)
plt.style.use('seaborn-darkgrid')
plt.xlabel('% start position')
plt.ylabel('% counts')
if no_title:
    plt.title("Reads distribution in CDS")
else:
    plt.title("Reads distribution in CDS. {}".format(subset_name))

full_df = pd.DataFrame()
full_df['start_pos'] = [i for i in range(101)]
full_df.index = [i for i in range(101)]
samples = []
for f in sorted(glob.glob(os.path.join(indir, "*.tsv"))):
    sample = os.path.basename(f).replace('.tsv', '')
    df = pd.read_csv(f, sep="\t", header=None, names=['transcript', 'start', 'seq'])
    df = pd.merge(trans_df, df, on='transcript', how='inner')
    df = df.loc[df['start'] > df['5utr_len']]
    df = df.loc[df['start'] < df['5utr_len'] + df['cds_len']]
    if len(df) == 0:
        print("No reads in CDS in {}".format(gene, f))
        continue
    df['start'] = df['start'] - df['5utr_len']
    df['start_pos'] = df['start'] / df['cds_len']
    if x_axis == 'perc':
        df['start_pos'] = df['start_pos'].round(2) * 100
    df['start_pos'] = df['start_pos'].astype(int)
    df1 = df[['start_pos']]
    df1['counts'] = 1

    df1['counts'] = df1.groupby(['start_pos'])['counts'].transform('sum')
    df1 = df1.drop_duplicates('start_pos')
    df1 = df1.sort_values(by="start_pos")
    df1.index = df1['start_pos']

    total = 0
    df1.index = df1['start_pos']
    for index, row in df1.iterrows():
        total += row['counts']
        df1.loc[index, 'counts'] = total

    full_df[sample] = df1['counts']
    samples.append(sample)

full_df = full_df.dropna()
full_df = full_df.fillna(0)
color = 0
for sample in samples:
    max_reads = full_df[sample].max()
    max_perc = full_df[sample].astype(float) / max_reads
    max_perc = max_perc * 100
    full_df[sample] = max_perc
    max_perc = max_perc.max()
    full_df.loc[full_df['start_pos'] == 1, sample] = full_df[(full_df['start_pos'] == 1) | (full_df['start_pos'] == 0)][sample].sum()
    # adding 0%
    if len(full_df.loc[full_df['start_pos'] == 0]) == 0:
        full_df.loc[-1, :] = 0
        # full_df.loc[-1, 'start_pos'] = 0
        full_df.index = full_df.index + 1
        full_df.sort_index(inplace=True)
    # adding 100% - maybe doesn't look great
    #if len(full_df.loc[full_df['start_pos'] == 100]) == 0:
    #    full_df.loc[100, ['start_pos', sample]] = 100
    #else:
    #    full_df.loc[full_df['start_pos'] == 100, sample] = 100

    plt.plot(full_df['start_pos'], full_df[sample], color=palette(color), label="{} ({} reads)".format(sample, str(int(max_reads))))
    color += 1
  
outfile = os.path.join(outdir, "cumulative_reads_{}.pdf".format(subset_name))
print("Saving plot: {}".format(outfile))
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.7, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(outfile)
 
outfile = outfile.replace('.pdf', '.tsv')
print('Saving plot: {}'.format(outfile))
full_df.to_csv(outfile, sep="\t", index=False)
