import pandas as pd
import sys
import os
import glob
import matplotlib.pyplot as plt

# example: 18436
# quantify from a certain position in transcript (regardless CDS)

project_id = sys.argv[1]
bam_type = sys.argv[2]
gene_file = sys.argv[3]
subset_name = os.path.basename(gene_file).replace('.tsv', '').replace('.txt', '')
gene_df = pd.read_csv(gene_file, sep="\t", header=None, names=['gene', 'transcript', 'plot_pos'])
genes = gene_df['transcript'].tolist()

base_dir = "/icgc/dkfzlsdf/analysis/OE0532/"
project_dir = os.path.join(base_dir, project_id)
indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/tsv".format(bam_type))
#indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/psl".format(bam_type))
outdir = os.path.join(project_dir, "analysis/output/figures/cumulative_reads/{}".format(subset_name))

if not os.path.isdir(outdir):
    os.makedirs(outdir)

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

# we don't care if we pass a list of gene names, gene IDs or transcript IDs
trans_df1 = trans_df.loc[trans_df['gene_name'].isin(genes)]
if len(trans_df1) == 0:
    trans_df1 = trans_df.loc[trans_df['transcript'].isin(genes)]
    if len(trans_df1) == 0:
        trans_df1 = trans_df.loc[trans_df['gene_id'].isin(genes)]
trans_df = trans_df1


for trans in genes:
    start_pos = gene_df.loc[gene_df['transcript'] == trans, 'plot_pos'].max()
    max_pos = trans_df.loc[trans_df['transcript'] == trans, 'length'].max() - 1

    plt.rcParams["figure.figsize"] = (10,4)
    plt.style.use('seaborn-darkgrid')
    plt.xlabel('transcript length ({} - {})'.format(start_pos, max_pos))
    plt.ylabel('% counts')
    plt.title("Reads distribution. {}".format(trans))

    full_df = pd.DataFrame()
    full_df['start'] = range(start_pos, max_pos)
    samples = []
    total_reads = {}
    for f in sorted(glob.glob(os.path.join(indir, "*.tsv"))):
        sample = os.path.basename(f).replace('.tsv', '')
        df = pd.read_csv(f, sep="\t", header=None, names=['transcript', 'start', 'seq'])
        df = df.loc[df['transcript'] == trans]
        df = pd.merge(trans_df, df, on='transcript', how='inner')
        df = pd.merge(df, gene_df[['transcript', 'plot_pos']], on='transcript', how='inner')
  
        df = df.loc[df['start'] >= df['plot_pos']]
        df['start'] = df['start'].astype(int)
        sample_reads = len(df)
        total_reads[sample] = sample_reads
        df1 = df[['start']]
        df1['counts'] = 1
        df1['counts'] = df1.groupby(['start'])['counts'].transform('sum')

        df1 = df1.drop_duplicates()
        df1 = df1.sort_values(by="start")
        df1.index = df1['start']
        
        if len(df.loc[df['start'] == max_pos]) == 0:
            df1 = df1.append({'start': max_pos, 'counts': 0}, ignore_index=True)
        total = 0
        df1.index = df1['start']
        for index, row in df1.iterrows():
            total += row['counts']
            df1.loc[index, 'counts'] = total
        # full_df['s1'] = [0,5,25] = [0, 100%], full_df['s2'] = [0,10, 15] = [0,20,100%]
        full_df[sample] = df1['counts']
        samples.append(sample)

    full_df = full_df.reset_index()
    full_df = full_df.dropna()
    full_df = full_df.fillna(0)
    color = 0
    for sample in samples:
        max_reads = full_df[sample].max()
        max_reads2 = total_reads[sample]
        max_perc = full_df[sample].astype(float) / max_reads
        max_perc = max_perc * 100
        full_df[sample] = max_perc
        max_perc = max_perc.max()
        # full_df.loc[full_df['start'] == 1, sample] = full_df[(full_df['start'] == 1) | (full_df['start'] == 0)][sample].sum()
        # adding start_pos
        if len(full_df.loc[full_df['start'] == start_pos]) == 0:
            full_df.loc[-1, :] = 0
            full_df.loc[-1, 'start'] = start_pos
            full_df.index = full_df.index + 1
            full_df.sort_index(inplace=True)
        plt.plot(full_df['start'], full_df[sample], color=palette(color), label="{} ({} reads)".format(sample, str(int(max_reads))))
        color += 1
    outfile = os.path.join(outdir, "cumulative_reads_{}.pdf".format(trans))
    print("Saving plot: {}".format(outfile))
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.7, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(outfile)
 
    outfile = outfile.replace('.pdf', '.tsv')
    print('Saving plot: {}'.format(outfile))
    full_df.to_csv(outfile, sep="\t", index=False)
    plt.clf()
