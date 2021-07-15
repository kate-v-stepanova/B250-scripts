#!/usr/bin/env python
import pandas as pd
import glob
import os
import sys
import matplotlib.pyplot as plt
import itertools

project_id = sys.argv[1]
bam_type = sys.argv[2]
aa = sys.argv[3] # e.g. pro, asp, leu ...

window = 30
if len(sys.argv) >= 5:
    window = int(sys.argv[4])
x = 1
if len(sys.argv) >= 6:
    x = int(sys.argv[5])

genes = []
if len(sys.argv) >= 7:
    gene = sys.argv[6]
    if gene == "all_genes":
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

samples = []
if len(sys.argv) >= 9:
    if sys.argv[8] in ['all', 'all_samples']:
        samples = []
    else:
        samples = [sys.argv[8]]

max_y = None

codons_file = '/icgc/dkfzlsdf/analysis/OE0532/static/hg19/codons.txt'
codons_df = pd.read_csv(codons_file, sep="\t", header=None, usecols=[0,1])
codons_df = codons_df.loc[codons_df[1].str.lower() == aa]
codons = codons_df[0].tolist()
codons = codons + codons
codons = list(itertools.permutations(codons, x))
codons = [''.join(c) for c in codons]
codons = list(set(codons))

# select only coding genes
coding = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/coding_genes_names.txt"
coding_df = pd.read_csv(coding, skiprows=1, header=None, names=['gene_name'])
coding_df = coding_df.drop(coding_df.loc[coding_df['gene_name'].str.startswith('MT-')].index)

trans_file = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/transcriptLength.txt"
trans_df = pd.read_csv(trans_file, sep="\t")
if genes:
    trans_df = trans_df.loc[trans_df['gene_name'].isin(genes)]

trans_df = pd.merge(trans_df, coding_df, on='gene_name', how='inner')
trans_df = trans_df[['transcript', 'cds_len', '5utr_len']]
trans_df.columns = ['trans', 'cds_len', '5utr_len']

trans = trans_df['trans'].tolist()
indir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/tsv".format(project_id, bam_type)
outdir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/rpf_coverage_cumulative/{}".format(project_id, bam_type, subset_name)
plot_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}//analysis/output/figures/ext_diricore/rpf_coverage_cumulative/{}/{}".format(project_id, bam_type, subset_name)
os.makedirs(outdir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

p_offset_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/p_offset".format(project_id, bam_type)
if not os.path.isdir(p_offset_dir) or len(os.listdir(p_offset_dir)) == 0:
    print('P-sites offset info is required!!!. Dir is empty: {}'.format(p_offset_dir))
    print('Run the following command to calculate offsets: (plastid has to be installed)')
    print('\tfor f in $(ls /icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/alignments/toGenome/*_dedup.bam); do fn=$(basename $f); fn=${{fn%_toGenome_dedup.bam}};  echo "bsub -q medium  -R \"rusage[mem=10G]\" psite /icgc/dkfzlsdf/analysis/OE0532/static/{}/plastid_rois.txt  /icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/plastid/p_offsets/$fn --min_length 15 --max_length 32 --require_upstream --count_files $f --title \"$fn\""; done'.format(project_id, genome, project_id))
    sys.exit(1)
# We don't care about P-sites!!! Just take everything starting from DD

def process_sample(f):
    sample = os.path.basename(f).replace('.tsv', '')
    outfile = "{}/{}.tsv".format(outdir, sample)
    if os.path.exists(outfile):
        print('File exists and will not be overwritten: {}'.format(outfile))
        return
    df = pd.read_csv(f, sep="\t", header=None, names=['trans', 'start', 'seq'])
    df['start'] = df['start'].astype(int)
    # convert 1-based coordinates into 0-based coordinates
    df['start'] = df['start'] - 1
    df = df.loc[df['trans'].isin(trans)]
    # convert coordinates. 0 is the beginning of cds. 
    df = pd.merge(df, trans_df, on='trans', how='inner')
    df['start'] = df['start'] - df['5utr_len']
    # remove reads which are not in cds or last 30 nts of 5' UTR
    df = df.loc[df['start'] >= - window]
    # remove reads which are not in cds or first 30 nts of 3' UTR
    df = df.loc[df['start'] <= df['cds_len'] + window]
    # remove reads which are not in frame
    df = df.loc[df['start'] % 3 == 0]
    # get psite offset
    offset_file = "{}/{}_p_offsets.txt".format(p_offset_dir, sample)
    offset_df = pd.read_csv(offset_file, sep="\t")
    for i, row in offset_df.iterrows():
        l = int(row['length'])
        offset = int(row['p_offset'])
        df.loc[df['seq'].str.len() == l, 'p_codon'] = df.loc[df['seq'].str.len() == l, 'seq'].str[offset:offset+x*3] 
        df.loc[df['seq'].str.len() == l, 'psite_pos'] = df.loc[df['seq'].str.len() == l, 'start'] + offset

    # remove reads with length for which p-site offset was not defined
    df = df.dropna()
    df = df.loc[df['p_codon'].isin(codons)]
    
    # norm by cds length + 2 * window
    df['p_norm_pos'] = df['psite_pos'] / (df['cds_len'] + window * 2)
    import pdb; pdb.set_trace() 

    # group by start position (normalized)
 
    # group by start pos
    df['counts'] = 1
    df1 = df.groupby(['trans', 'start'])['counts'].sum().reset_index()
    import pdb; pdb.set_trace()

    # assing windows to each read. one row - one window. one read - multiple windows
    w_df = pd.merge(df, codon_df, on=['trans'], how='inner')
    # remove rows which don't fit into a window
    w_df = w_df.loc[w_df.apply(lambda row: row['start'] >= row['pos'] - window and row['start'] <= row['pos'] + window, axis=1)]
    # get position in a window
    w_df['window_pos'] = w_df['start'] - w_df['pos']
    w_df.columns = ['transcript', 'read_pos', 'counts', 'codon', 'codon_pos', 'window_pos']
    # save results
    print('Writing file: {}'.format(outfile))
    w_df.to_csv(outfile, sep="\t", header=True, index=False)

def plot_only(sample):
    pf = os.path.join(outdir, 'plot_{}.tsv'.format(sample))
    df = pd.read_csv(pf, sep='\t')
    plt_file = os.path.join(plot_dir, '{}.pdf'.format(sample))
    print('Writing file: {}'.format(plt_file))
    fig = df.plot(kind='line', x='pos', y='norm_counts')
    fig.set_title('Sample: {}. Amino Acid(s): {}. Gene(s): {}. # windows: {}'.format(sample, aa, subset_name, windows))
    plt.xlabel('Distance from codon')
    plt.ylabel("Average RPF 5' ends")
    if max_y is not None:
        plt.ylim(0, max_y)
    plt.savefig(plt_file)
        

def plot_sample(f):
    sample = os.path.basename(f).replace('.tsv', '')
    df = pd.read_csv(f, sep="\t")
    avg_df = None
    windows = 0
    for t in df['transcript'].unique().tolist():
        for pos in df.loc[df['transcript'] == t, 'codon_pos'].unique().tolist():
            window_df = df.loc[(df['transcript'] == t) & (df['codon_pos'] == pos)]
            # intrawindow normalization
            reads_per_window = window_df['counts'].sum()
            if reads_per_window <= 25:
                continue
            window_df['norm_counts'] = (window_df['counts'] / reads_per_window) * (window * 2 + 1)
            df1 = pd.DataFrame(columns=['pos', 'counts', 'norm_counts'])
            df1['pos'] = [i for i in range(window * -1, window + 1)]
            for i, row in window_df.iterrows():
                df1.loc[df1['pos'] == row['window_pos'], 'counts'] = row['counts']
                df1.loc[df1['pos'] == row['window_pos'], 'norm_counts'] = row['norm_counts']
            df1 = df1.fillna(0)
            if avg_df is None:
                avg_df = df1
            else:
                avg_df['counts'] = avg_df['counts'] + df1['counts']
                avg_df['norm_counts'] = avg_df['norm_counts'] + df1['norm_counts'] 
            windows += 1
    # get an average among all windows
    avg_df['counts'] = avg_df['counts'] / windows # divide by the number of windows
    avg_df['norm_counts'] = avg_df['norm_counts'] / windows # divide by the number of windows
    plot_values_file = os.path.join(outdir, 'plot_{}.tsv'.format(sample))
    print('Writing file: {}'.format(plot_values_file))
    avg_df.to_csv(plot_values_file, sep='\t', index=False)
    plt_file = os.path.join(plot_dir, '{}.pdf'.format(sample))
    print('Writing file: {}'.format(plt_file))
    fig = avg_df.plot(kind='line', x='pos', y='norm_counts')
    fig.set_title('Sample: {}. Amino Acid(s): {}. Gene(s): {}. # windows: {}'.format(sample, aa, subset_name, windows))
    plt.xlabel('Distance from codon')
    plt.ylabel("Average RPF 5' ends")
    if max_y is not None:
        plt.ylim(0, max_y)
    plt.savefig(plt_file)



if len(samples) == 1:
    f = os.path.join(indir, "{}.tsv".format(samples[0]))
    if not os.path.exists(f):
        print('File does not exist: {}'.format(f))
        exit(1)
    pf = os.path.join(outdir, 'plot_{}.tsv'.format(samples[0]))
    if os.path.exists(pf):
        print('Plot data already exists: {}'.format(pf))
        print('WARNING: will not overwrite the file')
        plot_only(samples[0])      
    else:
        process_sample(f)
        f2 = os.path.join(outdir, '{}.tsv'.format(samples[0]))
        plot_sample(f2)
else:
    for f in glob.glob(os.path.join(indir, "*.tsv")):
        sample = os.path.basename(f).replace('.tsv', '')
        if max_y is None:
            max_y = ''
        out_string = 'bsub -q long -R "rusage[mem=50G]" python {}'.format(' '.join(sys.argv[:8]) + ' {} '.format(sample) + str(max_y))
        print(out_string)
    exit(0)
