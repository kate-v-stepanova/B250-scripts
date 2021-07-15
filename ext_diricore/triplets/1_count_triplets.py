import pandas as pd
import glob
import os
import sys
import matplotlib.pyplot as plt

project_id = sys.argv[1]
genome = sys.argv[2]
minreads = int(sys.argv[3])
bam_type = sys.argv[4]

x = 1 # number of repetitions of codons. E.g. x=1 means we look for "ATG" and x=2 means we look for "ATGATG" and so on. 
# x can be theoretically any number, but for now it's intended to look either for duplets or triplets (x=2 and x=3)
if len(sys.argv) >= 6:
    x = int(sys.argv[5])

site = 12 # p_site can be: 15: a_site, 18: e_site
if len(sys.argv) >= 7:
    site = int(sys.argv[6])

window = 0
if len(sys.argv) >= 8:
    window = int(sys.argv[7])

amino_acids = ['ala', 'arg', 'asn', 'asp', 'cys', 'gln', 'glu', 'gly', 'his', 'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'stp', 'thr', 'trp', 'tyr', 'val']
if len(sys.argv) >= 9:
    amino_acids = sys.argv[8].lower().split(',')

#filename_preffix = ""
#if len(sys.argv) >= 10:
#    gene = sys.argv[9]
#    if os.path.isfile(gene):
#        print('Reading file: {}'.format(gene))
#        df = pd.read_csv(gene, sep="\t", header=None)
#        genes = df[0].tolist()
#        subset_name = os.path.basename(gene).replace('.txt', '').replace('.tsv', '')
#    else:
#        genes = [gene]
#        subset_name = gene
#    filename_preffix = "{}_".format(subset_name)

#if len(genes) > 1:
#    out_preffix = "/{}".format(subset_name)
#else:
#    out_preffix = ""
#

BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = os.path.join(BASE_DIR, project_id)
indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/tsv".format(bam_type))
outdir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/codons/{}x.window{}".format(bam_type, x, window))
plot_dir = os.path.join(project_dir, "analysis/output/figures/ext_diricore/{}/{}x.window".format(bam_type, x, window))

#outdir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/codons/{}x.window{}/{}/".format(bam_type, x, window, out_preffix))
#plot_dir = os.path.join(project_dir, "analysis/output/figures/ext_diricore/{}/{}x.window{}/".format(bam_type, x, window, out_preffix))
    
codons_dir = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/codons/{}x".format(genome, x)
transcripts = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/transcriptLength.txt".format(genome)
trans_df = pd.read_csv(transcripts, sep="\t")

gene_names = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/gene_names.txt".format(genome)
genes_df = pd.read_csv(gene_names, sep="|", header=None, names=['trans', 'gene_id', 'gene_name'])
#genes_df = genes_df.loc[genes_df['gene_name'].isin(genes)]
#selected_transcripts = genes_df['trans'].tolist()

# create output dirs
os.makedirs(outdir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

palette = plt.get_cmap('Set3') # max 12 samples
plt.rcParams["figure.figsize"] = (10,4)
plt.style.use('seaborn-darkgrid')
plt.xlabel('% start position')
plt.ylabel('% counts')

color = 0
for f in sorted(glob.glob(os.path.join(indir, '*.tsv'))):
    sample = os.path.basename(f).replace('.tsv', '')
    full_df = pd.DataFrame()
    df = pd.read_csv(f, sep="\t", header=None, names=['transcript', 'start', 'seq'])
    df['start'] = df['start'].astype(int) - 1 + site # samtools view produce 1-based coordinates (starting from 1)
    # remove transcripts with number of reads < minreads
    counts = df['transcript'].value_counts()
#    df = df.loc[df['transcript'].isin(counts[counts >= minreads].index)]
    total_reads = 0
    for aa in amino_acids:
        filename = os.path.join(codons_dir, '{}_codon_positions.tsv'.format(aa))
        if not os.path.isfile(filename):
            print("File not found! {}".format(filename))
            continue
        codons_df = pd.read_csv(filename, sep="\t", header=None, names=['transcript', 'codon', 'codon_pos'])
        # filter by selected transcripts
#        codons_df = codons_df.loc[codons_df['transcript'].isin(selected_transcripts)]
        # remove transcripts which are not in the codon_df 
        df1 = pd.merge(df, codons_df, on="transcript", how="inner")
        if len(df) == 0:
            continue
        df1['end'] = df1['start'] + df1['seq'].str.len()
        df1['codon_end'] = df1['codon_pos'] + df1['codon'].str.len()
        if window == 0:
            df1 = df1.loc[(df1['start'] <= df1['codon_pos']) & (df1['end'] >= df1['codon_end'])]
        else:
            df1 = df1.loc[(df1['start'] >= df1['codon_pos'] - window) & (df1['end'] <= df1['codon_pos'] + window)]
        df2 = df1.groupby(['transcript', 'codon'], as_index=False).size().reset_index()
        df2.columns = ['transcript', 'codon', 'counts']
        df2['aa'] = aa
        df2 = df2.loc[df2['counts'] >= minreads]
        full_df = full_df.append(df2, ignore_index=True)
        total_reads += df2['counts'].sum()
#    import pdb; pdb.set_trace()
    # full_df['counts'] = full_df['counts'] / total_reads # normalize to the total number of reads per sample
    # full_df['counts'] = full_df['counts'].round(2) * 100
    full_df['norm_counts'] = full_df['counts'] / len(df)
    if len(full_df) == 0:
        print('No reads for sample {}'.format(sample))
        continue
    outfile = os.path.join(outdir, "{}_{}.tsv".format(site, sample))
    print("Writing data file: {}".format(outfile))
    full_df.to_csv(outfile, sep="\t", index=False)
    full_df1 = full_df.groupby('codon', as_index=False).size().reset_index()
    full_df1.columns = ['codon', 'counts']
    plt.plot(full_df1['codon'], full_df1['counts'], color=palette(color), label="{}".format(sample))
    color += 1

ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.7, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plot_file = os.path.join(plot_dir, "global_codons.pdf")
plt.xticks(rotation=90)
plt.savefig(plot_file)
