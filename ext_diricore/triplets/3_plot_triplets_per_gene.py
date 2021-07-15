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

amino_acids = ['ala', 'arg', 'asn', 'asp', 'cys', 'gln', 'glu', 'gly', 'his', 'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'stp', 'thr', 'trp', 'tyr', 'val']
if len(sys.argv) >= 7:
    amino_acids = sys.argv[6].lower().split(',')

genes = []
if len(sys.argv) >= 8:
    gene = sys.argv[7]
    if os.path.isfile(gene):
         print('Reading file: {}'.format(gene))
         subset_name = os.path.basename(gene).replace('.tsv', '').replace('.txt', '')
         genes = pd.read_csv(gene, sep="\t", header=None)
         genes = genes[0].tolist()
    else:
        subset_name = sys.argv[7]
#        genes = [gene]

genome = "hg19"
gene_names = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/gene_names.txt".format(genome)
genes_df = pd.read_csv(gene_names, sep="|", header=None, names=['trans', 'gene_id', 'gene_name'])
if len(genes) != 0:
	genes_df = genes_df.loc[genes_df['gene_name'].isin(genes)]
else:
    genes = genes_df['gene_name'].tolist()
selected_transcripts = genes_df['trans'].tolist()

BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = os.path.join(BASE_DIR, project_id)
indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/codons/{}x.window{}".format(bam_type, x, window))
outdir = os.path.join(indir, subset_name)
plot_dir = os.path.join(project_dir, "analysis/output/figures/ext_diricore/{}/{}x.window{}/{}".format(bam_type, x, window, subset_name))

# create output dirs
os.makedirs(plot_dir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)

palette = plt.get_cmap('Set3') # max 12 samples
full_df = None
for aa in amino_acids:
    aa_df = None
    input_files = sorted(glob.glob(os.path.join(indir, '{}*.tsv'.format(site))))
    if len(input_files) == 0:
        print('No input files found: {}'.format(indir))
    samples = []
    for gene in genes:
        gene_df = None
        full_gene_df = None
        gene_transcripts = genes_df.loc[genes_df['gene_name'] == gene]['trans'].tolist()
        gene_samples = []
        for f in input_files:
            sample = os.path.basename(f).replace('.tsv', '')
            samples.append(sample)
            df = pd.read_csv(f, sep="\t")
            df = df.loc[df['aa'] == aa]
            if len(df) == 0:
                print('No reads for {} in sample {}'.format(aa, sample))
                continue
            # filter by transcripts
            df = df.loc[df['transcript'].isin(gene_transcripts)]
            # remove zero counts
            df = df.loc[df['norm_counts2'] > 0]
            if len(df) == 0:
                print('No reads for aa {}, gene {}, sample {}'.format(aa, gene, sample))
                continue
            gene_samples.append(sample)
            gene_filename = os.path.join(outdir, "{}_{}_{}.tsv".format(aa, gene, sample))
            print('Writing file: {}'.format(gene_filename))
            df.to_csv(gene_filename, sep="\t", header=True, index=False)
            #df1 = df.groupby('codon', as_index=False)['norm_counts'].sum().reset_index()
            df1 = df[['codon', 'norm_counts2']]
            df1.columns = ['codon', sample]
            df1 = df1.groupby('codon', as_index=False).sum()
            if gene_df is None:
                gene_df = df1
            else:
                gene_df = pd.merge(gene_df, df1[['codon', sample]], on='codon', how='outer')
            df['sample'] = sample
            if full_gene_df is None:
                full_gene_df = df
            else:
#                full_gene_df = pd.merge(full_gene_df, df, on=['transcript', 'codon'], how='outer')
                full_gene_df = full_gene_df.append(df, ignore_index=True)
        if gene_df is None:
            print('No reads for gene {}, aa {}'.format(gene, aa))
            continue
        gene_df = gene_df.fillna(0)
        samples = list(set(samples))
        gene_plotname = os.path.join(plot_dir, '{}_{}_{}_{}.pdf'.format(subset_name, aa, site, gene))
        width = len(gene_df) * len(samples) * 0.02
        plt.rcParams['figure.figsize'] = (max(width, 10), 6)

        gene_df.plot(kind='bar', x='codon', y=gene_samples)
        plt.style.use('seaborn-darkgrid')
        plt.xlabel('codon')
        plt.ylabel('norm counts')
        plt.gcf().subplots_adjust(bottom=0.2)
        plt.title('Counts per codon.\nAmino Acid: {}. Site: {}. Gene: {}'.format(aa.upper(), site, gene)) 
        ax = plt.subplot(111)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.xticks(rotation=90)
        print('Writing plot: {}'.format(gene_plotname))
        plt.savefig(gene_plotname)   
        plt.close()

        if aa_df is None:
#        aa_df = df[['codon', 'norm_counts', 'sample']]
            aa_df = gene_df
        else:
            aa_df = pd.concat([aa_df, gene_df]).groupby(['codon']).sum().reset_index()
    if aa_df is None:
        print('No reads for {}'.format(aa.upper()))
        continue
    ## save data file: get by aa from the global file and group by codon
    samples = list(set(samples))
    width = len(aa_df) * len(samples) * 0.02
    plt.rcParams["figure.figsize"] = (max(width, 10),6)

    aa_df = aa_df.fillna(0)
    aa_df.plot(kind='bar', x='codon', y=samples)
    plt.style.use('seaborn-darkgrid')
    plt.xlabel('codon')
    plt.ylabel('norm counts')
    plt.gcf().subplots_adjust(bottom=0.2)

    plt.title('Norm counts per codon.\nAmino Acid: {}. Site: {}. Subset: {}'.format(aa, site, subset_name)) 
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plot_file = os.path.join(plot_dir, "{}.{}x.site{}.{}.pdf".format(aa, x, site, bam_type))
    print('Writing file: {}'.format(plot_file))
    plt.xticks(rotation=90)

    plt.savefig(plot_file)   
    plt.close()

    full_aa_df = aa_df.rename(columns={'codon': 'aa'})
    full_aa_df['aa'] = aa
    full_aa_df = full_aa_df.groupby(['aa'], as_index=False).agg('sum')
    full_aa_df['aa'] = aa

    if full_df is None:
        full_df = full_aa_df
    else:
        full_df = full_df.append(full_aa_df).reset_index()

if full_df is None:
    print('No reads for {}'.format(subset_name))

full_df = full_df.fillna(0)
plt.rcParams["figure.figsize"] = (max(width, 10),6)
full_df.plot(kind='bar', x='aa', y=samples)
plt.style.use('seaborn-darkgrid')
plt.xlabel('amino acid')
plt.ylabel('norm counts')
plt.gcf().subplots_adjust(bottom=0.2)

plt.title('Counts per amino acid: {}.\nSite: {}. Subset: {}'.format(aa, site, subset_name)) 
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plot_file = os.path.join(plot_dir, "per_aa.{}.{}x.site{}.{}.pdf".format(subset_name, x, site, bam_type))
print('Writing file: {}'.format(plot_file))
plt.xticks(rotation=90)

plt.savefig(plot_file)   
plt.close()

