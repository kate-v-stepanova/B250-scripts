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
         subset_name = sys.argv[8]
         genes = [gene]


gene_names = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/gene_names.txt".format(genome)
genes_df = pd.read_csv(gene_names, sep="|", header=None, names=['trans', 'gene_id', 'gene_name'])
if len(genes) != 0:
	genes_df = genes_df.loc[genes_df['gene_name'].isin(genes)]
selected_transcripts = genes_df['trans'].tolist()


BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = os.path.join(BASE_DIR, project_id)
indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/codons/{}x.window{}".format(bam_type, x, window))
outdir = os.path.join(indir, subset_name)
plot_dir = os.path.join(project_dir, "analysis/output/figures/ext_diricore/{}/{}x.window{}/{}".format(bam_type, x, window, subset_name))

# create output dirs
os.makedirs(plot_dir, exist_ok=True)

palette = plt.get_cmap('Set3') # max 12 samples
for aa in amino_acids:
    full_df = None
    samples = []
    norm_df = None
    input_files = sorted(glob.glob(os.path.join(indir, '{}*.tsv'.format(site))))
    if len(input_files) == 0:
        print('No input files found: {}'.format(indir))
    for f in input_files:
        sample = os.path.basename(f).replace('.tsv', '')
        samples.append(sample)
        df = pd.read_csv(f, sep="\t")
        df = df.loc[df['aa'] == aa]
        if len(df) == 0:
            print('No reads for {} in sample {}'.format(aa, sample))
            continue
        # filter by transcripts
        df = df.loc[df['transcript'].isin(selected_transcripts)]
        # remove zero counts
        df = df.loc[df['norm_counts'] > 0]
        # plot per gene
        if len(genes) > 1:
            for gene in genes:
                gene_trans = genes_df.loc[genes_df['gene_name'] == gene]
                df2 = df.loc[df['transcript'].isin(gene_trans)]
                if len(df2) == 0:
                    print('No reads for gene {}, sample {}'.format(gene, sample))
                    continue
                gene_filename = os.path.join(outdir, "{}_{}_{}.tsv".format(aa, gene, sample))
                df2.to_csv(gene_filename, sep="\t", header=True, index=False)
                

 
        df1 = df[['codon', 'norm_counts']]
        df1.columns = ['codon', sample]
        
        if full_df is None:
#        full_df = df[['codon', 'norm_counts', 'sample']]
            full_df = df1
        else:
            full_df = pd.merge(full_df, df1, on='codon', how='outer')
#            full_df = full_df.append(df[['codon', 'norm_counts', 'sample']], ignore_index=True)
        df2 = df[['codon', 'counts']]
        df2.columns = ['codon', sample]
        if norm_df is None:
            norm_df = df2
        else:
            norm_df = pd.merge(norm_df, df2, on='codon', how='outer')
    if full_df is None:
        continue
    width = len(full_df) * len(samples) * 0.02
    plt.rcParams["figure.figsize"] = (max(width, 10),6)

    full_df = full_df.fillna(0)
    full_df.plot(kind='bar', x='codon', y=samples)
    plt.style.use('seaborn-darkgrid')
    plt.xlabel('codon')
    plt.ylabel('norm counts')
    plt.gcf().subplots_adjust(bottom=0.2)

    plt.title('Counts per codon.\nAmino Acid: {}. Site: {}. Subset: {}'.format(aa, site, subset_name)) 
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plot_file = os.path.join(plot_dir, "{}{}.{}x.site{}.{}.pdf".format(filename_preffix, aa, x, site, bam_type))
    print('Writing file: {}'.format(plot_file))
    plt.xticks(rotation=90)

    plt.savefig(plot_file)   
    plt.close()

    plt.rcParams["figure.figsize"] = (max(width, 10),6)

    full_df = norm_df.fillna(0)
    full_df.plot(kind='bar', x='codon', y=samples)
    plt.style.use('seaborn-darkgrid')
    plt.xlabel('codon')
    plt.ylabel('counts per codon')
    plt.gcf().subplots_adjust(bottom=0.2)

    plt.title('Counts per codon.\nAmino Acid: {}. Site: {}. Subset: {}'.format(aa, site, subset_name)) 
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plot_file = os.path.join(plot_dir, "not_norm.{}{}.{}x.site{}.{}.pdf".format(filename_preffix, aa, x, site, bam_type))
    print('Writing file: {}'.format(plot_file))
    plt.xticks(rotation=90)

    plt.savefig(plot_file)   
    plt.close()




#    import pdb; pdb.set_trace()
