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


BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = os.path.join(BASE_DIR, project_id)
indir = os.path.join(project_dir, "analysis/output/ext_diricore/{}/codons/{}x.window{}".format(bam_type, x, window))
plot_dir = os.path.join(project_dir, "analysis/output/figures/ext_diricore/{}/{}x.window{}".format(bam_type, x, window))

# create output dirs
os.makedirs(plot_dir, exist_ok=True)

palette = plt.get_cmap('Set3') # max 12 samples
for aa in amino_acids:
    full_df = None
    samples = []
    norm_df = None

    data_files = glob.glob(os.path.join(indir, '{}*.tsv'.format(site)))
    if len(data_files) == 0:
        print('No *.tsv files found in {}'.format(indir))
        exit()
   
    for f in sorted(data_files):
        sample = os.path.basename(f).replace('.tsv', '').replace('{}_'.format(site), '')
        samples.append(sample)
        df = pd.read_csv(f, sep="\t", header=None, names=['codon', 'counts', 'aa', 'norm_counts'])
        df = df.loc[df['aa'] == aa]
        df = df.loc[df['norm_counts'] > 0]
        
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
    
    width = len(full_df) * len(samples) * 0.02
    plt.rcParams["figure.figsize"] = (max(width, 10),6)
    full_df = full_df.fillna(0)
    full_df.plot(kind='bar', x='codon', y=samples)
    plt.style.use('seaborn-darkgrid')
    plt.xlabel('codon')
    plt.ylabel('norm counts')
    plt.gcf().subplots_adjust(bottom=0.2)

    plt.title('Counts per codon.\nAmino Acid: {}. Site: {}.'.format(aa, site)) 
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plot_file = os.path.join(plot_dir, "sum.{}.{}x.site{}.{}.pdf".format(aa, x, site, bam_type))
    print('Writing file: {}'.format(plot_file))
    plt.xticks(rotation=90)

    plt.savefig(plot_file)   
    plt.close()

    plt.rcParams["figure.figsize"] = (max(width, 10),6)

    dff = pd.DataFrame(columns=['sample', 'norm_counts'])
    for col in full_df.columns:
        if col != 'codon':
             dff = dff.append({'sample': col, 'norm_counts': full_df[col].astype(float).sum()}, ignore_index=True)
    import pdb; pdb.set_trace()
    full_df = norm_df.fillna(0)
    full_df.plot(kind='bar', x='codon', y=samples)
    plt.style.use('seaborn-darkgrid')
    plt.xlabel('codon')
    plt.ylabel('counts per codon')
    plt.gcf().subplots_adjust(bottom=0.2)

    plt.title('Counts per codon.\nAmino Acid: {}. Site: {}.'.format(aa, site)) 
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plot_file = os.path.join(plot_dir, "not_norm.{}.{}x.site{}.{}.pdf".format(aa, x, site, bam_type))
    print('Writing file: {}'.format(plot_file))
    plt.xticks(rotation=90)

    plt.savefig(plot_file)   
    plt.close()




#    import pdb; pdb.set_trace()
