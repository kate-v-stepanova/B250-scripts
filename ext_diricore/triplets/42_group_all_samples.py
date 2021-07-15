import pandas as pd
import glob
import sys
import os

project_id = sys.argv[1]
bam_type = sys.argv[2]
x = sys.argv[3]
site = sys.argv[4]
window = sys.argv[5]

indir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/codons/{}x.window{}".format(project_id, bam_type, x, window)

full_df = None
for f in glob.glob(os.path.join(indir, 'gene_name.{}*.tsv'.format(site))):
    df = pd.read_csv(f, sep="\t")
    import pdb; pdb.set_trace()
    sample = os.path.basename(f).replace('{}_'.format(site), '').replace('gene_name.', '').replace('.tsv', '')
    df.columns = ['aa',  'gene', sample,  'norm_{}'.format(sample)]
    df = df[['aa', 'gene', 'norm_{}'.format(sample)]]
    if full_df is None:
        full_df = df
    else:
        import pdb; pdb.set_trace()
        full_df = pd.merge(full_df, df, on=['aa', 'gene'], how='outer')
        full_df = full_df.fillna(0)

outfile = os.path.join(indir, 'all_samples.{}.tsv'.format(site))
print("Writing file {}".format(outfile))
full_df.to_csv(outfile, sep="\t", index=False)
