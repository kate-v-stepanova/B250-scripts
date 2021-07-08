import pandas as pd
import os
import sys

project_id = sys.argv[1]
BASE_DIR = os.getenv('BASE_DIR')

bc_split_file = "{}/{}/analysis/output/bc_split_stats.txt".format(BASE_DIR, project_id)
rrna_file = "{}/{}/analysis/output/alignment_stats/rrna_stats.txt".format(BASE_DIR, project_id)
alignment_file = "{}/{}/analysis/output/alignment_stats/alignment_stats.txt".format(BASE_DIR, project_id)
dedup_file = "{}/{}//analysis/output/alignment_stats/dedup_stats.txt".format(BASE_DIR, project_id)

outfile = '{}/{}//analysis/output/alignment_stats/all_stats.txt'.format(BASE_DIR, project_id)

bc_df = pd.read_csv(bc_split_file, sep="\t", usecols=[0,1])
rrna_df = pd.read_csv(rrna_file, sep="\t")
align_df = pd.read_csv(alignment_file, sep="\t", header=None, names=['sample', 'align_reads'])
dedup_df = pd.read_csv(dedup_file, sep="\t", header=None, names=['sample', 'dedup_reads'])

bc_df.columns = ['sample', 'bc_reads']
bc_df = bc_df.loc[~bc_df['sample'].str.endswith('unmatched')] # discard unmatched
bc_df = bc_df.loc[bc_df['sample'] != 'total'] # discard total

df = pd.merge(bc_df, rrna_df, how='outer')
df = pd.merge(df, align_df, how='outer')
df = pd.merge(df, dedup_df, how='outer')

print('Writing file: {}'.format(outfile))
df.to_csv(outfile, sep="\t", index=False)
