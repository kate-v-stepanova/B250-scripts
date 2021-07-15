import pandas as pd
import sys
import glob
import os

BASE_DIR = os.getenv('BASE_DIR')

project_id = sys.argv[1] # 17246
bam_type = sys.argv[2]
min_reads = sys.argv[3]

subset = ''
if len(sys.argv) >= 5:
    subset = sys.argv[4]
if subset:    
    indir = os.path.join(BASE_DIR, project_id, "analysis/output/diff_expr/", "{}_reads{}/{}".format(bam_type, min_reads, subset))
else:
    indir = os.path.join(BASE_DIR, project_id, "analysis/output/diff_expr/", "{}_reads{}".format(bam_type, min_reads))

if not os.path.isdir(indir):
    print('Dir does not exist: {}'.format(indir))
    exit()

if len(glob.glob(indir)) == 0:
    print('Dir is empty: {}'.format(indir))
    exit()

coding_df = None
contrasts = []
for f in glob.glob(os.path.join(indir, 'coding_diff_expr_*.tsv')):
    contrast = os.path.basename(f).replace('coding_diff_expr_', '').replace('.tsv', '')
    contrasts.append(contrast)
    df = pd.read_csv(f, sep="\t")
    df = df[['gene', 'log2FoldChange']]
    df.columns = ['gene', contrast]
    if coding_df is None:
        coding_df = df
    else:
        coding_df = pd.merge(coding_df, df, how='outer', on='gene')


full_df = None
for f in glob.glob(os.path.join(indir, 'diff_expr_*.tsv')):
    contrast = os.path.basename(f).replace('diff_expr_', '').replace('.tsv', '')
    df = pd.read_csv(f, sep="\t")
    df = df[['gene', 'log2FoldChange']]
    df.columns = ['gene', contrast]
    if full_df is None:
        full_df = df
    else:
        full_df = pd.merge(full_df, df, how='outer', on='gene')


full_df = full_df.fillna(0)
coding_df = coding_df.fillna(0)

full_df = full_df[['gene'] + contrasts]
coding_df = coding_df[['gene'] + contrasts]

# sorting by sum of fold changes per gene
idx1 = full_df[contrasts].sum(axis=1).abs().sort_values(ascending=False).index.tolist()
idx2 = coding_df[contrasts].sum(axis=1).abs().sort_values(ascending=False).index.tolist()
full_df = full_df.loc[idx1]
coding_df = coding_df.loc[idx2]

outfile1 = os.path.join(indir, 'FC_diff_expr.tsv')
outfile2 = os.path.join(indir, 'coding_FC_diff_expr.tsv')
outfile3 = os.path.join(indir, "top50_coding_FC_diff_expr.tsv")

print('Writing file: {}'.format(outfile1))
print('Writing file: {}'.format(outfile2))
print('Writing file: {}'.format(outfile3))

full_df.to_csv(outfile1, sep="\t", index=False)
coding_df.to_csv(outfile2, sep="\t", index=False)
coding_df.head(50).to_csv(outfile3, sep="\t", index=False)
