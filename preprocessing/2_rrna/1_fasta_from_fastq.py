import pandas as pd
import os
import glob
import sys

project_id = sys.argv[1]

sample = None
if len(sys.argv) >= 3:
    sample = sys.argv[2]

BASE_DIR = os.getenv('BASE_DIR')
indir = '{}/{}/analysis/input/fastq'.format(BASE_DIR, project_id)
outdir = "{}/{}/analysis/output/fastq_to_fasta".format(BASE_DIR, project_id)

os.makedirs(outdir, exist_ok=True)

def fastq_to_fasta(f):
   sample = os.path.basename(f).replace('.fastq.gz', '')
   df = pd.read_csv(f, sep="\n", header=None)
   idx = range(1, len(df), 4)
   df = df.iloc[idx]
   df.columns = ['seq']
   df['counts'] = 1
   df = df.groupby('seq').agg({'counts': 'sum'}).reset_index()
   df['header'] = ">" + df['seq'] + '_' + df['counts'].astype(str)
   outfile = "{}/{}.fa".format(outdir, sample)
   print('Writing file: {}'.format(outfile))
   df[['header', 'seq']].to_csv(outfile, sep="\n", index=False, header=False)

if sample is None:
    for f in sorted(glob.glob('{}/*.fastq.gz'.format(indir))):
        sample = os.path.basename(f).replace('.fastq.gz', '')
        if 'unmatched' not in sample:
            print('bsub -q long  -R "rusage[mem=30G]" python {} {}'.format(' '.join(sys.argv), sample))
#        else:
#            print('Skipping unmatched: {}'.format(f))
else:
    fastq_to_fasta('{}/{}.fastq.gz'.format(indir, sample))
    
