import pandas as pd
import numpy as np
import glob
import sys
import os

project_id = sys.argv[1]
bam_type = sys.argv[2]

indir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/tsv".format(project_id, bam_type)
outdir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/fasta".format(project_id, bam_type)
os.makedirs(outdir, exist_ok=True)

for f in glob.glob("{}/*.tsv".format(indir)):
    df = pd.read_csv(f, header=None, sep="\t")
    df = df[[0,2]]
    df[0] = ">" + df[0]
    outfile = f.replace('.tsv', '.fasta').replace(indir, outdir)
    print("Writing: {}".format(outfile))
    df.to_csv(outfile, sep="\n", index=False, header=False)

