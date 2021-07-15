import pandas as pd
import glob
import os
import sys
import itertools

x = int(sys.argv[1]) # number of repetitions of codons. E.g. x=1 means we look for "ATG" and x=2 means we look for "ATGATG" and so on. 
# x can be theoretically any number, but for now it's intended to look either for duplets or triplets (x=2 and x=3)
genome = 'hg19'
if len(sys.argv) == 3:
    genome = sys.argv[2]

BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
transcripts = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/transcriptLength.txt".format(genome)
trans_df = pd.read_csv(transcripts, sep="\t")
cds_file = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/transcripts_cds.fa".format(genome)

cds = pd.read_csv(cds_file, sep=">", header=None)
seq = cds[0].dropna()
trans = cds[1].dropna()
cds = pd.DataFrame()
cds['transcript'] = trans.tolist()
cds['cds'] = seq.tolist()
trans_df = pd.merge(trans_df, cds, on="transcript", how="inner")

codons_file = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/codons.txt".format(genome)
codons = pd.read_csv(codons_file, sep="\t", header=None, names=["codon", "AA", "aa", "amino_acid"])

outdir = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/codons/{}x".format(genome, x)
os.makedirs(outdir, exist_ok=True)
def find_next(row):
    codon = row['codon']
    seq = row['cds'][row['pos']+len(codon):]
    return seq.find(codon) + row['pos'] + len(codon)

amino_acids = codons['AA'].unique()
for aa in amino_acids:
    aa_df = codons.loc[codons['AA'] == aa]
    codons1 = set(itertools.permutations(aa_df['codon'].tolist() * x, r=x))
    codons2 = [''.join(cod) for cod in codons1]
    full_df = pd.DataFrame()
    for codon in codons2:
        codon_df = cds.copy()
        codon_df['codon'] = codon
        codon_df['pos'] = codon_df['cds'].str.find(codon)
        # get all transcripts where position of the occurence is the position of actual codon 
        # pos % len(codon) == 0. e.g. pos = 3 and codon = "CCA" 
        # Means that CCA starts at position 3 and pos 3 is a start of a codon: (0 1 2) (3...
        # e.g. pos = 4. Means that CCA doesnt start at the first position of the codon (0 1 2) (3 4 ..) 
        # we keep only those rows where the occurence is at the codon first position
        # in-frame codons
        codon_df = codon_df.loc[codon_df['pos'] % len(codon) == 0] # first occurence
        df = codon_df.copy()
        ## df will be empty if no occurences are found
        while len(df) != 0:
            # find next occurence of the codon
            df['pos2'] = df.apply(find_next, axis=1)
            # keep only those rows where the codon found at the first position of a codon (in-frame?)
            df = df.loc[df['pos2'] % len(codon) == 0]
            df['pos'] = df['pos2']
            # and append to the codons_df
            codon_df = codon_df.append(df.drop('pos2', axis=1), ignore_index=True)
        full_df = full_df.append(codon_df.drop('cds', axis=1), ignore_index=True)
    outfile = os.path.join(outdir, "{}_codon_positions.tsv".format(aa.lower()))
    print("Writing file: {}".format(outfile))

    full_df.to_csv(outfile, sep="\t", header=False, index=False)
