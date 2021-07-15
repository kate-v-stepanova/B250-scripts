#!/usr/bin/env python
import pandas as pd
import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt

def get_reads_per_frame(df):
    # 5' utr reads: if psite pos < 5' utr len
    utr5_df = df.loc[df['psite_pos'] < 0]
    utr5_reads = len(utr5_df)
    # 3' utr + cds reads
    cds_utr3_df = df.loc[df['psite_pos'] >= 0]

    # cds reads: if psite pos < cds len
    cds_df = cds_utr3_df.loc[cds_utr3_df['psite_pos'] < cds_utr3_df['cds_len']]
    cds_reads = len(cds_df)
    # 3' utr reads: if psite pos >= cds len
    utr3_df = cds_utr3_df.loc[cds_utr3_df['psite_pos'] >= cds_utr3_df['cds_len']]
    utr3_reads = len(utr3_df)

    # count reads per frame for 5' UTR reads
    utr5_frame0 = len(utr5_df.loc[utr5_df['psite_pos'] % 3 == 0])
    utr5_frame1 = len(utr5_df.loc[utr5_df['psite_pos'] % 3 == 1])
    utr5_frame2 = len(utr5_df.loc[utr5_df['psite_pos'] % 3 == 2])

    # reads per frame for CDS
    cds_frame0 = len(cds_df.loc[cds_df['psite_pos'] % 3 == 0])
    cds_frame1 = len(cds_df.loc[cds_df['psite_pos'] % 3 == 1])
    cds_frame2 = len(cds_df.loc[cds_df['psite_pos'] % 3 == 2])

    # reads per frame for 3' UTR
    utr3_frame0 = len(utr3_df.loc[utr3_df['psite_pos'] % 3 == 0])
    utr3_frame1 = len(utr3_df.loc[utr3_df['psite_pos'] % 3 == 1])
    utr3_frame2 = len(utr3_df.loc[utr3_df['psite_pos'] % 3 == 2])

    return {
             'utr5_frame0': utr5_frame0,
             'utr5_frame1': utr5_frame1,
             'utr5_frame2': utr5_frame2,

             'cds_frame0': cds_frame0,
             'cds_frame1': cds_frame1,
             'cds_frame2': cds_frame2,
    
             'utr3_frame0': utr3_frame0,
             'utr3_frame1': utr3_frame1,
             'utr3_frame2': utr3_frame2
           }   

project_id = sys.argv[1]
bam_type = sys.argv[2]
min_reads = 100
if len(sys.argv) >= 4:
    min_reads = int(sys.argv[3])

genome = 'hg19'
if len(sys.argv) >= 5:
    genome = sys.argv[4]

p_offset_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/p_offset".format(project_id, bam_type)
if not os.path.isdir(p_offset_dir) or len(os.listdir(p_offset_dir)) == 0:
    print('P-sites offset info is required!!!. Dir is empty: {}'.format(p_offset_dir))
    print('Run the following command to calculate offsets: (plastid has to be installed)')
    print('\tfor f in $(ls /icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/alignments/toGenome/*_dedup.bam); do fn=$(basename $f); fn=${{fn%_toGenome_dedup.bam}};  echo "bsub -q medium  -R \"rusage[mem=10G]\" psite /icgc/dkfzlsdf/analysis/OE0532/static/{}/plastid_rois.txt  /icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/plastid/p_offsets/$fn --min_length 15 --max_length 32 --require_upstream --count_files $f --title \"$fn\""; done'.format(project_id, genome, project_id))
    sys.exit(1)

reads_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}//analysis/output/ext_diricore/{}/tsv".format(project_id, bam_type)
codons_file = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/codons.txt".format(genome)
codons_df = pd.read_csv(codons_file, sep="\t", header=None, names=['codon', 'Aa', 'aa', 'AA'])

trans_file = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/transcriptLength.txt".format(genome)
trans_file = "/icgc/dkfzlsdf/analysis/OE0532/static/{}/annotation_dt_ribowaltz.tsv".format(genome)

trans_df = pd.read_csv(trans_file, sep="\t")
output_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites".format(project_id, bam_type)
os.makedirs(output_dir, exist_ok=True)
counts_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/alignments/reads_per_gene/normalized/{}/".format(project_id, bam_type)
dotplot_dir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ext_diricore/{}/psites_dotplot".format(project_id, bam_type)
os.makedirs(dotplot_dir, exist_ok=True)

stats_df = pd.DataFrame(columns=['sample', 'total_reads', 'total_reads_after_threshold', 'reads_with_offset', 'reads_in_annotation', 'cds_reads', 'utr5_reads', 'utr3_reads', 'frame0_reads', 'frame1_reads', 'frame2_reads', 'psite_reads', 'asite_reads', 'esite_reads'])
frame_df = pd.DataFrame(columns=['sample', 'utr5_frame0', 'utr5_frame1', 'utr5_frame2', 'cds_frame0', 'cds_frame1', 'cds_frame2', 'utr3_frame0', 'utr3_frame1', 'utr3_frame2'])
psite_df = None
asite_df = None
esite_df = None
for f in sorted(glob.glob("{}/*.tsv".format(reads_dir))):
    print('Processing file: {}'.format(f))
    sample = os.path.basename(f).replace('.tsv', '')
    df = pd.read_csv(f, sep="\t", header=None, names=['transcript', 'pos', 'seq'])
    # get normalized counts
    counts_file = "{}/{}.tsv".format(counts_dir, sample)
    counts_df = pd.read_csv(counts_file, sep="\t")
    # aggregate by gene
    gene_df = counts_df.groupby('gene_name').agg({'counts': 'sum', 'cpm': 'sum', 'tpm': 'sum', 'rpkm': 'sum'}).reset_index()
    # discard genes with too few reads (by min_reads)
    gene_df = gene_df.loc[gene_df['counts'] >= min_reads]
    # remove from counts_df 
    genes = gene_df['gene_name'].unique().tolist()
    counts_df = counts_df.loc[counts_df['gene_name'].isin(genes)]
    # converting 1-based coordinates to 0-based coordinates
    df['pos'] = df['pos'] - 1
    total_reads = len(df)
    # discard reads with too few genes (based on counts_df)
    df = df.loc[df['transcript'].isin(counts_df['transcript'].tolist())]
    total_reads_after_threshold = len(df)
    # calculate length of the read
    df['length'] = df['seq'].str.len()
    # getting offset info
    offset_file = "{}/{}_p_offsets.txt".format(p_offset_dir, sample)
    offset_df = pd.read_csv(offset_file, sep="\t")
    lens = offset_df['length'].unique().tolist()
    # get rid of the reads for which offset is not defined
    df = df.loc[df['length'].isin(lens)]
    reads_with_offset = len(df)
    # annotation
    df = pd.merge(df, trans_df, how='left')
    reads_in_annotation = len(df.loc[~df['cds_len'].isna()])
    df = df.dropna()
    # transform coordinates to cds coordinates
    df['cds_pos'] = df['pos'] - df['5utr_len'] # if negative, then the read is from 5' UTR 
    # set psite, asite and esite pos for all reads
    for i, row in offset_df.iterrows(): # should be a few rows ~5-6
        length = row['length']
        p_offset = row['p_offset']
        a_offset = p_offset + 3
        e_offset = p_offset - 3
        # calculate positions of p_site, e_site and a_site
        df.loc[df['length'] == length, 'psite_pos'] = df.loc[df['length'] == length, 'cds_pos'].astype(int) + p_offset
        df.loc[df['length'] == length, 'asite_pos'] = df.loc[df['length'] == length, 'cds_pos'].astype(int) + a_offset
        df.loc[df['length'] == length, 'esite_pos'] = df.loc[df['length'] == length, 'cds_pos'].astype(int) + e_offset

        # calculate the codon of p_site, e_site and a_site
        df.loc[df['length'] == length, 'p_codon'] = df.loc[df['length'] == length, 'seq'].str[p_offset:p_offset+3]
        df.loc[df['length'] == length, 'a_codon'] = df.loc[df['length'] == length, 'seq'].str[a_offset:a_offset+3]
        df.loc[df['length'] == length, 'e_codon'] = df.loc[df['length'] == length, 'seq'].str[e_offset:e_offset+3]

    # get reads per frame
    reads_per_frame = get_reads_per_frame(df)
    reads_per_frame['sample'] = sample
    frame_df = frame_df.append([reads_per_frame])

    # utr5_reads
    utr5_reads = reads_per_frame['utr5_frame0'] + reads_per_frame['utr5_frame1'] + reads_per_frame['utr5_frame2']
    # cds_reads
    cds_reads = reads_per_frame['cds_frame0'] + reads_per_frame['cds_frame1'] + reads_per_frame['cds_frame2']
    # utr5_reads
    utr3_reads = reads_per_frame['utr3_frame0'] + reads_per_frame['utr3_frame1'] + reads_per_frame['utr3_frame2']

    # discard reads which map to 5'UTR
    df = df.loc[df['psite_pos'] >= 0]
    utr5_reads = reads_in_annotation - len(df)
    cds_with_utr3 = len(df)
    # discard reads which map to 3' UTR
    df = df.loc[df['psite_pos'] <= df['cds_len']] # considering end of psite
    cds_reads = len(df)
    utr3_reads = cds_with_utr3 - cds_reads
    cds_reads = len(df)
    # get number of frame-1 reads
    frame1_reads = len(df.loc[df['psite_pos'] % 3 == 1])
    # get number of frame-2 reads
    frame2_reads = len(df.loc[df['psite_pos'] % 3 == 2])
    # selecting only in-frame reads (frame-0)
    df = df.loc[df['psite_pos'] % 3 == 0]
    inframe_reads = len(df)
    # p-site codons
    p_df = df[['transcript', 'seq', 'length', 'cds_pos', 'psite_pos', 'p_codon']]
    p_df = p_df.loc[p_df['psite_pos'] >= 0] # discard reads where p-site not in cds
    psite_reads = len(p_df)
    p_df = pd.merge(p_df, codons_df, left_on='p_codon', right_on='codon', how='inner')
    ## Methionine codon and Start codon are the same - ATG. 
    ## To distinguish between them, we need to check the coordinates of p-site
    ## If the coordinate is 0, then it's a Start codon. Otherwise - Methionine
    p_df.loc[(p_df['codon'].str.startswith('ATG')) & (p_df['psite_pos'] == 0), ['Aa', 'aa', 'AA']] = ['Str', '*', 'Start']
    # merge with counts_df (to get Gene names)
    p_df = pd.merge(p_df, counts_df, on='transcript', how='inner')
    # for each codon - get number of unique GENES (not transcripts) per codon (considering Methionine and Start codon separetely). 
    # and divide norm_counts per codon by number of genes per codon.
    genes_per_codon = p_df.groupby(['codon', 'Aa'])['gene_name'].nunique().reset_index()
    genes_per_codon.columns = ['codon', 'Aa', 'genes_per_codon']
    # count p-sites for each codon
    p_df['p_counts'] = 1
    p_df = p_df.groupby(['gene_name', 'Aa', 'codon'])['p_counts'].count().reset_index()
    # merge with gene_df (to get gene counts)
    p_df = pd.merge(p_df, gene_df, on='gene_name', how='inner')
    # normalize to the total number of reads per codon per gene
    p_df['tpm_counts'] = p_df['p_counts'] / p_df['tpm']
    p_df['cpm_counts'] = p_df['p_counts'] / p_df['cpm']
    p_df['rpkm_counts'] = p_df['p_counts'] / p_df['rpkm']
    # normalize to the number of genes per codon (after we sum up, this gives us average)
    p_df = pd.merge(p_df, genes_per_codon, on=['codon', 'Aa'], how='inner')
    p_df['tpm_counts'] = p_df['tpm_counts'] / p_df['genes_per_codon']
    p_df['cpm_counts'] = p_df['cpm_counts'] / p_df['genes_per_codon']
    p_df['rpkm_counts'] = p_df['rpkm_counts'] / p_df['genes_per_codon']
    # aggregate psites for dotplot
    p_df_dotplot = p_df.groupby(['gene_name', 'Aa', 'codon']).agg({'p_counts': 'sum', 'tpm_counts': 'sum', 'cpm_counts': 'sum', 'rpkm_counts': 'sum'}).reset_index()
    outfile = "{}/psite_{}.tsv".format(dotplot_dir, sample)
    print('Writing data (per gene): {}'.format(outfile))
    p_df_dotplot.to_csv(outfile, sep="\t", index=False)
    # aggregate psites for heatmaps
    p_df = p_df.groupby(['Aa', 'codon']).agg({'p_counts': 'sum', 'tpm_counts': 'sum', 'cpm_counts': 'sum', 'rpkm_counts': 'sum'}).reset_index()
    # a-site codons
    a_df = df[['transcript', 'seq', 'length', 'cds_pos', 'asite_pos', 'a_codon']]
    a_df = a_df.loc[a_df['asite_pos'] >= 0] # discard reads where a-site not in cds
    asite_reads = len(a_df)
    a_df = pd.merge(a_df, codons_df, left_on='a_codon', right_on='codon', how='inner')
    ## Methionine codon and Start codon are the same - ATG. 
    ## To distinguish between them, we need to check the coordinates of p-site
    ## If the coordinate is 0, then it's a Start codon. Otherwise - Methionine
    a_df.loc[(a_df['codon'] == 'ATG') & (a_df['asite_pos'] == 0), ['Aa', 'aa', 'AA']] = ['Str', '*', 'Start']
    # merge with counts_df (to get Gene names)
    a_df = pd.merge(a_df, counts_df, on='transcript', how='inner')
    # for each codon - get number of unique GENES per codon. and divide norm_counts per codon by number of genes per codon.
    genes_per_codon = a_df.groupby(['codon', 'Aa'])['gene_name'].nunique().reset_index()
    genes_per_codon.columns = ['codon', 'Aa', 'genes_per_codon']
    # count a-sites for each codon
    a_df['a_counts'] = 1
    a_df = a_df.groupby(['gene_name', 'Aa', 'codon'])['a_counts'].count().reset_index()
    # merge with gene_df (to get gene counts)
    a_df = pd.merge(a_df, gene_df, on='gene_name', how='inner')
    # normalize to the total number of reads per codon per gene
    a_df['tpm_counts'] = a_df['a_counts'] / a_df['tpm']
    a_df['cpm_counts'] = a_df['a_counts'] / a_df['cpm']
    a_df['rpkm_counts'] = a_df['a_counts'] / a_df['rpkm']
    # normalize to the number of genes per codon (after we sum up, this gives us average)
    a_df = pd.merge(a_df, genes_per_codon, on=['codon', 'Aa'], how='inner')
    a_df['tpm_counts'] = a_df['tpm_counts'] / a_df['genes_per_codon']
    a_df['cpm_counts'] = a_df['cpm_counts'] / a_df['genes_per_codon']
    a_df['rpkm_counts'] = a_df['rpkm_counts'] / a_df['genes_per_codon']
    # dotplot
    a_df_dotplot = a_df.groupby(['gene_name', 'Aa', 'codon']).agg({'a_counts': 'sum', 'tpm_counts': 'sum', 'cpm_counts': 'sum', 'rpkm_counts': 'sum'}).reset_index()
    outfile = "{}/asite_{}.tsv".format(dotplot_dir, sample)
    print('Writing data (per gene, asite): {}'.format(outfile))
    a_df_dotplot.to_csv(outfile, sep="\t", index=False)
    # aggregated counts (for heatmap)
    a_df = a_df.groupby(['Aa', 'codon']).agg({'a_counts': 'sum', 'tpm_counts': 'sum', 'cpm_counts': 'sum', 'rpkm_counts': 'sum'}).reset_index()

    # e-site codons
    e_df = df[['transcript', 'seq', 'length', 'cds_pos', 'esite_pos', 'e_codon']]
    e_df = e_df.loc[e_df['esite_pos'] >= 0] # disccard reads where e-site not in cds
    esite_reads = len(e_df)
    e_df  = pd.merge(e_df, codons_df, left_on='e_codon', right_on='codon', how='inner')
    ## Methionine codon and Start codon are the same - ATG. 
    ## To distinguish between them, we need to check the coordinates of p-site
    ## If the coordinate is 0, then it's a Start codon. Otherwise - Methionine
    e_df.loc[(e_df['codon'] == 'ATG') & (e_df['esite_pos'] == 0), ['Aa', 'aa', 'AA']] = ['Str', '*', 'Start']
    # merge with counts_df (to get Gene names)
    e_df = pd.merge(e_df, counts_df, on='transcript', how='inner')
    # for each codon - get number of unique GENES per codon. and divide norm_counts per codon by number of genes per codon.
    genes_per_codon = e_df.groupby(['codon', 'Aa'])['gene_name'].nunique().reset_index()
    genes_per_codon.columns = ['codon', 'Aa', 'genes_per_codon']
    # get e-sites for each dodon
    e_df['e_counts'] = 1
    e_df = e_df.groupby(['gene_name', 'Aa', 'codon'])['e_counts'].count().reset_index()
    # merge with gene_df (to get gene counts)
    e_df = pd.merge(e_df, gene_df, on='gene_name', how='inner')
    # normalize (divide by the norm_counts per gene)
    e_df['tpm_counts'] = e_df['e_counts'] / e_df['tpm']
    e_df['cpm_counts'] = e_df['e_counts'] / e_df['cpm']
    e_df['rpkm_counts'] = e_df['e_counts'] / e_df['rpkm']
    # normalize to the number of reads per codon (after we sum up, this gives us average)
    e_df = pd.merge(e_df, genes_per_codon, on=['codon', 'Aa'], how='inner')
    e_df['tpm_counts'] = e_df['tpm_counts'] / e_df['genes_per_codon']
    e_df['cpm_counts'] = e_df['cpm_counts'] / e_df['genes_per_codon']
    e_df['rpkm_counts'] = e_df['rpkm_counts'] / e_df['genes_per_codon']
    # dotplot
    e_df_dotplot = e_df.groupby(['gene_name', 'Aa', 'codon']).agg({'e_counts': 'sum', 'tpm_counts': 'sum', 'cpm_counts': 'sum', 'rpkm_counts': 'sum'}).reset_index()
    outfile = "{}/esite_{}.tsv".format(dotplot_dir, sample)
    print('Writing data (per gene, esite): {}'.format(outfile))
    e_df_dotplot.to_csv(outfile, sep="\t", index=False)
    # aggregated counts (for heatmap)
    e_df = e_df.groupby(['Aa', 'codon']).agg({'e_counts': 'sum', 'tpm_counts': 'sum', 'cpm_counts': 'sum', 'rpkm_counts': 'sum'}).reset_index()

    # append stats
    stats_df = stats_df.append([{'sample': sample, 'total_reads': total_reads, 'total_reads_after_threshold': total_reads_after_threshold, 'reads_with_offset': reads_with_offset, 'reads_in_annotation': reads_in_annotation, 'utr5_reads': utr5_reads, 'utr3_reads': utr3_reads, 'cds_reads': cds_reads, 'frame0_reads': inframe_reads, 'frame1_reads': frame1_reads, 'frame2_reads': frame2_reads, 'psite_reads': psite_reads, 'asite_reads': asite_reads, 'esite_reads': esite_reads}])

    # rename cols
    cols = ['aa', 'codon', sample, 'tpm_{}'.format(sample), 'cpm_{}'.format(sample), 'rpkm_{}'.format(sample)]
    p_df.columns = cols
    a_df.columns = cols
    e_df.columns = cols

#    # normalize to the inframe reads
#    p_df['{}_norm'.format(sample)] = p_df[sample] / inframe_reads
#    a_df['{}_norm'.format(sample)] = a_df[sample] / inframe_reads
#    e_df['{}_norm'.format(sample)] = e_df[sample] / inframe_reads

    # append to the main df
    psite_df = p_df if psite_df is None else  pd.merge(psite_df, p_df, on=['aa', 'codon'], how='outer')
    asite_df = a_df if asite_df is None else  pd.merge(asite_df, a_df, on=['aa', 'codon'], how='outer')
    esite_df = e_df if esite_df is None else  pd.merge(esite_df, e_df, on=['aa', 'codon'], how='outer')

# sort by codon
psite_df = psite_df.sort_values(['aa', 'codon'])
esite_df = esite_df.sort_values(['aa', 'codon'])
asite_df = asite_df.sort_values(['aa', 'codon'])

stats_df = stats_df[['sample', 'total_reads', 'total_reads_after_threshold', 'reads_with_offset', 'reads_in_annotation', 'utr5_reads', 'utr3_reads', 'cds_reads', 'frame0_reads', 'frame1_reads', 'frame2_reads', 'psite_reads', 'asite_reads', 'esite_reads']]
frame_df = frame_df[['sample', 'utr5_frame0', 'utr5_frame1', 'utr5_frame2', 'cds_frame0', 'cds_frame1', 'cds_frame2', 'utr3_frame0', 'utr3_frame1', 'utr3_frame2']]

stats_file = "{}/stats.tsv".format(output_dir)
frame_file = "{}/reads_per_frame.tsv".format(output_dir)
psite_file = "{}/psites.tsv".format(output_dir)
asite_file = "{}/asites.tsv".format(output_dir)
esite_file = "{}/esites.tsv".format(output_dir)

print('Writing stats: {}'.format(stats_file))
stats_df.to_csv(stats_file, sep="\t", index=False)
print('Writing reads per frame: {}'.format(frame_file))
frame_df.to_csv(frame_file, sep="\t", index=False)
print('Writing psites: {}'.format(psite_file))
psite_df.to_csv(psite_file, sep="\t", index=False)
print('Writing asites: {}'.format(asite_file))
asite_df.to_csv(asite_file, sep="\t", index=False)
print('Writing esites: {}'.format(esite_file))
esite_df.to_csv(esite_file, sep="\t", index=False)    


