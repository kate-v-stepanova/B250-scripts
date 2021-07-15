#!/bin/bash
project_id=$1
bam_type=$2 # recommended: all_unique

indir="$BASE_DIR/$project_id/analysis/output/alignments/toTranscriptome"

if [[ $bam_type == "all_unique" ]]; then
    bam_pattern="_toTranscriptome_dedup.bam"
elif [[ $bam_type == "all" ]]; then
    bam_pattern="_toTranscriptome.bam"
elif [[ $bam_type == "hq_unique" ]]; then
    bam_pattern="_toTranscriptome.hqmapped_dedup.bam"
else
    echo "No settings for bam_type: $bam_type. Open the script and add it!!!"
    exit
fi

outdir="$BASE_DIR/$project_id/analysis/output/ext_diricore/$bam_type/tsv"
mkdir -p $outdir
 
for f in $(ls $indir/*$bam_pattern); do 
    fn=$(basename $f)
    fn=${fn%$bam_pattern}.tsv; 
    fn=$outdir/$fn
    if [[ ! -f $fn ]]; then
        samtools view $f | cut -f3,4,10 > $fn
        echo "File written: $fn"
    else
        echo "File exists. Skipping $fn"
    fi
done;
