#!/bin/bash

project_id=$1
bam_type=$2

indir="$BASE_DIR/$project_id/analysis/output/alignments/toGenome"
positions="$BASE_DIR/$project_id/analysis/input/metadata/positions.bed"
outdir="$BASE_DIR/$project_id/analysis/output/alignments_subset/$bam_type"
mkdir -p $outdir

if [[ $bam_type -eq 'all_unqiue' ]]; then
  bam_pattern='*_dedup.bam'
elif [[ $bam_type -eq 'all' ]]; then
  bam_pattern='*me.bam'
fi


echo $bam_pattern
for f in $(ls $indir/$bam_pattern); do
  fn=$(basename $f)
  echo "samtools view -b -L $positions $f > ${outdir}/${fn}"
done
