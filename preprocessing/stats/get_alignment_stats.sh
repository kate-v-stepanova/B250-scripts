#!/bin/bash

dataset_id=$1

PROJECT_DIR="$BASE_DIR/$dataset_id"
INDIR="$PROJECT_DIR/analysis/output/alignments/toGenome"
OUTDIR="$PROJECT_DIR/analysis/output/alignment_stats"
clean_dir="$PROJECT_DIR/analysis/output/clean"

mkdir -p ${OUTDIR}

# Get ALL stats
rm -f ${OUTDIR}/alignment_stats.txt
touch ${OUTDIR}/alignment_stats.txt
for inf in $(ls $INDIR/*_toGenome.bam); do
    echo -e "Doing file $inf"
    bn=$(basename $inf)
    bn=${bn%"_toGenome.bam"}
    c=$(samtools view -c $inf)
    echo -e "${bn}\t$c" >> ${OUTDIR}/alignment_stats.txt
done
echo "Created file: ${OUTDIR}/alignment_stats.txt"

# Get ALL unique stats
rm -f ${OUTDIR}/dedup_stats.txt
touch ${OUTDIR}/dedup_stats.txt
for inf in $(ls $INDIR/*_toGenome_dedup.bam); do
    bn=$(basename $inf)
    bn=${bn%"_toGenome_dedup.bam"}
    c=$(samtools view -c $inf)
    echo -e "${bn}\t$c" >> ${OUTDIR}/dedup_stats.txt
done
echo "Created file: ${OUTDIR}/dedup_stats.txt"
