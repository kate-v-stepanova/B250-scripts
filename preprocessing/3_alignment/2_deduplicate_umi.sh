#!/bin/bash

set -e;
set -u;

dataset_id=$1
project_id=$dataset_id
genome_dir="toGenome"
if [[ $# -ge 3 ]]; then
   if [[ $3 == "--transcriptome" ]]; then
       genome_dir="toTranscriptome"
   elif [[ $3 == "--toTranscriptome" ]]; then
       genome_dir="toTranscriptome"
   elif [[ $3 == "--totranscriptome" ]]; then
       genome_dir="toTranscriptome"
   elif [[ $3 == "--trans" ]]; then
       genome_dir="toTranscriptome"
   elif [[ $3 == "trans" ]]; then
       genome_dir="toTranscriptome"
   elif [[ $3 == "toTranscriptome" ]]; then
       genome_dir="toTranscriptome"
   fi
fi
echo $genome_dir    

bam_type="hq"
# can be: hq, all
if [[ $# -ge 2 ]]; then
   if [[ $2 == "all" ]]; then
       bam_type="all"
   fi
fi 
echo $bam_type

mem="10G"
if [[ $# -ge 4 ]]; then
    mem=$4
fi


queue="long"
if [[ $# -ge 5 ]]; then
    queue=$5
fi

PROJECT_DIR="$BASE_DIR/$dataset_id"
INDIR="$PROJECT_DIR/analysis/output/alignments/$genome_dir"
OUTDIR="$PROJECT_DIR/analysis/output/alignments/$genome_dir"
script_dir="$BASE_DIR/tmp/deduplicating/$project_id"

mkdir -p $script_dir
mkdir -p $OUTDIR/logs


if [[ $bam_type == "hq" ]]; then
    bam_pattern="_${genome_dir}.hqmapped.bam"
    bam_prefix="_${genome_dir}.hqmapped_dedup.bam"
else
    bam_pattern="_${genome_dir}.bam"
    bam_prefix="_${genome_dir}_dedup.bam"
fi

for f in $(ls ${INDIR}/*$bam_pattern); do
    samplename=$(basename $f);
    samplename=${samplename%"$bam_pattern"}
    script_file="$script_dir/${samplename}${bam_pattern}.sh"
    ind_file="${OUTDIR}/${samplename}${bam_pattern}.bai";
    if [[ ! -f $ind_file ]]; then
        echo "module load samtools" > $script_file
        echo "samtools index ${f}" >> $script_file
    fi
    outfile=${OUTDIR}/${samplename}${bam_prefix}
    samplename=$(basename $f)
    samplename=${samplename%".bam"}
    echo "if [[ -f $outfile ]]; then echo file exist: $f; else umi_tools dedup --log2stderr -L ${OUTDIR}/logs/${samplename}_umidedup_log.txt -I ${f} -S $outfile --output-stats=${OUTDIR}/logs/${samplename}_dedup.log; fi" >> $script_file
    chmod +x $script_file
done;

if ls $script_dir/*${bam_pattern}.sh 1> /dev/null 2>&1; then
    for f in $(ls $script_dir/*${bam_pattern}.sh); do
        echo "bsub -q $queue -R \"rusage[mem=$mem]\" $f"
    done
else
    echo "Files already exist"
    ls $OUTDIR/*${bam_pattern}
fi
