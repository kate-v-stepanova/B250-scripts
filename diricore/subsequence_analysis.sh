#!/bin/bash
set -e;
set -u;

dataset_id=$1
species=$2;
projectname=$dataset_id

if [[ $# -ge 3 ]]; then
  minreads=$3;
else
  minreads=15
fi

bam_type="hq_unique"
# can be: hq, hq_unique, all, all_unique
if [[ $# -ge 4 ]]; then
   bam_type=$4
fi
echo "$bam_type"

subset='all_samples'
if [[ $# -ge 5 ]]; then
    subset=$5
fi

plots_only=0
if [[ $# -ge 6 ]]; then
    plots=$6
    if [[ $plots -eq "plots_only" ]]; then
        plots_only=1
    fi
fi
echo "$plots_only"

y_limits=1
if [[ $# -eq 7 ]]; then
    y_limits=$7
fi

PROJECT_DIR="$BASE_DIR/$dataset_id"

OUTDIR="$PROJECT_DIR/analysis/output/subsequence_data";
INDIR="$PROJECT_DIR/analysis/output/alignments/toGenome"
PLOTDIR="$PROJECT_DIR/analysis/output/figures/subsequence_shift_plots/$bam_type/";
if [[ $subset == 'all_samples' ]]; then
    CONTRASTS="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts.tsv";
    SAMPLENAMES="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv";
else
    CONTRASTS="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts_$subset.tsv";
    SAMPLENAMES="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames_$subset.tsv";
    PLOTDIR="$PLOTDIR/$subset"
fi

DIRICORE_DIR="$BASE_DIR/software/diricore"
INDEXDATAFN="$BASE_DIR/static/${species}/subseq_index_data.pkl.gz";

all_frame_file="${OUTDIR}/${projectname}.subsequence_data.all.${minreads}.hdf5"
all_dedup_frame_file="${OUTDIR}/${projectname}.subsequence_data.all.dedup.${minreads}.hdf5"
hq_frame_file="${OUTDIR}/${projectname}.subsequence_data.hq.${minreads}.hdf5"
hq_dedup_frame_file="${OUTDIR}/${projectname}.subsequence_data.hq.dedup.${minreads}.hdf5"

if [ $bam_type == "all" ]; then
    frame_file=$all_frame_file
    bam_pattern="_toGenome.bam"
elif [ $bam_type == "all_unique" ]; then
    frame_file=$all_dedup_frame_file
    bam_pattern="_toGenome_dedup.bam"
elif [ $bam_type == "hq" ]; then
    frame_file=$hq_frame_file
    bam_pattern="_toGenome.hqmapped.bam"
else
    frame_file=$hq_dedup_frame_file
    bam_pattern="_toGenome.hqmapped_dedup.bam"
fi

###
mkdir -p ${OUTDIR}
mkdir -p "${PLOTDIR}"

if [[ -f $frame_file ]]; then
    if [ "$plots_only" -eq 0 ]; then
        echo "File exists: $frame_file"
        echo "Delete file or use option 'plots_only'"
        exit
    fi
fi

if [ $plots_only -eq 0 ]; then
    # rm -f $OUTDIR/*
    rm -f $frame_file;
    echo "Extracting subsequences ($bam_type)";
    for bamfn in $(ls $INDIR/*$bam_pattern); do
        b=$(basename $bamfn);
        b=${b%"$bam_pattern"};
        $DIRICORE_DIR/bin/extract_subsequences.py \
           -v \
           run \
           -o $frame_file \
           -f 0 \
           ${INDEXDATAFN} \
           ${b},${bamfn}
    done;
    echo "Done. Generated file: ${frame_file}";
fi


echo "Creating plots"
if [[ -f $frame_file ]]; then
    # create subsequence shift plots
    mkdir -p ${PLOTDIR}
    $DIRICORE_DIR/bin/plot_subsequence_shifts.py \
    -o ${PLOTDIR}/${projectname}.$bam_type.m${minreads}. \
    -m $minreads \
    --sample-names ${SAMPLENAMES} \
    --contrasts ${CONTRASTS} \
    --y-limits ${y_limits} \
    ${frame_file}
    echo "Created plots in ${PLOTDIR}"
else
    echo "ERROR: file $frame_file not found!! Skipping"
fi
###

