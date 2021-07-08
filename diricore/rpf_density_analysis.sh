#!/bin/bash

set -e;
set -u;

dataset_id=$1

PROJECT_DIR="$BASE_DIR/$dataset_id"
OUTDIR="$PROJECT_DIR/analysis/output/rpf_5p_density";
INDIR="$PROJECT_DIR/analysis/output/alignments/toGenome"


species=$2;
projectname=$dataset_id;
project_id=$dataset_id

if [[ $# -ge 3 ]]; then
  minreads=$3
else
  minreads=100
fi


DIRICORE_DIR="$BASE_DIR/software/diricore"
INDEXDATA_FILE="$BASE_DIR/static/${species}/transcript_data.hdf5";
MAPS_FILE="$BASE_DIR/static/${species}/codon_regions.width_61.hdf5";
MAPSSTART_FILE="$BASE_DIR/static/${species}/codon_regions.START_Other_ATG.width_61.hdf5";

# Checking that files are present
ls $INDEXDATA_FILE
ls $MAPS_FILE
ls $MAPSSTART_FILE

###
#  setup
of_hq_unique="${OUTDIR}/${projectname}.txcoord_counts.hq.dedup.${minreads}.hdf5";
of_hq="${OUTDIR}/${projectname}.txcoord_counts.hq.${minreads}.hdf5";
of_all_unique="${OUTDIR}/${projectname}.txcoord_counts.all.dedup.${minreads}.hdf5";
of_all="${OUTDIR}/${projectname}.txcoord_counts.all.${minreads}.hdf5";

bam_type="hq_unique"
# can be: hq, hq_unique, all, all_unique
if [[ $# -ge 4 ]]; then 
    bam_type=$4;
    echo "bam_type: $bam_type"
fi
echo "$bam_type"

subset="all_samples"
if [[ $# -ge 5 ]]; then
    subset=$5
fi

plots_only=0
if [[ $# -ge 6 ]]; then
  plots_only=$6
  if [[ "$plots_only" == "plots_only" ]]; then
      plots_only=1
  fi
fi

PLOTDIR="$PROJECT_DIR/analysis/output/figures//rpf_5p_density_plots//$bam_type"
if [[ $subset == "all_samples" ]]; then
    CONTRAST_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts.tsv";
    SAMPLENAME_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames.tsv";
else
    CONTRAST_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_contrasts_${subset}.tsv";
    SAMPLENAME_FILE="$PROJECT_DIR/analysis/input/metadata/rpf_density_samplenames_${subset}.tsv";
    PLOTDIR="$PLOTDIR/$subset/"
fi

echo $PLOTDIR
mkdir -p ${OUTDIR}
mkdir -p ${PLOTDIR}

if [ "$bam_type" == "hq" ]; then
    outfile="$of_hq"
    bam_pattern="_toGenome.hqmapped.bam"
elif [ "$bam_type" == "all_unique" ]; then
    outfile=$of_all_unique
    bam_pattern="_toGenome_dedup.bam"
elif [ $bam_type == "all" ]; then
    outfile=$of_all
    bam_pattern="_toGenome.bam"
else
    outfile=$of_hq_unique
    bam_pattern="_toGenome.hqmapped_dedup.bam"
fi

echo $plots_only

if [[ -f $outfile ]]; then
    if [ "$plots_only" -eq 0 ]; then
        echo "File exists: $outfile"
        echo "Delete the file or use option 'plots_only'"
        exit
    fi
fi

if [ "$plots_only" -eq 0 ]; then
 # map RPFs to transcriptome coordinates
 echo "Mapping RPFs to transcripttome coordinates ($bam_type)"
 bam_files=$(ls ${INDIR}/*$bam_pattern);
 echo $bam_files;
 $DIRICORE_DIR/bin/map_rpfs_to_transcriptome_positions.py -t "${INDEXDATA_FILE}" -o ${outfile} -b <(ls -1 ${INDIR}/*$bam_pattern | sort -V | while read fn; do b=$(basename $fn); b=${b%"$bam_pattern"}; echo -e "${b}\t${fn}"; done)
 echo "Mapping done. Created file: $outfile"
fi

output="\
Ala\tGCA,GCC,GCG,GCT
Arg\tCGA,CGC,CGG,CGT,AGA,AGG
Asn\tAAC,AAT
Asp\tGAC,GAT
Cys\tTGC,TGT
Gln\tCAA,CAG
Glu\tGAA,GAG
Gly\tGGA,GGC,GGG,GGT
His\tCAC,CAT
Ile\tATA,ATC,ATT
Leu\tCTA,CTC,CTG,CTT,TTA,TTG
Lys\tAAA,AAG
Met\tATG
Phe\tTTC,TTT
Pro\tCCA,CCC,CCG,CCT
Ser\tTCA,TCC,TCG,TCT,AGC,AGT
Thr\tACA,ACC,ACG,ACT
Trp\tTGG
Tyr\tTAC,TAT
Val\tGTA,GTC,GTG,GTT
"

echo "Generating RPF density shift plots ($bam_type)"
echo -ne "$output" | while read aa codongroupstr; do
    if [[ $subset == "all_samples" ]]; then
        plot_file="${PLOTDIR}/${projectname}.${bam_type}.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    else
         plot_file="${PLOTDIR}/${projectname}.${subset}.${bam_type}.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    fi

    if [[ ! -f $plot_file ]]; then
      codons=$(echo $codongroupstr | sed 's/,/ /g');

      python $DIRICORE_DIR/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${plot_file}" \
        -m ${minreads} \
        "${outfile}" \
        ${MAPS_FILE} \
        ${codons}
      echo "Shift plots done ($bam_type). Created file: $plot_file"
    else
      echo "File exists, skipping: $plot_file"
    fi
done

# RPF density at special codons
# (START/other ATG)
echo "RPF density at special codons ($bam_type)"
echo -ne "\
ATG_split\tSTART_ATG,Other_ATG
" | while read aa codongroupstr; do
    if [[ $subset == "all_samples" ]]; then
        plot_file="${PLOTDIR}/${projectname}.$bam_type.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    else
         plot_file="${PLOTDIR}/${projectname}.${subset}.$bam_type.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    fi

    if [[ ! -f $plot_file ]]; then
      codons=$(echo $codongroupstr | sed 's/,/ /g');
      python $DIRICORE_DIR/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${plot_file}" \
        -m ${minreads} \
        "${outfile}" \
        ${MAPSSTART_FILE} \
        ${codons}
      echo "Special codons done ($bam_type). Created file: $plot_file"
    else
      echo "File exists, skipping: $plot_file"
    fi
done
