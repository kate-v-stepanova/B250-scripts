#!/bin/bash

set -e

project_id=$1
BASE_PATH="/icgc/dkfzlsdf/analysis/OE0532"
PROJECT_PATH="$BASE_PATH/$project_id"
DIRICORE_PATH="$BASE_PATH/software/diricore"
blat_path="$DIRICORE_PATH/programs/blat_for_linux"
PSL_PARSER="$BASE_PATH/software/ext_diricore/parse_psl.py"

genome="hg19"
if [[ $# -ge 2 ]]; then
   genome=$2
fi

bam_type="all" # all or all_unique. Do not use hq
if [[ $# -ge 3 ]]; then
    bam_type=$3
fi

GENOME_PATH="$BASE_PATH/static/$genome/longest_trans/longest_transcripts.2bit"
INDIR="$PROJECT_PATH/analysis/output/ext_diricore/$bam_type/fasta"
OUTDIR="$PROJECT_PATH/analysis/output/ext_diricore/$bam_type/psl"
script_dir="$BASE_PATH/tmp/$project_id/blat_alignment"
mkdir -p $OUTDIR
mkdir -p $script_dir

for fasta_file in $(ls $INDIR/*.fasta); do
    sample_name=$(basename $fasta_file)
    sample_name=${sample_name%.fasta}
    psl_file="$OUTDIR/${sample_name}_aligned.with_header.psl"
    psl_no_header="$OUTDIR/${sample_name}.no_header.psl"
    script_file="$script_dir/${sample_name}.sh"
    
    echo "#!/bin/bash" > $script_file
    echo "set -e" >> $script_file
    echo "echo Alignment of $fasta_file" >> $script_file
    echo "if [[ ! -f $psl_file ]]; then $blat_path -stepSize=5 -repMatch=4096 -minScore=0 -minIdentity=0 $GENOME_PATH $fasta_file $psl_file; else echo File exists. Skipping $psl_file; fi" >> $script_file
    #$blat_path -minIdentity=80 $GENOME_PATH $fasta_file $psl_file

    # remove first 5 lines (header)
    echo "tail -n +6 $psl_file > $psl_no_header" >> $script_file
    output_file="$OUTDIR/${sample_name}.parsed_psl.txt"
    echo "cat $psl_no_header | cut -f14,16,17,15,9 > $output_file" >> $script_file
    #echo "echo Alignment done. Processing $psl_no_header" >> $script_file
    #output_file="$OUTDIR/${sample_name}.parsed_psl.txt"
    #echo "python $PSL_PARSER $psl_no_header > $output_file" >> $script_file
    echo "echo Created $output_file" >> $script_file
    chmod +x $script_file
    echo "bsub -q long $script_file"
done;
