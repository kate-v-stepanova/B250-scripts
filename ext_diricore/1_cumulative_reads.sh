project_id=$1
bam_type=$2

BASEDIR="/icgc/dkfzlsdf/analysis/OE0532"
project_dir="$BASEDIR/$project_id"
indir="$project_dir/analysis/output/alignments/toTranscriptome"
outdir="$project_dir/analysis/output/ext_diricore/$bam_type/tsv"

mkdir -p $outdir

if [[ $bam_type == "hq" ]]; then
    bam_pattern="toTranscriptome.hqmapped.bam"
elif [[ $bam_type == "all_unique" ]]; then
    bam_pattern="toTranscriptome_dedup.bam"
elif [[ $bam_type == "all" ]]; then
    bam_pattern="toTranscriptome.bam"
else
    bam_type="hq_unique"
    bam_pattern="toTranscriptome.hqmapped_dedup.bam"
fi

module load samtools
for f in $(ls $indir/*$bam_pattern); do
    fn=$(basename $f);
    fn=${fn%_$bam_pattern}.tsv
    outfile=$outdir/$fn
#    echo "Writing $outfile"
    echo "samtools view $f | cut -f3,4,10 > $outfile"
done
