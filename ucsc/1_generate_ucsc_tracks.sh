#!/bin/bash

set -e
set -u

dataset_id=$1
if [ "$#" -ge 2 ]; then
   genome=$2
else
   genome="hg19"
fi

module load bedtools

bam_type="hq_unique"
bam_pattern="_toGenome.hqmapped_dedup.bam"
# can be: hq, hq_unique, all, all_unique
if [ $# -ge 3 ]; then
   bam_type=$3
fi

if [[ $bam_type == "hq" ]]; then
   bam_pattern="_toGenome.hqmapped.bam"
elif [[ $bam_type == "all" ]]; then
   bam_pattern="_toGenome.bam"
elif [[ $bam_type == "all_unique" ]]; then
   bam_pattern="_toGenome_dedup.bam"
fi

PROJECT_DIR="$BASE_DIR/$dataset_id"
SCRIPT_DIR="$BASE_DIR/tmp/ucsc/$dataset_id"

chrom_sizes="$BASE_DIR/static/$genome/$genome.chrom.sizes"
#INDIR="$PROJECT_DIR/analysis/output/tophat_out"
INDIR="$PROJECT_DIR/analysis/output/alignments/toGenome"
OUTDIR="$PROJECT_DIR/analysis/output/ucsc_tracks/$bam_type"

BDG2BW="$BASE_DIR/software/ucsc/bdg2bw.sh"

mkdir -p $SCRIPT_DIR
mkdir -p $OUTDIR
for f in $(ls $INDIR/*"$bam_pattern"); do
  # Getting bdg files
  samplename=$(basename $f)
  samplename=${samplename%"$bam_pattern"}
  script_file="$SCRIPT_DIR/${samplename}.sh"
  echo "#!/bin/bash" > $script_file
  echo "set -e;" >> $script_file
  echo "module load bedtools" >> $script_file
  echo "f=${f}" >> $script_file
  echo "chrom_sizes=${chrom_sizes}" >> $script_file
  echo "samplename=${samplename}" >> $script_file
  echo 'echo "Calculating normalization factor for $f"' >> $script_file;
  echo "fact=$\(samtools view -c $f\) &&" >> $script_file;
  echo 'fact=$\(echo "scale=6; 1000000.0 / $fact" | bc\)' >> $script_file;
  plus_file="$OUTDIR/${samplename}_plus.bdg";
  minus_file="$OUTDIR/${samplename}_minus.bdg"
  echo "plus_file=${plus_file}" >> $script_file
  echo "minus_file=${minus_file}" >> $script_file
  echo 'samtools view -b $f | tee >\(bedtools genomecov -ibam /dev/stdin -g $chrom_sizes -scale $fact -bg -split -strand + | sort -k1,1 -k2,2n > $plus_file\) | bedtools genomecov -ibam /dev/stdin -g $chrom_sizes -scale $fact -bg -split -strand - | sort -k1,1 -k2,2n > $minus_file' >> $script_file
  # Converting to bw
  echo 'echo "Converting to bigwig "' >> $script_file
  echo "$BDG2BW $plus_file $chrom_sizes" >> $script_file
  echo "$BDG2BW $minus_file $chrom_sizes" >> $script_file
  echo 'echo "Done $samplename"' >> $script_file
  sed 's/\\//g' $script_file > ${script_file}_1
  mv ${script_file}_1 ${script_file}
  chmod +x $script_file
  echo "bsub -q long -R \"rusage[mem=30G]\"" $script_file
done
