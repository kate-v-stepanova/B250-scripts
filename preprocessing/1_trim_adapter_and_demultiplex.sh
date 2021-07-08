#!/bin/bash

# exit on error:
set -e

project_id=$1
project_path="$BASE_DIR/$project_id"
# Check if merged file is present
merged_path="$project_path/analysis/input/merged"
echo "Checking if $merged_path is not empty:"
if [[ $(ls -A $merged_path) ]]; then
    echo "Fine"
else
    echo "Dir is empty. Exiting"
    exit
fi

subset=""
if [[ $# -ge 3 ]]; then
    subset=$3
fi

demultiplexed_path="$project_path/analysis/output/demultiplexed/$subset/"
umi_extract_path="$project_path/analysis/output/umi_extract/$subset/"
umi_extract_logs="${umi_extract_path}/logs"
fastq_path="$project_path/analysis/input/fastq"


bc_path="$project_path/analysis/input/metadata/"

cutadapt_trimming_stats="$project_path/analysis/output/cutadapt_trimming_stats.txt"
bc_split_stats="$project_path/analysis/output/bc_split_stats.txt"
stats_dir="$project_path/analysis/output"

adapter_sequence="AGATCGGAAGAGCACACGTCTGAA"
# adapter_sequence="TTCAGACGTGTGCTCTTCCGATCT" # reverse

bc_pattern="NNNNNNNNNN" # random sequence 5nt long + barcode
barcode_splitter="$BASE_DIR/software/bin/fastx_barcode_splitter.pl"

# Creating directory structure (if not exist)
mkdir -p $demultiplexed_path
mkdir -p $umi_extract_logs
mkdir -p $fastq_path

# remove adapter
min_len=30
if [[ $# -ge 2 ]]; then
   min_len=$2
fi

# should be just 1 file per subset. Or multiple files if no subset specified. 
# And if just 1 file and no subset, it should also work (*.fastq and bc_file_* should be named the same way)
for f in $(ls $merged_path/$subset.fastq.gz); do 
    echo "Processing $f"
    fn=$(basename $f);
    fn=${fn%.fastq.gz}
    trimmed_file="$merged_path/${fn}_trimmed.fastq.gz"
    if [[ -f $trimmed_file ]]; then
        echo "File exists! $trimmed_file"
    else
        cutadapt_trimming_stats="$stats_dir/${fn}_cutadapt_trimming_stats.txt"
        echo "Removing adapter: gzip -dc $f | cutadapt -u 3 -O 7 -m $min_len -j 10 -a $adapter_sequence --discard-untrimmed -o $trimmed_file - > $cutadapt_trimming_stats"
        gzip -dc $f | cutadapt -u 3 -O 5 -m $min_len -j 10 -a $adapter_sequence --discard-untrimmed -o $trimmed_file - > $cutadapt_trimming_stats
        echo "Removing adapter done. Outfiles: $trimmed_file. Stats file: $cutadapt_trimming_stats"
    fi
    bc_split_stats="$stats_dir/${fn}_bc_split_stats.txt"
    bc_file="$bc_path/bc_file_${fn}.txt"

    echo "Demultiplexing: cat $trimmed_file | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_ --suffix .fastq --eol > $bc_split_stats"
    zcat $trimmed_file | $barcode_splitter --bcfile $bc_file --prefix $demultiplexed_path/dem_  --suffix .fastq --eol > $bc_split_stats
    echo "Demultiplexing done"
    mv $demultiplexed_path/dem_unmatched.fastq $demultiplexed_path/dem_unmatched_${fn}.fastq
done

# remove umi (any random sequence of 5-nt length)
echo "Removing UMIs"
for f in $demultiplexed_path/*.fastq; do
    fb=$(basename ${f});
  b=${fb%%.*};
  if [[ "$b" != "*unmatched*" ]]; then
      echo "Processing ${f}";
      umi_tools extract --extract-method=string --3prime --bc-pattern=$bc_pattern --log=$umi_extract_logs/${b}.log --log2stderr --stdin=$demultiplexed_path/${fb} --stdout=$umi_extract_path/${b}_umi_extracted.fastq.gz
      echo "Creating a symlink to $fastq_path/${b#dem_}.fastq.gz"
      ln -s $umi_extract_path/${b}_umi_extracted.fastq.gz $fastq_path/${b#"dem_"}.fastq.gz
      echo "Done ${f}"
  fi
done;

