#!/bin/bash

set -e
set -u

if [[ "$#" -eq 1 ]]; then
    dataset_id=$1
else
    echo "Usage example: $0 14548
	where 14548 is dataset_id" 
    exit 
fi;

username=$(whoami)

# Create directory structure
mkdir -p $BASE_DIR/$dataset_id/analysis/input/fastq
mkdir -p $BASE_DIR/$dataset_id/analysis/input/metadata
mkdir -p $BASE_DIR/$dataset_id/analysis/input/merged

mkdir -p $BASE_DIR/$dataset_id/analysis/output/demultiplexed
mkdir -p $BASE_DIR/$dataset_id/analysis/output/umi_extract/logs

touch $BASE_DIR/$dataset_id/analysis/input/metadata/bc_file.txt
touch $BASE_DIR/$dataset_id/analysis/input/metadata/rpf_density_contrasts.tsv
touch $BASE_DIR/$dataset_id/analysis/input/metadata/rpf_density_samplenames.tsv

# cp -R /midterm/$dataset_id/data/ $analysis_path/$dataset_id
echo "sftp -r ad+$username@ftp4midterm.dkfz.de:0$dataset_id/data/* $BASE_DIR/$dataset_id"
sftp -r ad+$username@ftp4midterm.dkfz.de:0$dataset_id/data/* $BASE_DIR/$dataset_id


