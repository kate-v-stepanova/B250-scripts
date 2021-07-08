#!/bin/bash

BIN_PATH=$BASE_DIR/software/bin

# bedClip installed via conda install to diricore environment
# bedtools are in diricore/programs
# bedGraphToBigWig also in diricore/programs


# check commands: slopBed, bedGraphToBigWig and bedClip

#which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
#which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
#which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }

# end of checking

if [ $# -lt 2 ];then
    echo "Need 2 parameters! <bedgraph> <chrom info>"
    exit
fi

F=$1
G=$2

$BIN_PATH/bedtools slop -i ${F} -g ${G} -b 0 | $BIN_PATH/bedClip stdin ${G} ${F}.clip

LC_COLLATE=C sort -k1,1 -k2,2n ${F}.clip > ${F}.sort.clip

$BIN_PATH/bedGraphToBigWig ${F}.sort.clip ${G} ${F/bdg/bw}

#rm -f ${F}.clip ${F}.sort.clip
