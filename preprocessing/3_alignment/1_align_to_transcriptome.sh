#!/bin/bash

set -e
set -u

dataset_id=$1
genome=$2
memory="40G"
if [[ $# -ge 3 ]]; then
    memory=$3
fi

PROJECT_DIR="$BASE_DIR/$dataset_id"
INDIR="$PROJECT_DIR/analysis/output/clean"
STAR_GENOME_DIR="$BASE_DIR/static/$genome"
OUTDIR=$PROJECT_DIR/analysis/output/alignments # after preprocessing
script_dir="$BASE_DIR/tmp/faster_alignment/$dataset_id"

mkdir -p $script_dir
mkdir -p $OUTDIR/toGenome
mkdir -p $OUTDIR/toTranscriptome
mkdir -p $OUTDIR/reads_per_gene
mkdir -p $OUTDIR/logs

for f in $(ls ${INDIR}/*.fastq.gz); do
    b=$(basename ${f});
    b=${b%%.*};
    mkdir -p $script_dir
    outfile="$script_dir/$b.sh"
    echo "#!/bin/bash" > $outfile
    echo "set -e" >> $outfile # exit on error, don't continue
    # align
    echo "module load STAR" >> $outfile
    echo "STAR --genomeDir $STAR_GENOME_DIR --runThreadN 100 --readFilesIn ${f} --outFileNamePrefix $OUTDIR/${b}_ --outSAMtype BAM Unsorted --readFilesCommand zcat  --quantMode TranscriptomeSAM GeneCounts --outSAMmapqUnique 0" >> $outfile
	# sort transcriptome
    echo "mv $OUTDIR/${b}_Aligned.toTranscriptome.out.bam ${OUTDIR}/toTranscriptome/${b}_toTranscriptome.bam" >> $outfile
    # samtools sort: option -T is important!! this is where the tmp will be written. Otherwise there will be problems with:
	#  1) running 2 datasets with equal sample names, the bam files will be mixed up!
    #  2) hard drive quota will be exceeded.
    echo "module load samtools/1.6" >> $outfile
	sorted="$OUTDIR/toTranscriptome/$b.sorted.trans.bam"
    echo "samtools sort $OUTDIR/toTranscriptome/${b}_toTranscriptome.bam -o $sorted -T $OUTDIR/trans.$b" >> $outfile
	echo "mv $sorted $OUTDIR/toTranscriptome/${b}_toTranscriptome.bam" >> $outfile

	# sort genome
	echo "mv ${OUTDIR}/${b}_Aligned.out.bam ${OUTDIR}/toGenome/${b}_toGenome.bam" >> $outfile
	sorted="$OUTDIR/toGenome/$b.sorted.bam"
	echo "samtools sort $OUTDIR/toGenome/${b}_toGenome.bam -o $sorted -T $OUTDIR/gen.$b" >> $outfile
	echo "mv $sorted $OUTDIR/toGenome/${b}_toGenome.bam" >> $outfile

	# mv log files
	echo "mv $OUTDIR/${b}*ReadsPerGene.out.tab $OUTDIR/reads_per_gene/" >> $outfile
	echo "mv $OUTDIR/${b}*.out $OUTDIR/logs/" >> $outfile
	echo "mv $OUTDIR/${b}*.tab $OUTDIR/logs/" >> $outfile

	chmod +x $outfile
	echo "bsub -q long  -R \"rusage[mem=$memory]\" $outfile"
done

