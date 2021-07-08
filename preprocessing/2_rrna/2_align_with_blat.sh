project_id=$1
genome=$2

indir="$BASE_DIR/$project_id/analysis/output/fastq_to_fasta/"
ref_file="$BASE_DIR/static/$genome/rRNA_genes.2bit"
outdir="$BASE_DIR/$project_id/analysis/output/rrna/blat_results/"
blat_bin="$BASE_DIR/software/bin/blat_for_linux"

for f in $(ls $indir/*.fa); do 
    fn=$(basename $f); 
    fn=${fn%.fa}; 
    echo "bsub -q long -R \"rusage[mem=50G]\" $blat_bin -stepSize=5 -repMatch=4096 -minScore=0 -minIdentity=0 $ref_file $f $outdir/${fn}.psl"; 

done
