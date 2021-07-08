


for f in $(ls $BASE_DIR/20252/analysis/output/fastq_to_fasta/*.fa); do 
    fn=$(basename $f); 
    fn=${fn%.fa}; 
    echo "bsub -q long -R \"rusage[mem=50G]\" /icgc/dkfzlsdf/analysis/OE0532/software/diricore/programs/blat_for_linux -stepSize=5 -repMatch=4096 -minScore=0 -minIdentity=0 $BASE_DIR/static/mm9/rRNA_genes.2bit $f $BASE_DIR/20252/analysis/output/rrna/blat_results/${fn}.psl"; 

done
