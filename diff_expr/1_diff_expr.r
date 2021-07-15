.libPaths(c('/omics/groups/OE0532/internal/RStudio_3.4', .libPaths()))

# install DESeq2 - has to be exactly this version
# wget https://bioconductor.org/packages/3.6/bioc/src/contrib/Archive/DESeq2/DESeq2_1.18.0.tar.gz
# install.packages("/omics/groups/OE0532/internal/RStudio_3.4/DESeq2_1.18.0.tar.gz", repos = NULL, type="source")
library(DESeq2) 
library(stats)
library(data.table)
library(dplyr)

BASE_DIR = Sys.getenv('BASE_DIR')

# R 3.4.2
args = commandArgs(trailingOnly=TRUE)
project_id=args[1]
bam_type = args[2]
min_reads = 1
if (length(args) >= 3) {
  min_reads = args[3]
}

indir <- paste(BASE_DIR, project_id, "/analysis/output/alignments/reads_per_gene/tsv/", bam_type, sep="")
diff_exp_dir = paste(BASE_DIR, project_id, "/analysis/output/diff_expr/", bam_type, "_reads", min_reads, sep="")
gene_names = paste(BASE_DIR, "/static/hg19/gene_names.txt", sep="")


# Get contrasts
contrasts_file <- paste(BASE_DIR, project_id, "/analysis/input/metadata/rpf_density_contrasts.tsv", sep="")
if (length(args) >= 4) {
  subset = args[4]
  contrasts_file <- paste(BASE_DIR, project_id, "/analysis/input/metadata/rpf_density_contrasts_", subset, ".tsv", sep="")
}
contrasts_dt <- fread(contrasts_file, header=F, sep="\t")
colnames(contrasts_dt) <- c("sample", "control", "color")
contrasts_dt$color <- NULL


selected_genes = c()
if (length(args) >= 5) {
   genes_file = args[5]
   subdir = basename(genes_file)
   subdir = gsub('.tsv', '', subdir)
   subdir = gsub('.txt', '', subdir)
   selected_genes = fread(genes_file, header = F)$V1
   diff_exp_dir = paste(diff_exp_dir, subdir, sep="/")
}

dir.create(diff_exp_dir, recursive=T, showWarnings = F)

# Get gene names
gene_names_dt <- as.data.frame(fread(gene_names, sep="|", header=F))
colnames(gene_names_dt) <- c("transcript_id", "gene_id", "gene_name")
#gene_names_dt$transcript_id <- NULL # discard
gene_names_dt <- gene_names_dt[!duplicated(gene_names_dt$gene_id), ]
rownames(gene_names_dt) <- gene_names_dt$gene_id

# get coding genes
path_to_coding_genes = paste(BASE_DIR, "/static/hg19/coding_genes.txt", sep="")
coding_genes <- as.vector(scan(path_to_coding_genes, character(), quote = "", skip=1))
gene_names_dt <- gene_names_dt[gene_names_dt$gene_id %in% coding_genes,]
coding_genes <- gene_names_dt$gene_name

selected_samples = c()

counts_df <- NULL
samples = c()
for (counts_file in list.files(indir, pattern="*.tsv")) {
  samplename = basename(counts_file)
  samplename = gsub(".tsv", "", counts_file)
  if (!(samplename %in% samples)) {
    samples = c(samples, samplename)
  }
  counts <- fread(paste(indir, counts_file, sep="/"), header = T, sep = "\t", col.names=c('gene', samplename))
  
  if (is.null(counts_df)) {
    counts_df = counts
  } else {
    counts_df = merge(counts_df, counts, by="gene")
  }
}

if (length(selected_samples) == 0) {
  selected_samples = samples
}


if (length(selected_genes) != 0) {
  counts_df = counts_df[counts_df$gene %in% selected_genes,]
}

# row names from genes
counts_df <- as.data.frame(counts_df)
counts_df[is.na(counts_df)] <- 0
rownames(counts_df) <- counts_df$gene
counts_df$gene <- NULL


# # remove genes with 0 counts in all samples
counts_df = counts_df[ rowSums(counts_df) > 0,]
# counts_df = counts_df[rowSums(counts_df >= min_reads) == ncol(counts_df),]
# counts_df <- counts_df[rowSums(counts_df >= min_reads) >= 1, ]

colData = NULL
colData$v1 <- as.vector(selected_samples)
colData <- as.data.frame(colData)
colData$v1 <- as.factor(colData$v1)
rownames(colData) = as.vector(selected_samples)

designFormula <- colData
dds <- DESeqDataSetFromMatrix(countData = counts_df, colData = colData, design = as.formula(designFormula))
dds <- DESeq(dds) # for normalization


print("Differential expression")
all_de_results = NULL
all_de_results$a <- counts_df$samplename
list_of_contrasts = c()
for(i in 1:nrow(contrasts_dt)) {
  row <- contrasts_dt[i,]
  contrast_name <- paste(row$sample, "__vs__", row$control, sep="")
  if (row$sample %in% colnames(dds) && row$control %in% colnames(dds)) {
    contrast_name <- paste(row$sample, "__vs__", row$control, sep="")
    list_of_contrasts <- c(list_of_contrasts, contrast_name)
    de_outfile = paste(diff_exp_dir, "/diff_expr_", contrast_name, ".tsv", sep="")
    DEresults = results(dds, contrast = c("v1", row$sample, row$control))
    DEresults <- DEresults[order(DEresults$pvalue),]
    de_results <- as.data.frame(DEresults)
    de_results$gene <- rownames(de_results)
    
    de_counts <- counts_df[c(row$sample, row$control)]
    # the number of columns with row sum should be 2? 
    # just removing the rows if in all columns counts are less than min_reads
    de_counts <- de_counts[rowSums(de_counts >= min_reads) >= ncol(de_counts),] 
    de_results[de_results$gene %in% rownames(de_counts),]
    
    coding_de <- de_results[de_results$gene %in% coding_genes,]
  
    # Save all DE results
    print(paste("Writing file:", de_outfile))
    write.table(de_results, file=de_outfile, sep="\t", quote = F, row.names = F)
    #de_results$contrast <- contrast_name
    all_de_results[[contrast_name]] = de_results
    
    # Save coding DE results
    coding_outfile = paste(diff_exp_dir, "/coding_diff_expr_", contrast_name, ".tsv", sep="")
    print(paste("Writing file:", coding_outfile))
    write.table(coding_de, file = coding_outfile, sep="\t", quote = F, row.names = F)
  } else {
    print(paste('No data for contrast', contrast_name))
  }
}
all_de_results$a <- NULL

