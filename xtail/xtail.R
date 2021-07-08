library(xtail) # devtools::install_github('FelixErnst/xtail-1')
library(data.table)

# R/4.0.0 

RP_id = '21221'
RNA_id = '21221_RNA'
bam_type = 'all_unique'

args = commandArgs(trailingOnly=TRUE)

RP_id = toString(args[1])
RNA_id = toString(args[2])
bam_type = args[3]

subset = ''
if (length(args) >= 4) {
  subset = paste('_', args[4], sep="")
}

BASE_DIR = Sys.getenv('BASE_DIR')

contrast_file <- paste(BASE_DIR, RP_id, '/analysis/input/metadata/xtail_contrasts', subset, '.tsv', sep="")
contrasts_df <- as.data.frame(fread(contrast_file, sep="\t"))

RP_dir = paste(BASE_DIR, RP_id, '/analysis/output/ext_diricore/', bam_type, '/tsv', sep="")
RNA_dir = paste(BASE_DIR, RNA_id, '/analysis/output/ext_diricore/', bam_type, '/tsv', sep="")

gene_file = paste(BASE_DIR, '/static/hg19/gene_names.txt', sep="")
gene_df = as.data.frame(fread(gene_file, sep='|', header = F))
colnames(gene_df) <- c('trans', 'gene_id', 'gene')
gene_df$gene <- NULL

# RP counts
rp_df = NULL
samplenames <- c()
for (f in list.files(path=RP_dir, pattern = "*.tsv", full.names = T)) {
  sample <- tools::file_path_sans_ext(basename(f))
  samplenames = c(samplenames, sample)
  df = as.data.frame(fread(f))
  df$V2 <- NULL
  df$V3 <- NULL
  colnames(df) <- c('trans')
  df$sample <- 1
  df <- merge(df, gene_df)
  df$trans <- NULL
  df <- aggregate(sample ~ gene_id, df, sum)
  df[[sample]] <- df$sample
  df$sample <- NULL
  df$gene_id <- sapply(strsplit(df$gene_id,".", fixed = T), `[`, 1)
  if (is.null(rp_df)) {
    rp_df <- df
  } else {
    rp_df <- merge(rp_df, df, all=T)
  }
}
rp_df[is.na(rp_df)] <- 0
rownames(rp_df) <- rp_df$gene_id
rp_df$gene_id <- NULL


# RNA-seq counts
rna_df = NULL
for (f in list.files(path=RNA_dir, pattern = "*.tsv", full.names = T)) {
  sample <- tools::file_path_sans_ext(basename(f))
  df = as.data.frame(fread(f))
  df$V2 <- NULL
  df$V3 <- NULL
  colnames(df) <- c('trans')
  df$sample <- 1
  df <- merge(df, gene_df)
  df$trans <- NULL
  df <- aggregate(sample ~ gene_id, df, sum)
  df[[sample]] <- df$sample
  df$sample <- NULL
  df$gene_id <- sapply(strsplit(df$gene_id,".", fixed = T), `[`, 1)
  if (is.null(rna_df)) {
    rna_df <- df
  } else {
    rna_df <- merge(rna_df, df, all=T)
  }
}

rna_df[is.na(rna_df)] <- 0
rownames(rna_df) <- rna_df$gene_id
rna_df$gene_id <- NULL

for (i in 1:nrow(contrasts_df)) {
  row <- contrasts_df[i,]
  condition <- colnames(contrasts_df)
  samples <- as.character(row)
  # select samples
  rna <- rna_df[,samples]
  rp <- rp_df[,samples]
  # delete rows with too few reads
  keep_rp <- rowSums(rp[,1:length(samples)]) >= length(samples) * 2
  rp <- rp[keep_rp,]
  keep_rna <- rowSums(rna[,1:length(samples)]) >= length(samples) * 2
  rna <- rna[keep_rna,]
  
  xtail_results <- xtail(rna, rp, condition, bins=1000)
}

