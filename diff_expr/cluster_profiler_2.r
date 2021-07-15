.libPaths(c('/omics/groups/OE0532/internal/RStudio_3.6.2', .libPaths()))

library(stats)
library(data.table)
library(org.Hs.eg.db) # BiocManager::install('org.Hs.eg.db')
# library(org.Mm.eg.db) # BiocManager::install('org.Mm.eg.db')
library(DOSE)         # BiocManager::install('DOSE')
library(enrichplot)   # BiocManager::install('enrichplot')
library(ggplot2)
library(clusterProfiler)  # BiocManager::install('clusterProfiler')
library(ggpubr)
library(gridExtra)

# R 3.6.2
args = commandArgs(trailingOnly=TRUE)
project_id=args[1]
bam_type = args[2]

min_reads = 20
min_reads = 50
min_reads = 100
min_reads = as.numeric(args[3])

threshold = 2
threshold = -2
threshold = 1.5
threshold = -1.5
threshold = as.numeric(args[4])


genome = 'mm10'
genome = 'hg19'
if (length(args) >= 5) {
  genome = args[5]
}
if (genome == 'mm10') {
  org_db = 'org.Mm.eg.db'
  org_db_var = org.Mm.eg.db
} else {
  org_db = 'org.Hs.eg.db'
  org_db_var = org.Hs.eg.db
}

prefix = ''
subset = ''
if (length(args) >= 6) {
  subset = args[6]
  prefix = paste('_', subset, sep='')
}
# before: get diff_expr files
# 

BASE_DIR = Sys.getenv('BASE_DIR')

# genes_dir <- paste(BASE_DIR, project_id, "/analysis/output/cluster_profiler/", sep="")
plot_dir <- paste(BASE_DIR, project_id, "/analysis/output/figures/go_term/", bam_type, "/reads", min_reads, "_fc", threshold, "/", subset, "/", sep="")
data_dir <- plot_dir # paste(BASE_DIR, project_id, "/analysis/output/go_term/", bam_type, "/reads", min_reads, "_fc", threshold, "/", subset, "/", sep="")
diff_expr_file = paste(BASE_DIR, project_id, "/analysis/output/diff_expr/", bam_type, "_reads", min_reads, "/FC_diff_expr.tsv", sep="")

de_results = as.data.frame(fread(diff_expr_file, sep="\t", header=T))

## to specify order of columns

contrasts_file <- paste(BASE_DIR, project_id, "/analysis/input/metadata/rpf_density_contrasts", prefix, ".tsv", sep="")
contrasts_dt <- fread(contrasts_file, header=F, sep="\t")
colnames(contrasts_dt) <- c("sample", "control", "color")
contrasts_dt$color <- NULL

all_contrasts <- c()
for (i in 1:nrow(contrasts_dt)) {
  row <- contrasts_dt[i,]
  sample = row$sample
  ctrl = row$control
  contrast <- paste(sample, ctrl, sep="__vs__")
  if (contrast %in% colnames(de_results)) {
    all_contrasts <- c(all_contrasts, contrast) 
  }
}

# all_contrasts <- colnames(de_results)
# all_contrasts <- all_contrasts[2:length(all_contrasts)] # remove 'gene' column

gene_names = paste(BASE_DIR, "/static/", genome, "/gene_names.txt", sep="")
gene_names_dt <- as.data.frame(fread(gene_names, sep="|", header=F))
colnames(gene_names_dt) <- c('trans', 'gene_id', 'gene')

## get a table with genes in each GO term category
dir.create(plot_dir, recursive = T, showWarnings = F)
dir.create(data_dir, recursive = T, showWarnings = F)


for (contrast in all_contrasts) {
  print(contrast)
  selected_de <- de_results[,c('gene', contrast)]
  if (threshold < 0) {
 	 selected_de <- selected_de[selected_de[[contrast]] <= threshold,]
  } else {
     selected_de <- selected_de[selected_de[[contrast]] >= threshold,]
  }
  gene_list <- unique(selected_de$gene)
  if (length(gene_list) == 0) {
    print(paste('No genes for contrast', contrast))
  }  else {
    entrez_ids <- as.data.frame(mapIds(org_db_var, gene_list, 'ENTREZID', 'SYMBOL'))
    colnames(entrez_ids) <- c('entrez_id')
    entrez_ids$gene <- rownames(entrez_ids)
    
    selected_de <- merge(selected_de, entrez_ids)
    selected_de <- selected_de[!is.na(selected_de$entrez_id),]
    
    gene_ids <- as.vector(selected_de$entrez_id)
    
    # how to make a custom geneList: 
    # https://github.com/YuLab-SMU/DOSE/wiki/how-to-prepare-your-own-geneList
    geneList = as.numeric(selected_de[[contrast]])
    names(geneList) = gene_ids
    geneList = sort(geneList, decreasing = TRUE)
    
    ego <- enrichGO(gene_ids, OrgDb=org_db, pvalueCutoff=0.05, qvalueCutoff=0.05, ont="BP", readable=TRUE)
    # ego <- enrichGO(gene_ids, OrgDb=org_db, pvalueCutoff=0.05, ont="BP")
    if (is.null(ego)) {
      print(paste('No GO terms for contrast', contrast))
    } else {
      ego_df <- as.data.frame(ego)
      ego_df <- ego_df[,c('ID', 'Description', 'Count', 'geneID')]
      outfile = paste(data_dir, contrast, ".tsv", sep="")
      print(paste('Writing file:', outfile))
      write.table(ego_df, outfile, row.names = F, quote = F)
      
      # 12.2 Dot plot
      p1 <- dotplot(ego, showCategory=30) + ggtitle(gsub('__', ' ', contrast))
      outplot <- paste(plot_dir, "/", contrast, ".pdf", sep="")
      print(paste('12.2 Dot plot. Saving plot:', outplot))
      ggsave(outplot, p1, width = max(7, max(nchar(ego$Description)) / 10 + 3), height = 7)
    }
  }
}


# Single plot for all comparisons
# 13.2 https://bioconductor.statistik.tu-dortmund.de/packages/3.3/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.htm
# PS: instead of "plot" use "dotplot" (not like in the book), otherwise - stupid error

full_df <- NULL
groups <- c()
all_ids <- c()
clusters = c()

for (i in 1:length(all_contrasts)) {
  contrast = all_contrasts[i];
  print(contrast)
  
  clusters <- c(clusters, paste('C', i, sep=""))
  selected_de <- de_results[,c('gene', contrast)]
  # selected_de <- selected_de[abs(selected_de[[contrast]]) >= threshold,]
  if (threshold < 0) {
    selected_de <- selected_de[selected_de[[contrast]] <= threshold,]
  } else {
    selected_de <- selected_de[selected_de[[contrast]] >= threshold,]
  }
  gene_list <- unique(selected_de$gene)
  if (length(gene_list) == 0) {
    print(paste('No genes for contrast', contrast))
  }  else {
    entrez_ids <- as.data.frame(mapIds(org_db_var, gene_list, 'ENTREZID', 'SYMBOL'))
    colnames(entrez_ids) <- c('entrez_id')
    entrez_ids$gene <- rownames(entrez_ids)
    
    selected_de <- merge(selected_de, entrez_ids)
    selected_de <- selected_de[!is.na(selected_de$entrez_id),]
    
    gene_ids <- as.vector(selected_de$entrez_id)
    
    # geneList = as.numeric(selected_de[[contrast]])
    # names(geneList) = gene_ids
    # geneList = sort(geneList, decreasing = TRUE)
    
    full_df[[paste('C', i, sep="")]] <- gene_ids
    groups <- c(groups, replicate(length(gene_ids), paste('C', i, sep="")))
    all_ids <- c(all_ids, gene_ids)
  }
}


mydf <- data.frame(Entrez=all_ids, group = groups)
xx.formula <- compareCluster(Entrez~group, data=mydf, fun='enrichGO', OrgDb=org_db, ont='BP', pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE)

# go_df1 <- as.data.frame(xx.formula)
# go_df1 <- go_df1[!(go_df1$ID %in% c('GO:0051607', 'GO:0009615', 'GO:0034341')),]
# dotplot(go_df1) 

go_plot1 <- dotplot(xx.formula, showCategory=30) + 
  ggtitle(paste('Cluter Profiler. Dataset', project_id, ' ', bam_type, '. Min reads: ', min_reads, sep=""))
outplot <- paste(plot_dir, "/dotplot_reads",  min_reads, "_fc", threshold, "_", prefix, ".pdf", sep="")
print(paste('Saving plot:', outplot))
ggsave(outplot, go_plot1, width=16, height = 20)

labels <- data.frame(label=clusters,contrast = all_contrasts)
labels_plot <- tableGrob(labels)
labels_file <- paste(plot_dir, prefix, "_annotation_all_comparisons.jpeg", sep="")
print(paste('Saving plot annotation:', labels_file))
ggsave(labels_file, labels_plot)



go_df <- as.data.frame(xx.formula)
# replace entrez_id with gene name
# go_df$gene <- ""
# for(i in 1:nrow(go_df)) {
#   genes = go_df[i,]$geneID
#   genes <- strsplit(genes, "/")[[1]]
#   gene_names <- as.data.frame(mapIds(org_db_var, genes, 'SYMBOL', 'ENTREZID'))
#   colnames(gene_names) <- c('name')
#   go_df[i,]$gene <- paste(as.vector(gene_names$name), collapse = '/')
# } 
# go_df$geneID <- NULL

go_df$contrast = ""
for (i in 1:length(all_contrasts)) {
  contrast = all_contrasts[i]
  cluster = paste('C', i, sep="")
  go_df1 <- go_df[go_df$Cluster == cluster,]
  if (nrow(go_df1) > 0) {
    go_df[go_df$Cluster == cluster,]$contrast = contrast
  }
}

go_df <- go_df[,c('contrast', 'ID', 'Description', 'Count', 'GeneRatio', 'p.adjust',  'geneID')]

outfile = paste(data_dir, "/go_all_comparisons", prefix, ".tsv", sep="")
print(paste('Writing file:', outfile))
write.table(go_df, outfile, row.names = F, quote = F, sep="\t")

