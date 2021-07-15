.libPaths(c('/omics/groups/OE0532/internal/RStudio_3.6.2', .libPaths()))

library(ggplot2)
library(data.table)

BASE_DIR = Sys.getenv('BASE_DIR')

## EXAMPLE: Rscript highlight_genes.r 21221_RNA /icgc/dkfzlsdf/analysis/OE0532/21221_RNA/analysis/input/metadata/te_genes.txt

project_id = '21221'
bam_type = 'all_unique'
genes_file = '/omics/groups/OE0532/internal/from_snapshot//21221_RNA/analysis/input/metadata/RPG_small.txt'
genes_file = '/omics/groups/OE0532/internal/from_snapshot//21221_RNA/analysis/input/metadata/RPG_large.txt'

args = commandArgs(trailingOnly=TRUE)

project_id = toString(args[1])
bam_type = args[2]
genes_file = args[3]

subset = basename(genes_file)
subset = gsub('.txt', '', subset)

genes_df = as.data.frame(fread(genes_file, header = F))
colnames(genes_df) <- c('gene')

indir = paste(BASE_DIR, project_id, '/analysis/output/ribo_diff/te/', bam_type, sep="")
outdir = paste(indir, subset, sep="/")
dir.create(outdir, recursive=T, showWarnings = F)

for (f in list.files(path=indir, pattern = "*.tsv", full.names = T)) {
  fn = basename(f)
  contrast = sub('_plot_data', '', tools::file_path_sans_ext(basename(f)))
  df = as.data.frame(fread(f))
  df$log2RNA = log2(df$cntRnaMean)
  df <- df[df$log2RNA > 0,]
  df$log2FC_TE <- df$logFoldChangeTE
  df$logFoldChangeTE <- NULL
  # p-value threshold
  sig_df = df[df$padj <= 0.1,]
  non_sig_df = df[df$padj > 0.1,]
  # remove genes
  sig_df = sig_df[!sig_df$gene_id %in% genes_df$gene,]
  non_sig_df = non_sig_df[!non_sig_df$gene_id %in% genes_df$gene,]
  highlight_df = df[df$gene_id %in% genes_df$gene,]
  # titles
  sig_title = paste('Significant genes \n(pval <= 0.1):', length(sig_df$gene_id))
  tested_title = paste('Tested genes:', length(df$gene_id))
  selected_title = paste('Selected genes:', length(highlight_df$gene_id))
  # define type (color + legend)
  df$type = tested_title
  df[df$gene_id %in% sig_df$gene_id,]$type <- sig_title
  df[df$gene_id %in% highlight_df$gene_id,]$type <- selected_title
  # sort by type
  df = df[order(df$type),]
  final_plot <- ggplot(aes(y=log2FC_TE, x=log2RNA), data=df[!df$gene_id %in% highlight_df$gene_id,], label=gene_id) + 
    geom_point(aes(colour=type), size=1) +
    #geom_text(data=highlight_df, aes(y=log2FC_TE, x=log2RNA, label=gene_id, hjust=0, vjust=0)) + 
    geom_point(aes(y=log2FC_TE, x=log2RNA), data=highlight_df, size=2, color='red') +
    theme(legend.title=element_blank()) +
    scale_color_manual(values=c("orange", "darkgrey", "red")) +
    labs(title = paste('Translation Efficiency Change,', contrast),
         subtitle = paste('Genes: ', subset, " (", length(genes_df$gene), ")", sep=""))
  # output
  plot_file = paste(outdir, sub("_plot_data.tsv", ".scatter-highlight.pdf", fn), sep="/")
  print(paste('Writing file:', plot_file))
  ggsave(plot_file, final_plot, width=10, height=6)
  outfile = paste(outdir,  sub("_plot_data.tsv", ".scatter-highlight.tsv", fn), sep="/")
  print(paste('Writing file:', outfile))
  write.table(highlight_df, quote = F, sep="\t", row.names = F, file=outfile)
}
