#!/usr/bin/env Rscript

# R/3.5.1
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(reshape2)

# get all_stats.txt file:
## 1. Get alignment stats: `bsub -q medium -R "rusage[mem=30G]" /icgc/dkfzlsdf/analysis/OE0532/software/diricore/get_alignment_stats.sh 20910`
## 2. Aggregate stats into 1 file: `python /icgc/dkfzlsdf/analysis/OE0532/software/diricore/get_alignment_stats_2.py 20910`


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop(paste("Usage example:", args[0], "14548"))
} else if (length(args)>=1) {
  project_id=toString(args[1])
}

subset = ""
if (length(args) >= 2) {
  subset = paste(args[2], "_", sep="")
}
BASE_DIR = Sys.getenv('BASE_DIR')
output_file = paste(BASE_DIR, project_id, "/analysis/output/figures/Alignment_stats.pdf", sep="")
stats_file = paste(BASE_DIR, project_id, "/analysis/output/alignment_stats/all_stats.txt", sep="")
df = fread(stats_file, sep="\t", fill=T)
colnames(df) <- c("sample", "bc_reads", "rrna_reads", "align_reads", "unique_reads")
df$bc_reads <- NULL
df$duplicated_reads = df$align_reads - df$unique_reads
df$align_reads <- NULL
df <- df[,c('sample', 'rrna_reads', 'duplicated_reads', 'unique_reads')]

df2 <- reshape2::melt(df,id.vars='sample')
df2$sample <- factor(df2$sample, levels = rev(c(df$sample, df$variable))) # to keep original order


# Plot
print(paste("Plot: ", output_file, sep=""))
myplot <- ggplot(df2, aes(x=sample, y=value/1000000,fill=variable)) + 
  geom_bar(stat='identity')+ ggtitle(paste("Dataset: ", project_id)) +
  theme_bw(8) + coord_flip() + ylab('Million reads') + xlab(NULL)

ggsave(output_file, myplot,width = 5, height = 3)
myplot

