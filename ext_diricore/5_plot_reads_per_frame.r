# R/3.6.2
library(data.table)
library(ggplot2)
require(gridExtra)

args = commandArgs(trailingOnly=TRUE)
project_id=args[1]

bam_type = "all_unique"
if (length(args) >= 2) {
  bam_type = args[2] # can be: hq, hq_unique, all, all_unique
}

infile = paste('/icgc/dkfzlsdf/analysis/OE0532/', project_id, '/analysis/output/ext_diricore/', bam_type, '/psites/reads_per_frame.tsv', sep='')
outdir = paste('/icgc/dkfzlsdf/analysis/OE0532/', project_id, '/analysis/output/figures/ext_diricore/', bam_type, '/reads_per_frame/', sep='')

dir.create(outdir, recursive=T, showWarnings = F)

df = as.data.frame(fread(infile, header=T, sep="\t"))

samples <- df$sample

for (sample in samples) {
  plot_file = paste(outdir, sample, '.pdf', sep='')
  sample_df = df[df$sample == sample,]
  # get total reads per region
  utr5 = sample_df[['utr5_frame0']] + sample_df[['utr5_frame1']] + sample_df[['utr5_frame2']]
  utr3 = sample_df[['utr3_frame0']] + sample_df[['utr3_frame1']] + sample_df[['utr3_frame2']]
  cds = sample_df[['cds_frame0']] + sample_df[['cds_frame1']] + sample_df[['cds_frame2']]
  
  # normalize
  sample_df[['utr5_frame0']] = sample_df[['utr5_frame0']] / utr5
  sample_df[['utr5_frame1']] = sample_df[['utr5_frame1']] / utr5
  sample_df[['utr5_frame2']] = sample_df[['utr5_frame2']] / utr5
  
  sample_df[['utr3_frame0']] = sample_df[['utr3_frame0']] / utr3
  sample_df[['utr3_frame1']] = sample_df[['utr3_frame1']] / utr3
  sample_df[['utr3_frame2']] = sample_df[['utr3_frame2']] / utr3
  
  sample_df[['cds_frame0']] = sample_df[['cds_frame0']] / cds
  sample_df[['cds_frame1']] = sample_df[['cds_frame1']] / cds
  sample_df[['cds_frame2']] = sample_df[['cds_frame2']] / cds
  
  df2 <- data.frame(c(sample, sample, sample, sample, sample, sample, sample, sample, sample),
                    c(0, 1, 2, 0, 1, 2, 0, 1, 2),
                    c("5' UTR", "5' UTR", "5' UTR", 'CDS', 'CDS', 'CDS', "3' UTR", "3' UTR", "3' UTR"),
                    c(sample_df$utr5_frame0, sample_df$utr5_frame1, sample_df$utr5_frame2, 
                      sample_df$cds_frame0, sample_df$cds_frame1, sample_df$cds_frame2, 
                      sample_df$utr3_frame0, sample_df$utr3_frame1, sample_df$utr3_frame2))
  colnames(df2) <- c('sample', 'frame', 'region', 'percentage')
  
  
  ylim = round(max(df2$percentage), digits = 1) + 0.1
  utr5_plot <- ggplot(df2[df2$region == "5' UTR",], aes(x = frame, y = percentage)) + 
    geom_bar(stat = "identity") + theme_bw(base_size = 20) + 
    labs(x = "Frame", y = "P-site signal (%)", subtitle ="5' UTR", title='') +
    coord_cartesian(ylim=c(0, ylim))
  cds_plot <- ggplot(df2[df2$region == 'CDS',], aes(x = frame, y = percentage)) + 
    geom_bar(stat = "identity") + theme_bw(base_size = 20) + 
    labs(x = "Frame", y = "", subtitle = 'CDS', title=sample) +
    coord_cartesian(ylim=c(0, ylim))
  utr3_plot = ggplot(df2[df2$region == "3' UTR",], aes(x = frame, y = percentage)) + 
    geom_bar(stat = "identity") + theme_bw(base_size = 20) + 
    labs(x = "Frame", y='', subtitle="3' UTR", title='') +
    coord_cartesian(ylim=c(0, ylim))
    
  plot <- ggpubr::ggarrange(utr5_plot, cds_plot, utr3_plot, nrow=1, ncol=3)
  print(paste('Writing plot: ', plot_file))
  ggpubr::ggexport(plot, filename=plot_file)
}
