library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(spatstat)

args = commandArgs(trailingOnly=TRUE)
project_id = args[1]
bam_type = args[2]
x = as.numeric(args[3])
site = as.numeric(args[4])
aa = args[5]
window = 0
if (length(args) >= 6) {
  window = args[6]
}
subset_name = ''
if (length(args) >= 7) {
  subset_name = paste('_', args[7], sep="")
}

contrasts = paste("/icgc/dkfzlsdf/analysis/OE0532/", project_id, "/analysis/input/metadata/rpf_density_contrasts", subset_name, ".tsv", sep="")

contrasts_df = as.data.frame(fread(contrasts, sep="\t", header = F))
colnames(contrasts_df) <- c('sample', 'ctrl', 'color')
contrasts_df$color <- NULL

BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
project_dir = paste(BASE_DIR, project_id, sep="/")
indir = paste(project_dir, "/analysis/output/ext_diricore/", bam_type, "/codons_per_pos/", x, "x.window", window, sep="")
outdir = paste(indir, subset_name, sep="/")
plot_dir = paste(project_dir, "/analysis/output/figures/ext_diricore/", bam_type, "/", x, "x.window", window, "/", subset_name, sep="")

# create output dirs
dir.create(outdir, recursive=T, showWarnings = F)
dir.create(plot_dir, recursive=T, showWarnings = F)

# get 2 groups of transcripts - with >= 1 DD sites per transcript and with 1 site per transcript
codon_file = paste("/icgc/dkfzlsdf/analysis/OE0532/static/hg19/codons/", x, "x/", aa, "_codon_positions.tsv", sep="")
codon_df = as.data.frame(fread(codon_file, sep="\t", header=F))
colnames(codon_df) <- c('transcript', 'codon', 'pos')


codon_df <- count(codon_df, c("transcript", "codon"))

group1 = codon_df[codon_df['freq'] > 1,]$transcript
group2 = codon_df[codon_df['freq'] == 1,]$transcript

# norm counts
norm_dir = paste("/icgc/dkfzlsdf/analysis/OE0532/", project_id, "/analysis/output/alignments/reads_per_gene/normalized/", bam_type, sep="")
if (!dir.exists(norm_dir) || length(list.files(norm_dir) == 0)) {
  print(paste('Normalized counts not found: ', norm_dir))
  print(paste('Get norm counts: python /icgc/dkfzlsdf/analysis/OE0532/software/scripts/normalize_counts.py', project_id, bam_type))
  stop()
} 

full_df <- NULL
for(i in 1:length(contrasts_df)) {
  row = contrasts_df[i,]
  sample = row$sample
  ctrl = row$ctrl
  contrast = paste(sample, '__vs__', ctrl, sep="")
  
  sample_file = paste(indir, "/", site, "_", sample, ".tsv", sep="")
  sample_df = as.data.frame(fread(sample_file, sep="\t"))
  
  ctrl_file = paste(indir, "/", site, "_", ctrl, ".tsv", sep="")
  ctrl_df = as.data.frame(fread(ctrl_file, sep="\t"))
  
  sample_df = sample_df[sample_df['aa'] == aa,]
  ctrl_df = ctrl_df[ctrl_df['aa'] == aa,]
  
  sample_df <- select(sample_df, c('transcript', 'codon_pos', 'counts'))
  ctrl_df <- select(ctrl_df, c('transcript', 'codon_pos', 'counts'))
  
  norm_sample = paste(norm_dir, "/", sample, '.tsv', sep="")
  norm_ctrl = paste(norm_dir,  "/", ctrl, '.tsv', sep="")
  
  norm_s_df = as.data.frame(fread(norm_sample, sep="\t"))
  norm_c_df = as.data.frame(fread(norm_ctrl, sep="\t"))
  
  norm_s_df <- select(norm_s_df, c('transcript', 'rpkm'))
  norm_c_df <- select(norm_s_df, c('transcript', 'rpkm'))
  
  sample_df = merge(sample_df, norm_s_df, by='transcript', how='inner')
  ctrl_df = merge(ctrl_df, norm_c_df, by='transcript', how='inner')
  
  sample_df['rpkm'] = sample_df['counts'] * sample_df['rpkm']
  ctrl_df['rpkm'] = ctrl_df['counts'] * ctrl_df['rpkm']
  
  colnames(sample_df) <- c('transcript', 'codon_pos', paste(sample, '_raw', sep=""), paste(sample, '_norm', sep=""))
  colnames(ctrl_df) <- c('transcript', 'codon_pos', paste(ctrl, '_raw', sep=""), paste(ctrl, '_norm', sep=""))
  
  
  df = merge(sample_df, ctrl_df, by=c('transcript', 'codon_pos'))
  
  # threshold: >= 3 reads in at least one condition
  df = df[(df[[paste(sample, '_raw', sep="")]] >= 3) & (df[paste(ctrl, '_raw', sep="")] >=3),]
  
  # calculate ratio
  df[contrast] = df[paste(sample,'_norm', sep="")] / df[paste(ctrl, '_norm', sep="")]
  df[contrast] = log2(df[[contrast]])
  df = df[order(df[[contrast]]),]
  
  
  # df1 = df[df[['transcript']] %in% group1,] # >1 DD site per transcript
  # df2 = df[df[['transcript']] %in% group2,] # =1 DD site per transcript
  df1[['log2(rpkm)']] <- df1[[contrast]]
  df2[['log2(rpkm)']] <- df2[[contrast]]
  
  df[df[['transcript']] %in% group1, 'group'] = paste(contrast, "(>1 DD site)") # >1 DD site per transcript
  df[df[['transcript']] %in% group2, 'group'] = paste(contrast, "(=1 DD site)") # =1 DD site per transcript
  
  df['log2(rpkm)'] <- df[contrast]
  df['cdf'] = ecdf(df[['log2(rpkm)']])(df[['log2(rpkm)']])
  df <- select(df, c('transcript', 'codon_pos', 'log2(rpkm)', 'cdf', 'group'))
  
  if (is.null(full_df)) {
    full_df <- df
  } else {
    full_df <- rbind(full_df,df)
  }
}

plot = ggplot(data=full_df, aes(x=`log2(rpkm)`, y=cdf, color=group)) + geom_line() +
    scale_color_discrete()
  
plot_file = paste(plot_dir, "/", aa, subset_name, ".", x, "x.site", site, ".", bam_type, ".pdf", sep="")
print(paste('Writing file:', plot_file))
ggsave(plot_file, plot, width=10)  
