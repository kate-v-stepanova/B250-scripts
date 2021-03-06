#!/usr/bin/env Rscript

# R/3.5.1
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
project_id = ''
if (length(args)==0) {
  stop(paste("Usage example:", args[0], "14548"))
} else if (length(args)>=1) {
   project_id=toString(args[1])
}

subset = ""
if (length(args) >= 2) {
  subset = paste(args[2], "_", sep="")
}

BASE_DIR=Sys.getenv('BASE_DIR')
PROJECT_DIR=paste(BASE_DIR, project_id, sep="/")

INDIR=paste(PROJECT_DIR, "analysis/output/alignment_stats", sep="/")
OUTDIR=paste(PROJECT_DIR, "analysis/output/figures", sep="/")
BC_SPLIT_FILE=paste(PROJECT_DIR, "/analysis/output/", subset, "bc_split_stats.txt", sep="")
cutadapt_file=paste(PROJECT_DIR, "/analysis/output/", subset, "cutadapt_plot_stats.txt", sep="")
# hq_file = paste(INDIR, "hq_stats.txt", sep="/")
# hq_unique = paste(INDIR, "hq_unique_stats.txt", sep="/")
# all_stats = paste(INDIR, "all_stats.txt", sep="/")
# all_unique = paste(INDIR, "all_unique_stats.txt", sep="/")
# rrna_file = paste(INDIR, "rrna_stats.txt", sep="/")
# trna_file = paste(INDIR, "trna_stats.txt", sep="/")

dir.create(OUTDIR, recursive=T)

print("Parsing cutadapt_plot_stats.txt")
if (file.exists(cutadapt_file)) {
  lines <- read.csv(file=cutadapt_file, sep="\t", colClasses=c("NULL", NA), header=F)
  no_adapt = lines[1,] - lines[2,]
  too_short = strtoi(lines[3,])
  passed = strtoi(lines[4,])
  mydt <- data.table(
    'Reads' = c('No_adapt','Too_short','Passed'),
    'Counts' = c(no_adapt, too_short, passed))
  
  myplot <- ggplot(mydt,aes(x=project_id,y=Counts/1000000, fill=Reads))
  myplot <- myplot + geom_bar(stat='identity',position='stack')+
    theme_bw()+xlab(NULL)+ ylab('Milion Reads')
  if (subset != '') {
    myplot <- myplot + ggtitle(gsub('_', ' ', subset))
  }
  print(paste('Saving ', OUTDIR, '/', subset, 'Cutadapt_stats.pdf', sep=""))
  ggsave(paste(OUTDIR, '/', subset, 'Cutadapt_stats.pdf', sep=""), myplot, width = 2.5, height = 4)
} else {
  print(paste('File does not exist:', cutadapt_file))
}

# BC stats
print("BC stats")
if (file.exists(BC_SPLIT_FILE)) {
  
  mydt <- fread(BC_SPLIT_FILE, sep="\t", fill=T)
  mydt <-mydt[!(mydt$Barcode=="total"),]
  mydt$Barcode <- factor(mydt$Barcode, levels = mydt$Barcode) # to keep original order
  myplot <- ggplot(mydt, aes(x=project_id, y=Count/1000000, fill=Barcode)) + 
    geom_bar(stat='identity',position='stack', colour='white')+
    scale_fill_manual(values=c(brewer.pal(12,'Set3'), brewer.pal(8,'Accent'), brewer.pal(9,'Set1'), brewer.pal(8,'Dark2')))  +#'gray30', '#00cc99', '#ff6666', '#669999', '#006666'))+
    xlab(NULL) + ylab('Million reads') + theme_bw()
  # If Get an error:
  ### Error: Insufficient values in manual scale. 17 needed but only 14 provided.
  # then add some colors manually:  scale_fill_manual(values=c(brewer.pal(12,'Set3'),'gray30', '#00cc99', '#cc00cc', #862d59', '#006666'))+
  if (subset != '') {
    myplot <- myplot + ggtitle(gsub('_', ' ', subset))
  }
  print(paste("Saving ", OUTDIR, '/', subset, 'BCsplit_stats.pdf', sep=""))
  ggsave(paste(OUTDIR, '/', subset, 'BCsplit_stats.pdf', sep=""), myplot, width = 3.5, height = 4.5) 
} else {
  print(paste('File does not exist:', BC_SPLIT_FILE))
}

# samples = mydt$Barcode
# samples <- samples[!grepl('unmatched$', samples)]
# mydt$Location <- NULL
# # mydt$sample <- mydt$Barcode
# # mydt$Barcode <- NULL
# colnames(mydt) <- c('sample', 'bc_split')
# # Alignment stats
# # rRNA stats
# print("rRNA stats")
# if(!file.exists(rrna_file)) {
#   mydt$rrna <- 0 
# } else {
#   mydt <- merge(mydt, fread(rrna_file, col.names = c('sample', 'rrna')), by='sample')
# }
# 
# # tRNA
# print("tRNA stats")
#   if(!file.exists(trna_file)) {
#   mydt$trna <- 0 
# } else {
#   dt <- fread(trna_file, col.names = c('sample', 'trna'))
#   mydt <- merge(mydt,dt, by='sample')
# }

# # HQ unique
# print("HQ unique")
# if(!file.exists(hq_unique)) {
#   mydt$hq_unique <- 0 
# } else {
#   dt <- fread(hq_unique, col.names=c('sample', 'hq_unique'))
#   mydt <- merge(mydt,dt, by='sample')
# }
# 
# # HQ with dup
# print("HQ stats")
# if(!file.exists(hq_file)) {
#   mydt$hq_with_dup <- 0 
# } else {
#   dt <- fread(hq_file, col.names=c('sample', 'hq_with_dup'))
#   mydt <- merge(mydt, dt, by="sample")
#   mydt$hq_with_dup <- mydt$hq_with_dup - mydt$hq_unique # hq_with_dup
# }
# 
# 
# # All unique
# print("All unique")
# if(!file.exists(all_unique)) {
#   mydt$all_unique <- 0 
# } else {
#   dt <- fread(all_unique, col.names = c('sample', 'all_unique'))
#   mydt <- merge(mydt,dt, by='sample')
#   mydt$lq_unique <- mydt$all_unique - mydt$hq_unique
#   mydt$all_unique <- NULL
# }
# 
# # LQ with dup
# print("LQ stats")
# if(!file.exists(all_stats)) {
#   mydt$all <- 0 
#   mydt$lq_with_dup <- 0
# } else {
#   dt <- fread(all_stats, col.names=c('sample', 'all'))
#   mydt <- merge(mydt,dt, by='sample')
#   mydt$lq_with_dup <- mydt$all - mydt$hq_with_dup - mydt$hq_unique - mydt$lq_unique
#   mydt$all <- NULL
# }
# 
# # # total reads
# # print("total reads")
# # dt <- fread(BC_SPLIT_FILE, sep="\t", fill=T, col.names = c('sample', 'bc_split', 'file'))
# # dt$file <- NULL
# # mydt <- merge(mydt, dt, by='sample')
# 
# col_order = c("sample", "rrna", "trna", "lq_with_dup", "lq_unique", "hq_with_dup", "hq_unique", "bc_split")
# mydt <- as.data.frame(mydt)[, col_order]
# rownames(mydt) <- mydt$sample
# mydt <- as.data.frame(mydt)[samples,]
# 
# # Create stats
# print("Create diricore stats")
# write.table(mydt, paste(INDIR, '/diricore_stats.txt',sep=""), sep='\t', quote=F, row.names=F)
# 
# # Check
# rowSums(mydt[,-1]) == mydt$bc_split # True
# mydt$bc_split <- NULL
# # Arrange and save
# print("Arrange and save")
# mysum <- reshape2::melt(mydt,id.vars='sample')
# mysum$sample <- factor(mysum$sample, levels = c(mysum$sample, mysum$variable)) # to keep original order
# 
# # Plot
# print(paste("Plot: ", OUTDIR, "/Diricore_stats.pdf", sep=""))
# myplot <- ggplot(mysum, aes(x=sample, y=value/1000000,fill=variable)) + 
#       geom_bar(stat='identity')+ ggtitle("Alignment stats") +
#       theme_bw(8) + coord_flip() + ylab('Million reads') + xlab(NULL)
# ggsave(paste(OUTDIR, '/Diricore_stats.pdf', sep=""), myplot,width = 5, height = 3)
# 
