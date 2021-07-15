library(ggplot2)
library(data.table)
library(RColorBrewer)
library(tools)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
project_id=args[1]
bam_type = args[2]
x = args[3]
site = args[4] # 12: psite, 15: asite, 18: esite

# all by default
aa = c("ala","arg","asn","asp","cys","gln","glu","gly","his","ile","leu","lys","met","phe","pro","ser","stp","thr","trp","tyr","val") 
if (length(args) >= 5) {
  aa = strsplit(args[5], ",") # amino acid(s) in format: pro,asp,leu  
}

# pretty purple: '#ff8c66',
colors = c("#cc00cc", '#668cff',  '#00cc99', '#262626', '#b366ff', '#ff6666', '#66ffff', '#006666', '#558000', '#005580', '#862d59', '#4d0000', '#ff0000', '#4d0066', '#663300', '#006666', '#669999', '#ff99ff', '#f58231', '#000075', '#a9a9a9')
# '#b366ff',
colors = c('#ff6666', '#cc00cc',
           '#911eb4',
           '#4363d8', '#808080', '#3cb44b', 
           '#ffe119',  
           '#f58231', '#46f0f0', '#bcf60c', 
           '#fabebe', '#008080', '#e6beff', 
           '#9a6324', '#fffac8', '#800000', 
           '#aaffc3', '#808000', '#ffd8b1', 
           '#000075', '#00cc99', 
           '#ff8c66')

indir = paste("/icgc/dkfzlsdf/analysis/OE0532/", project_id, "/analysis/output/ext_diricore/", bam_type, "/codons/", x, "x", sep="")
plot_dir = paste("/icgc/dkfzlsdf/analysis/OE0532", project_id, "analysis/output/figures/ext_diricore", bam_type, paste(x, "x", sep=""), sep="/")

alignment_stats = paste("/icgc/dkfzlsdf/analysis/OE0532/", project_id, "/analysis/output/alignment_stats/", bam_type, "_stats.txt", sep="")
reads_df = as.data.frame(read.csv(alignment_stats, sep="\t", header=F, col.names = c("sample", "reads")))

samplenames = c()
full_df <- NULL
for (filename in list.files(indir)) {
  infile = paste(indir, filename, sep="/")
  samplename = file_path_sans_ext(filename)
  samplename = gsub(paste(site, "_*", sep=""), "\\1", samplename)
  if (startsWith(filename, toString(site))) {
    samplenames <- c(samplenames, samplename)
    df <- as.data.frame(read.csv(infile, sep="\t", header=F, col.names = c('codon', 'counts', 'aa', 'norm_counts')))
    df <- df[df$aa %in% aa,]
    df$sample <- samplename
    # df$counts <- NULL
    # df$aa <- NULL
    # colnames(df) = c('codon', samplename)
    if (is.null(full_df)) {
      full_df <- df
    } else {
      full_df <- rbind(full_df, df)
      # merge(x = df1, y = df2, by = "CustomerId", all = TRUE)
      # full_df <- merge(x=full_df, y=df, by='codon', all=TRUE)
    }
    # frame_reads = sum(df$counts)
    # sample_reads = reads_df[reads_df$sample==samplename,]$reads[1]
    # df$counts <- round((df$counts / frame_reads) * 100, 2)
    # df$codon <- factor(df$codon, levels=df$codon)
    # wid = max(length(df$codon) * 0.3, 3)
    # p1 <- ggplot(data=df, aes(x=codon, y=counts, fill=aa)) +
    #   geom_bar(stat="identity", width=0.75) + #scale_fill_manual(values=colors, name="amino acids") +
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", legend.direction="horizontal") +
    #   scale_fill_discrete(name = "amino acids") + ylab('counts')
    # 
    # p2 <- ggplot(data=df, aes(x=codon, y=norm_counts, fill=aa)) +
    #   geom_bar(stat="identity", width=0.75) + #scale_fill_manual(values=colors, name="amino acids") +
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1),  legend.position="none") +
    #   scale_fill_discrete(name = "amino acids") + ylab('norm counts')
    # n_rows = 2
    # height = 14
    # if (length(df$codon) <= 20) {
    #   n_rows = 1
    #   height = 7
    # }
    # title = paste(samplename, " (pos", site, ")\nFrame reads: ", frame_reads, "\nTotal sample reads: ", sample_reads, sep="")
    # p = grid.arrange(p1, p2, nrow=n_rows, top=title)
    # 
    # 
    # imgfile = paste(plot_dir, '/', aa, '.', x, "x.", "site", site, ".", samplename, ".", bam_type, ".jpeg", sep="")
    # print(paste("Saving plot:", imgfile))
    # ggsave(imgfile, device = NULL, plot=p, width=wid, height=height, units="in")
  }
}

# full_df$counts <- round((full_df$counts / frame_reads) * 100, 2)
# full_df$codon <- factor(full_df$codon, levels=full_df$codon)
n_codons = length(unique(full_df$codon[1]))
wid = max(n_codons * 0.3 * length(samplenames), 3)


factor(full_df$codon,levels=samplenames,ordered=FALSE)


df1 <- melt(full_df, id.vars='codon')
ggplot(data.m, aes(codon, value)) +   
  geom_bar(aes(fill = variable), position = "dodge", stat="identity")

colnames(df1) <- c('codon', 'sample', 'norm_counts')



ggplot(data=df1, aes(x=codon, y=norm_counts, fill=sample)) +
  geom_bar(stat="identity", width=0.75, aes(fill=sample)) + #scale_fill_manual(values=colors, name="amino acids") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom", legend.direction="horizontal") +
  scale_fill_discrete(name = "amino acids") + ylab('norm counts')


p1 <- ggplot(data=full_df, aes(x=codon, y=counts, fill=aa)) +
  geom_bar(stat="identity", width=0.75) + #scale_fill_manual(values=colors, name="amino acids") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", legend.direction="horizontal") +
  scale_fill_discrete(name = "amino acids") + ylab('counts')

p2 <- ggplot(data=full_df, aes(x=codon, y=norm_counts, fill=aa)) +
  geom_bar(stat="identity", width=0.75) + #scale_fill_manual(values=colors, name="amino acids") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  legend.position="none") +
  scale_fill_discrete(name = "amino acids") + ylab('norm counts')
n_rows = 2
height = 14
if (length(full_df$codon) <= 20) {
  n_rows = 1
  height = 7
}
title = paste(samplename, " (pos", site, ")\nFrame reads: ", frame_reads, "\nTotal sample reads: ", sample_reads, sep="")
p = grid.arrange(p1, p2, nrow=n_rows, top=title)


imgfile = paste(plot_dir, '/', aa, '.', x, "x.", "site", site, ".", samplename, ".", bam_type, ".jpeg", sep="")
print(paste("Saving plot:", imgfile))
ggsave(imgfile, device = NULL, plot=p, width=wid, height=height, units="in")
