#!/usr/bin/env Rscript

# Input: Three-column TSV file of the following format:
# Sample  EHDN call count Dataset/platform/library type   Predicted ancestry
# 7-0314-003A     2321    MSSNG/Illumina HiSeq X/PCR-free EUR
# 2-1594-001      1866    MSSNG/Illumina HiSeq X/PCR-based        OTH
# 3-0325-101      2374    MSSNG/Illumina HiSeq X/PCR-free EUR

# Output: Boxplots showing the distribution of call counts for each "Dataset/platform/library type" and (for MSSNG only) "Predicted ancestry"


suppressMessages(library(BTlib))
suppressMessages(library(extrafont))
suppressMessages(library(ggplot2))
suppressMessages(library(OneR))

args = commandArgs(TRUE)
file = args[1]
outfile = args[2]

data = read.table(file, header=TRUE, check.names=FALSE, sep="\t")
data_predicted_ancestry = subset(data, `Dataset/platform/library type` == "MSSNG/Illumina HiSeq X/PCR-free")

text_size=20

p1 =    ggplot(data, aes_string(x="`Dataset/platform/library type`", y="`EHDN call count`")) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size=10)) +
        theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        theme(panel.background=element_rect(fill = "white")) +
        theme(panel.border=element_rect(color="black", fill=NA)) +
        theme(panel.grid.major.y = element_line(colour = "lightgrey", size=0.15)) +
        scale_y_continuous(breaks=seq(from=0, to=100000, by=2000)) +
        theme(text=element_text(size=text_size))

p2 =    ggplot(data_predicted_ancestry, aes_string(x="`Predicted ancestry`", y="`EHDN call count`")) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
        theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        theme(panel.background=element_rect(fill = "white")) +
        theme(panel.border=element_rect(color="black", fill=NA)) +
        theme(panel.grid.major.y = element_line(colour = "lightgrey", size=0.15)) +
        scale_y_continuous(breaks=seq(from=0, to=100000, by=2000)) +
        theme(text=element_text(size=text_size))

suppressMessages(ggsave(paste(outfile, ".platform_PCR.pdf", sep=""), plot=p1))
suppressMessages(ggsave(paste(outfile, ".predicted_ancestry.pdf", sep=""), plot=p2))
