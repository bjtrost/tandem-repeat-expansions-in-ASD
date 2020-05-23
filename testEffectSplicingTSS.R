## this script is to test that expansions with potential disruption of TSS or splicing, 
## impact genes with higher selection pressure than other genes
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(cowplot)

load("gs.asd.2019.RData")
annotation <- fread("merged.expansion.coordinates_forAnno.txt.annovar.out_rev38.1mFmt.tsv", stringsAsFactors = F)
annotation <- data.frame(annotation)
annotation$key <- paste(annotation$chr, annotation$start, annotation$end, sep="#")
annotation$typeseq_priority <- factor(annotation$typeseq_priority, 
                                      levels = c("exonic", "splicing", "exonic;splicing", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_exonic;ncRNA_splicing", 
                                                 "UTR5", "UTR3", "intronic", "ncRNA_intronic", "upstream", "downstream", "intergenic"))
annotation <- annotation[order(annotation$typeseq_priority), ]
annotation <- annotation[!duplicated(annotation$key), ]

tss <- read.delim("tss.candidate.expansion.tsv", stringsAsFactors = F)
splicing <- read.delim("splicing.candidate.expansion.tsv", stringsAsFactors = F)

gnomad.constraint <- read.delim("gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

tss <- merge(tss, annotation[, c("key", "gene_symbol")], by.x = "repeatID", by.y = "key", all.x = T)
splicing <- merge(splicing, annotation[, c("key", "gene_symbol")], by.x = "repeatID", by.y = "key", all.x = T)

gnomad.constraint <- gnomad.constraint[order(gnomad.constraint$oe_lof_upper, decreasing = T), ]
gnomad.constraint <- gnomad.constraint[!duplicated(gnomad.constraint$gene), ]

tss.constraint <- gnomad.constraint$oe_lof_upper[gnomad.constraint$gene %in% tss$gene_symbol]
nontss.constraint <- gnomad.constraint$oe_lof_upper[!gnomad.constraint$gene %in% tss$gene_symbol]

splicing.constraint <- gnomad.constraint$oe_lof_upper[gnomad.constraint$gene %in% splicing$gene_symbol]
nonsplicing.constraint <- gnomad.constraint$oe_lof_upper[!gnomad.constraint$gene %in% splicing$gene_symbol]
non.constraint <-  gnomad.constraint$oe_lof_upper[(!gnomad.constraint$gene %in% splicing$gene_symbol) &
                                                    (!gnomad.constraint$gene %in% tss$gene_symbol)]

wilcox.test(tss.constraint, non.constraint, alternative = "less") #2e-5
wilcox.test(splicing.constraint, non.constraint, alternative = "less") #1.6e-6

boxplot(tss.constraint, splicing.constraint, non.constraint,
        names = c("TSS impacted", "Splice junction impacted", "Other genes"), ylim = c(-0.4, 2),
        ylab = "gnomAD o/e upper bound")
lines(x = c(2, 3), y = c(-0.2, -0.2))
lines(x = c(1, 3), y = c(-0.4, -0.4))
text(x = 2.5, y = -0.1, labels = parse(text=c("Wilcoxon~test~p==1.6~x~10^-6")), cex = 0.8)
text(x = 2, y = -0.3, labels = parse(text=c("Wilcoxon~test~p==2~x~10^-5")), cex = 0.8)

