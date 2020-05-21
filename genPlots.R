### This is the main script to generate figures related to statistical findings in the paper.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(data.table)
library(GenomicRanges)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(magrittr)
library(phastCons100way.UCSC.hg38)
library(ordinal)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

### figure 1e
gs <- phastCons100way.UCSC.hg38

bin.dt <- read.delim("genome.bin.1k.tsv", stringsAsFactors = F)
bin.dt <- bin.dt[bin.dt$seqnames %in% paste0("chr", c(1:22, "X", "Y")), ]
bin.g <- GRanges(bin.dt$seqnames, IRanges(bin.dt$start, bin.dt$end), "*")

tmp <- scores(gs, bin.g)

bin.dt$PhastCons <- tmp
expansion <- read.delim("merged.any.outlier.autosomal.sex.chr.nobadsamples.tsv", stringsAsFactors = F)
expansion.g <- GRanges(expansion$chr, IRanges(expansion$start, expansion$end), "*")
olap <- data.frame(findOverlaps(bin.g, expansion.g))
olap <- aggregate(subjectHits ~ queryHits, olap, length)
bin.dt$expan <- 0
bin.dt$expan[olap$queryHits] <- olap$subjectHits

phylop <- read.delim("subset_hg38.phyloP20way.wigFix.bed", header = F)
phylop.g <- GRanges(phylop$V1, IRanges(phylop$V2, phylop$V3), "*")
olap <- data.frame(findOverlaps(bin.g, phylop.g))
olap$score <- phylop$V5[olap$subjectHits]
olap <- aggregate(score ~ queryHits, olap, mean)
bin.dt$phylop <- NA
bin.dt$phylop[olap$queryHits] <- olap$score

fragile <- read.delim("fragileSites2019_revised.csv", 
                      stringsAsFactors = F)
fragile <- fragile[fragile$start.hg38. != "-", ]
fragile.g <- GRanges(fragile$chr.hg38., IRanges(as.numeric(fragile$start.hg38.), as.numeric(fragile$end.hg38.)), "*")
olap <- data.frame(findOverlaps(bin.g, fragile.g))
bin.dt$fragile <- 0
bin.dt$fragile[unique(olap$queryHits)] <- 1

known.str <- fread("simpleRepeat.txt")
known.str <- data.frame(known.str) ### 1031708
known.str.g <- GRanges(known.str$V2, IRanges(known.str$V3, known.str$V4), "*")
olap <- data.frame(findOverlaps(bin.g, known.str.g))
bin.dt$str <- 0
bin.dt$str[unique(olap$queryHits)] <- 1

bin.dt <- bin.dt[bin.dt$seqnames %in% paste0("chr", c(1:22, "X", "Y")), ]

features <- c("fragile", "GC", "str", "PhastCons", "phylop")
bin.dt[is.na(bin.dt)] <- 0
bin.dt$EHdn <- factor(bin.dt$expan > 0)
bin.dt$knownSTR <- factor(bin.dt$str > 0)
dt.out <- data.frame()
p <- list()

for(type in c("EHdn", "knownSTR"))
for(f in features){
  bin.dt$feat <- bin.dt[, f]#scale(bin.dt[, f])#
  
  lm <- glm(sprintf("%s ~ feat", type), data = bin.dt, family = binomial(link = "logit"))
  pvalue <- summary(lm)$coefficients["feat", "Pr(>|z|)"]
  
  conf <- confint.default(lm)
  dt.out <- rbind(dt.out, data.frame("feature" = f, "Odds ratio" = exp(lm$coefficients["feat"]),
                                     "OR.lower" = exp(conf["feat", 1]),
                                     "OR.upper" = exp(conf["feat", 2]),
                                     "type" = type,
                                     "pvalue" = pvalue, stringsAsFactors = F))
  
  p[[length(p) + 1]] <- ggplot(bin.dt, aes(x = feat, y = expan)) + geom_point(shape = 1) + ylab("expansions") + xlab(f) +
    geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2], lty = 2, color = "red")
}

dt.out$feature <- c("fragile sites", "GC content", "known STRs", 
                    "PhastCons", "phyloP")
dt.out$feature <- factor(dt.out$feature,
                         levels = c("GC content", "phyloP",  "PhastCons", "fragile sites", "known STRs"))

ggplot(dt.out[dt.out$feature != "known STRs", ], 
       aes(x = feature, y = Odds.ratio, fill = type)) + 
  geom_bar(stat = "identity", width = .5, position = position_dodge(width = .5), show.legend = F, color = "black") +
  geom_errorbar(aes(ymin = OR.lower, ymax = OR.upper), 
                position = position_dodge(width = .5), size = .4, width = .2, show.legend = F) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_classic() + ylab("Odds ratio") + xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        panel.border = element_rect(fill = NA)) +
    scale_fill_manual(values = c("#D6604D", "#4393C3")) +
  scale_y_continuous(breaks = c(0, 1))

######### figure 1d 
all.expansion <- read.delim("merged.EHdn.genotype.samples.without.outliers.tsv", stringsAsFactors = F)
all.g <- GRanges(all.expansion$chr, IRanges(all.expansion$start, all.expansion$end), "*")

bin.dt <- read.delim("genome.bin.1k.tsv", stringsAsFactors = F)
bin.g <- GRanges(bin.dt$seqnames, IRanges(bin.dt$start, bin.dt$end), "*")

ref.olap <- findOverlaps(all.g, bin.g)

dir <- "Definitions/"

feature.list <- list()
feature.list[["3UTR"]] <- read.delim(sprintf("%s3UTR_1b.txt", dir), stringsAsFactors = F)
feature.list[["5UTR"]] <- read.delim(sprintf("%s5UTR_1b.txt", dir), stringsAsFactors = F)
feature.list[["intergenic"]] <- read.delim(sprintf("%sintergenicExcludingCentromeresTelomeres_1b.txt", dir), stringsAsFactors = F)
feature.list[["downstream"]] <- read.delim(sprintf("%sdownstreamTxEnd_1b.txt", dir), stringsAsFactors = F)
feature.list[["upstream"]] <- read.delim(sprintf("%supstreamTxStart_1b.txt", dir), stringsAsFactors = F)
feature.list[["exon"]] <- read.delim(sprintf("%stranslatedExons_1b.txt", dir), stringsAsFactors = F)
feature.list[["intron"]] <- read.delim(sprintf("%sintrons_1b.txt", dir), stringsAsFactors = F)

dt.plot <- data.frame()
dt.all.plot <- data.frame()
all <- 0

for(i in 1:length(feature.list)){
  feat <- feature.list[[i]]
  feat.name <- names(feature.list)[i]
  feat.g <- GRanges(feat$seqnames, IRanges(feat$start, feat$end), "*")
  
  all <- all + sum(feat$width) 
  
  ### all
  olap <- data.frame(findOverlaps(feat.g, all.g))
  olap$size <- width(pintersect(feat.g[olap$queryHits], all.g[olap$subjectHits]))
  
  observed <- sum(olap$size)
  unobserved <- sum(feat$width) - observed
  proportion <- sum(olap$size)/sum(feat$width) * 100
  or <- (observed/unobserved)/(10475659/1625055396)
  p <- chisq.test(data.frame("X1" = c(observed,10475659), "X2"= c(unobserved, 1625055396)))$p.value
  dt.all.plot <- rbind(dt.all.plot, data.frame(feat.name, observed, unobserved, proportion, or, p))
}

genome.rare <- sum(width(exp.g))/all * 100
genome.all <- sum(width(all.g))/all * 100

dt.plot$feat.name <- factor(dt.plot$feat.name, levels = c("intergenic", "upstream", "5UTR", "exon", "intron", "3UTR", "downstream"))
dt.all.plot$feat.name <- gsub("UTR", "'UTR", dt.all.plot$feat.name)
dt.all.plot$feat.name <- factor(dt.all.plot$feat.name, levels = c("intergenic", "upstream", "5'UTR", "exon", "intron", "3'UTR", "downstream"))
dt.all.plot <- dt.all.plot[order(dt.all.plot$feat.name), ]

ggplot(dt.all.plot, aes(x = feat.name, y = proportion)) + 
  geom_bar(stat = "identity", fill = "black", width = .5) +
  geom_hline(yintercept = genome.all, color = "black", lty = 2) +
  ylab("proportion of features") + theme_classic() +
  xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), panel.border = element_rect(fill = NA),
                   axis.text.y = element_text(size = 14), 
                   axis.title.y = element_text(size = 14)) +
  scale_y_continuous(breaks = c(0, round(genome.all, digits = 2), 1))

library(cowplot)

#######################
######### figure 2c
rare.function <- read.delim("functional.any.motif.tsv", stringsAsFactors = F)
rare.function <- rare.function[rare.function$geneset != "intergenic", ]
rare.function$geneset <- factor(c("upstream", "downstream", "intron", "5'UTR", "3'UTR", "exon", "splicing"), 
                                levels = c("upstream", "5'UTR", "exon", "splicing", "intron", "3'UTR", "downstream"))
rare.function$fwer <- p.adjust(rare.function$all.aggr.glm.pvalue, method = "bonferroni")

ggplot(rare.function, aes(x = geneset, y = all.aggr.glm.OR, fill = all.aggr.glm.pvalue < 0.05)) + 
  geom_bar(stat="identity", show.legend = F, color = "black", width = .5) +
  geom_errorbar(aes(ymin = all.aggr.glm.lower, ymax = all.aggr.glm.upper), show.legend = F,
                position = position_dodge(width = .5), size = .4, width = .2) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_classic() + ylab("Odds Ratio") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(fill = NA)) +
  coord_cartesian(ylim = c(0, 7)) +
  scale_fill_manual(values = c("#4393C3", "#D6604D")) + 
  scale_y_continuous(breaks = c(0, 1, 3, 5)) +
  annotate("text", x = c(2.7, 4), y = c(4.7, 5), angle = 90,
           label = paste0("p==", c("2.3~x~10^-2",
                                   "2.3~x~10^-2")), parse = T,
           fontface = "italic", cex = 4)

######################### 2d
rare.expansion <- read.delim("tabulated.EHdn.outliers.overlapped.any.tsv", stringsAsFactors = F)
rare.expansion <- rare.expansion[rare.expansion$chr %in% paste0("chr", c(1:22)), ]
rare.expansion$size <- rare.expansion$end - rare.expansion$start
median.size <- median(rare.expansion$size)
iqr.size <- IQR(rare.expansion$size)

known.expansion <- read.delim("~/hpf/largeprojects/tcagstor/scratch/gpellecchia/RyanCircosPlot/Other/simpleRepeat.txt", stringsAsFactors = F, header = F)
known.expansion <- known.str[known.str$V2 %in% paste0("chr", c(1:22)), ]
detected.expansion <- read.delim("merged.any.expansions.autosomal.sex.chr.nobadsamples.tsv", stringsAsFactors = F)
detected.expansion <- detected.expansion[detected.expansion$chr %in% paste0("chr", c(1:22)), ]
detected.expansion$size <- detected.expansion$end - detected.expansion$start
boxplot(detected.expansion$size)

rare.expansion <- data.frame(reduce(GRanges(rare.expansion$chr, IRanges(rare.expansion$start, rare.expansion$end), "*")))
known.expansion <- data.frame(reduce(GRanges(known.expansion$V2, IRanges(known.expansion$V3, known.expansion$V4), "*")))
detected.expansion <- data.frame(reduce(GRanges(detected.expansion$chr, IRanges(detected.expansion$start, detected.expansion$end), "*")))

refflat <- read.delim("../data/hg38_refFlat.txt", stringsAsFactors = F, header = F)[, 1:8]
refflat$tis100kbstart <- ifelse(refflat$V4 == "+", refflat$V5-10000, refflat$V6-10000)
refflat$tis100kbend <- ifelse(refflat$V4 == "+", refflat$V5+10000, refflat$V6+10000)

refflat.g <- GRanges(refflat$V3, IRanges(refflat$tis100kbstart, refflat$tis100kbend))

common.expansion <- known.expansion
rare.expansion.g <- GRanges(rare.expansion$seqnames, IRanges(rare.expansion$start, rare.expansion$end), "*")
common.expansion.g <- GRanges(common.expansion$seqnames, IRanges(common.expansion$start, common.expansion$end), "*")
detected.expansion.g <- GRanges(detected.expansion$seqnames, IRanges(detected.expansion$start, detected.expansion$end), "*")

getdistance <- function(refflat.g, rare.expansion.g, refflat, rare.expansion){
  rare.olap <- data.frame(findOverlaps(rare.expansion.g, refflat.g))
  rare.olap$strain <- refflat$V4[rare.olap$subjectHits]
  rare.olap$tss <- ifelse(rare.olap$strain == "+", 
                          refflat$tis100kbstart[rare.olap$subjectHits] + 10000, 
                          refflat$tis100kbend[rare.olap$subjectHits] - 10000)
  
  rare.olap$expansion.start <- rare.expansion$start[rare.olap$queryHits]
  rare.olap$expansion.end <- rare.expansion$end[rare.olap$queryHits]
  rare.olap$mid.point <- (rare.olap$expansion.start + rare.olap$expansion.end) / 2
    
  rare.olap$distance <- rare.olap$tss - rare.olap$mid.point
  rare.olap$distance <- ifelse(rare.olap$strain == "+", rare.olap$distance, -rare.olap$distance)
  rare.olap <- rare.olap[order(abs(rare.olap$distance)), ]
  rare.olap <- rare.olap[!duplicated(rare.olap$queryHits), ]
  return(rare.olap)
}

rare.distance <- getdistance(refflat.g, rare.expansion.g, refflat, rare.expansion)
common.distance <- getdistance(refflat.g, common.expansion.g, refflat, common.expansion)
detected.distance <- getdistance(refflat.g, detected.expansion.g, refflat, detected.expansion)

rare.distance$rarity <- "rare expansions"
common.distance$rarity <- "known STRs"
detected.distance$rarity <- "all expansions"
  
wilcox.test(abs(rare.distance$distance), abs(detected.distance$distance), alternative = "less")$p.value #0.008
wilcox.test(abs(rare.distance$distance), abs(common.distance$distance), alternative = "less")$p.value #0.002

distance <- rbind(common.distance, rbind(rare.distance, detected.distance))

ggplot(distance, aes(x = distance, y = ..density.., color = rarity)) +
  geom_density(alpha = .15, adjust = 1/10) +
  geom_vline(xintercept = c(0), lty = 2) + xlab("Distance from TSS") + theme_classic() +
  theme(panel.border = element_rect(fill = NA), 
        legend.position = "none", 
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  geom_hline(yintercept = 0, size = 1, color = "white") +
  coord_cartesian(xlim = c(-5000, 5000)) + scale_x_continuous(breaks = c(-5000, 0, 5000)) +
  scale_y_continuous(breaks = c(0, 0.0002)) +
  scale_color_manual(values = c("#4393C3", "#338833", "#D6604D"))

###### figure 2e
refflat <- read.delim("hg38_refFlat.txt", stringsAsFactors = F, header = F)

#############################################################
introns <- read.delim("hg38_intron_refFlat.txt", stringsAsFactors = F)
appris <- rbind(read.delim("appris_data.principal.refseq108.hg38.txt", 
                     stringsAsFactors = F, header = F),
                read.delim("appris_data.principal.refseq109.hg38.txt", 
                           stringsAsFactors = F, header = F))
appris$V3 <- sapply(sapply(appris$V3, strsplit, "\\."), "[", 1)
introns <- introns[introns$isoform %in% appris$V3, ]
introns.g <- GRanges(introns$chr, IRanges(introns$start, introns$end), "*")


exons <- read.delim("hg38_exon_refFlat.txt", stringsAsFactors = F)

rare.expansion <- read.delim("tabulated.EHdn.outliers.overlapped.any.tsv", stringsAsFactors = F)
rare.expansion <- rare.expansion[rare.expansion$chr %in% paste0("chr", c(1:22)), ]
rare.expansion$size <- rare.expansion$end - rare.expansion$start
rare.expansion <- rare.expansion[rare.expansion$size < 10000, ]

detected.expansion <- read.delim("merged.any.expansions.autosomal.sex.chr.nobadsamples.tsv", stringsAsFactors = F)
detected.expansion <- detected.expansion[detected.expansion$chr %in% paste0("chr", c(1:22)), ]
detected.expansion$size <- detected.expansion$end - detected.expansion$start
detected.expansion <- detected.expansion[detected.expansion$size < 10000, ]

rare.expansion <- unique(rare.expansion[, c("chr", "start", "end")])
common.expansion <- unique(known.expansion[, c("seqnames", "start", "end")])
common.expansion <- common.expansion[common.expansion$end-common.expansion$start < 10000, ] #0.01%

detected.expansion <- unique(detected.expansion[, c("chr", "start", "end")])
names(common.expansion) <- names(rare.expansion)

rare.expansion.g <- GRanges(rare.expansion$chr, IRanges(rare.expansion$start, rare.expansion$end), "*")
common.expansion.g <- GRanges(common.expansion$chr, IRanges(common.expansion$start, common.expansion$end), "*")
detected.expansion.g <- GRanges(detected.expansion$chr, IRanges(detected.expansion$start, detected.expansion$end), "*")
common.expansion.g <- reduce(common.expansion.g)
common.expansion <- data.frame(common.expansion.g)
names(common.expansion) <- c("chr", "start", "end", "width", "strand")

getsplicedistance <- function(introns.g, rare.expansion.g, introns, rare.expansion, exons){
  introns.real.g <- GRanges(introns$chr, IRanges(introns$start, introns$end), "*")
  rare.olap <- data.frame(findOverlaps(rare.expansion.g, introns.real.g))
  rare.olap$isoform <- introns$isoform[rare.olap$subjectHits]
  dup <- rare.olap$queryHits[which(duplicated(rare.olap[, c("queryHits", "isoform")]))]
  
  
  exons.g <- GRanges(exons$chr, IRanges(exons$start, exons$end), "*")
  rare.exon <-  data.frame(findOverlaps(rare.expansion.g, exons.g))
  rare.exon$sizeOlap <- width(pintersect(rare.expansion.g[rare.exon$queryHits], 
                                         exons.g[rare.exon$subjectHits]))
  rare.exon$sizeExon <- width(rare.expansion.g[rare.exon$queryHits])
  rare.exon <- rare.exon[rare.exon$sizeOlap == rare.exon$sizeExon, ]
  
  dup <- union(dup, rare.exon$queryHits)
  
  rare.olap <- data.frame(findOverlaps(rare.expansion.g, introns.g))
  rare.olap <- rare.olap[!rare.olap$queryHits %in% rare.olap$queryHits[dup], ]
  rare.olap$strand <- introns$strand[rare.olap$subjectHits]
  rare.olap$chr <- rare.expansion$chr[rare.olap$queryHits]
  rare.olap$start.expansion <- rare.expansion$start[rare.olap$queryHits]
  rare.olap$end.expansion <- rare.expansion$end[rare.olap$queryHits]
  rare.olap$start.intron <- introns$start[rare.olap$subjectHits]
  rare.olap$end.intron <- introns$end[rare.olap$subjectHits]
  rare.olap$mid.point <- round((rare.olap$start.expansion + rare.olap$end.expansion) / 2)
  
  rare.olap$donor.distance <- ifelse(rare.olap$strand == "+", rare.olap$mid.point - rare.olap$start.intron, 
                                     rare.olap$end.intron - rare.olap$mid.point)
  rare.olap$acceptor.distance <- ifelse(rare.olap$strand == "+", rare.olap$end.intron - rare.olap$mid.point, 
                                        rare.olap$mid.point - rare.olap$start.intron)
  
  
  rare.olap$nearest.site <- ifelse(abs(rare.olap$donor.distance) < abs(rare.olap$acceptor.distance), "Donor", "Acceptor")
  rare.olap$distance <- ifelse(rare.olap$nearest.site == "Donor", rare.olap$donor.distance, rare.olap$acceptor.distance) %>% abs()
  
  rare.olap <- rare.olap[order(abs(rare.olap$distance)), ]
  rare.olap <- rare.olap[!duplicated(rare.olap$queryHits), ]
  
  rare.olap <- rare.olap[abs(rare.olap$distance) < 10000, ]

  rare.olap$distance[which(rare.olap$nearest.site == "Acceptor")] <-
    -rare.olap$distance[which(rare.olap$nearest.site == "Acceptor")]
  return(rare.olap)
}

rare.distance <- getsplicedistance(introns.g, rare.expansion.g, introns, rare.expansion, exons)
detected.distance <- getsplicedistance(introns.g, detected.expansion.g, introns, detected.expansion, exons)
common.distance <- getsplicedistance(introns.g, common.expansion.g, introns, common.expansion, exons)

rare.distance$rarity <- "rare expansions"
common.distance$rarity <- "known STRs"
detected.distance$rarity <- "all expansions"

wilcox.test(abs(rare.distance$distance), 
            abs(detected.distance$distance), alternative = "less")$p.value #0.0006
wilcox.test(abs(rare.distance$distance), 
            abs(common.distance$distance), alternative = "less")$p.value #1e-5

distance <- rbind(rare.distance, detected.distance, common.distance)

ggplot(distance, aes(x = distance, y = ..density.., color = rarity)) +
  geom_density(alpha = .15, adjust = 1/4) +
  geom_vline(xintercept = c(0), lty = 2) + xlab("Distance from splice junction") + theme_classic() +
  theme(panel.border = element_rect(fill = NA), 
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
  legend.position = "none") +  
  geom_hline(yintercept = 0, size = 1, color = "white") + 
  scale_x_continuous(breaks = c(-5000, 0, 5000)) +
  scale_color_manual(values = c("#4393C3", "#338833", "#D6604D"))

###### figure 2b
selected.columns <- c("all.aggr.glm.OR",  
                      "all.aggr.glm.upper", 
                      "all.aggr.glm.lower", 
                      "all.aggr.glm.pvalue"
                      )
global.genic <- read.delim("global.any.motif.genic.tsv", stringsAsFactors = F)[, selected.columns]
global.all <- read.delim(".global.any.motif.allexpansions.tsv", stringsAsFactors = F)[, selected.columns]
global.intergenic <- read.delim("functional.any.motif.tsv", stringsAsFactors = F)[8, selected.columns]
names(global.genic) <- c("OR", "upper", "lower", "pvalue")
names(global.all) <- names(global.genic)
names(global.intergenic) <- names(global.genic)

global.genic$type <- "Genic expansions"
global.all$type <- "All expansions"
global.intergenic$type <- "Intergenic expansions"

global <- rbind(rbind(global.genic[1,], global.intergenic), global.all[1, ])
global$minuslog10pvalue <- -log10(unlist(global$pvalue))
global$OR <- unlist(global$OR)
global$lower <- unlist(global$lower)
global$upper <- unlist(global$upper)
global$type <- factor(global$type, levels = c("All expansions", "Intergenic expansions", "Genic expansions"))
ggplot(global, aes(x = type, y = OR)) +
  geom_bar(stat = "identity", color = "black", fill = "#2166AC", width = .5) + 
  geom_hline(yintercept = 1, lty = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5), size = 0.5, width = 0.2) +
  theme_classic() + ylab("Odds ratio") + xlab("") + scale_y_continuous(breaks = c(0,1), limits = c(0, 1.8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(fill = NA)) +
annotate("text", x = c(3), y = c(1.78), label =
           paste0("p==", c("3.5~x~10^-03")), parse = T, fontface = "italic", cex = 4)

###### figure 2f
gs.name <- read.delim("geneset.name.txt", stringsAsFactors = F, header = F)
geneset <- read.delim("geneset.any.motif.tsv", stringsAsFactors = F)
geneset <- geneset[geneset$geneset != "sig.geneset", ]
geneset <- merge(geneset, gs.name, by.x = "geneset", by.y = "V2", all.x = T)

geneset <- geneset[order(geneset$all.aggr.glm.OR, decreasing = F), ]
geneset$V1 <- factor(geneset$V1, levels = geneset$V1)
geneset$minuslog10pvalue <- -log10(geneset$all.aggr.glm.pvalue)

ggplot(geneset, aes(x = minuslog10pvalue, y = all.aggr.glm.OR)) +
  geom_point(aes(color = aggr.all.fwer < 0.2), show.legend = F) +
  geom_errorbar(aes(ymin = all.aggr.glm.lower, 
                    ymax = all.aggr.glm.upper,
                    color = aggr.all.fwer < 0.2), size = .3, width = .1, show.legend = F) + 
  geom_hline(yintercept = 1, lty = 2) +
  annotate("line", x = c(1.2, 2.7,
                         geneset$minuslog10pvalue[geneset$aggr.all.fwer < 0.2]),
           y = c(0.4,   0.55,
                 geneset$all.aggr.glm.OR[geneset$aggr.all.fwer < 0.2]), 
           group = c("A", "A", "B", "B")) +
  annotate("text_repel", label = c("GO:Nervous system dev",
                                   "MPO: Cardiovascular & Muscle"), 
           x = c(1,  2.5), 
           y = c(0.3,   0.42), cex = 5) + coord_cartesian(ylim = c(0, 4.5)) +
  scale_color_manual(values = c("black", "#D6604D")) +
  theme_classic() + ylab("Odds ratio") + xlab("-log10(p-value)") + 
  theme(panel.border = element_rect(fill = NA),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), 
        axis.title.x = element_text(size = 14)) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4)) +
  scale_x_continuous(breaks = c(0, 3), limits = c(0, 3.5))

###### figure 3a
rare.expansion <- read.delim("Merged.Any.Expansion.with.outliers.EHdn.tsv", stringsAsFactors = F)[,
                                                c("chr.x", "repeatID", "samples", "gene_symbol", "insigset")]
rare.expansion <- rare.expansion[rare.expansion$chr.x %in% paste0("chr", 1:22), ]
sample.genic <- unique(unlist(sapply(rare.expansion$samples[!is.na(rare.expansion$gene_symbol)], strsplit, ";")))
sample.all <- unique(unlist(sapply(rare.expansion$samples, strsplit, ";")))
sample.gs <- unique(unlist(sapply(rare.expansion$samples[rare.expansion$insigset == T], strsplit, ";")))


all.samples <- read.delim("fam.data.kids.clean.tsv", stringsAsFactors = F)
all.samples <- all.samples[all.samples$Status == "AffectedKid", ]

all.f <- sum(all.samples$Sex == "female")
all.m <- sum(all.samples$Sex == "male")
all.rare.f <- sum(all.samples$Sex == "female" & all.samples$Sample.ID %in% sample.all)
all.rare.m <- sum(all.samples$Sex == "male" & all.samples$Sample.ID %in% sample.all)
all.genic.f <- sum(all.samples$Sex == "female" & all.samples$Sample.ID %in% sample.genic)
all.genic.m <- sum(all.samples$Sex == "male" & all.samples$Sample.ID %in% sample.genic)
all.gs.f <- sum(all.samples$Sex == "female" & all.samples$Sample.ID %in% sample.gs)
all.gs.m <- sum(all.samples$Sex == "male" & all.samples$Sample.ID %in% sample.gs)

test.rare <- fisher.test(data.frame("sex" = c(all.rare.f, all.rare.m), "all" = c(all.f - all.rare.f, all.m - all.rare.m)))
test.genic <- fisher.test(data.frame("sex" = c(all.genic.f, all.genic.m), "all" = c(all.f - all.genic.f, all.m - all.genic.m)))
test.gs <- fisher.test(data.frame("sex" = c(all.gs.f, all.gs.m), "all" = c(all.f - all.gs.f, all.m - all.gs.m)))

dt <- data.frame("Region" = c("All rare expansions", "Genic rare expansions", "Expansion in enriched gene-sets"),
                 "OR" = c(test.rare$estimate, test.genic$estimate, test.gs$estimate),
                 "OR.up" = c(test.rare$conf.int[2], test.genic$conf.int[2], test.gs$conf.int[2]),
                 "OR.low" = c(test.rare$conf.int[1], test.genic$conf.int[1], test.gs$conf.int[1]), 
                 "pvalue" = c(test.rare$p.value, test.genic$p.value, test.gs$p.value), stringsAsFactors = F)

dt$Region <- factor(dt$Region, levels = dt$Region)
ggplot(dt, aes(x = Region, y = OR)) + 
  geom_bar(stat = "identity", width = .5, color = "black", fill = "#2166AC") +
  geom_errorbar(aes(ymin = OR.low, ymax = OR.up), width = .3) +   geom_hline(yintercept = 1, lty = 2) + xlab("") +
  theme_classic() + ylab("Female/male odds ratio") + theme(panel.border = element_rect(fill = NA),
                                                           axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                                               axis.text.y = element_text(size = 14), 
                                               axis.title.y = element_text(size = 12)) +
  annotate("text", x = 3, y = 2.1, label = "p=0.11", fontface = "italic", size = 5) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(0, 2.5))

###### figure 3b
rare.expansion <- read.delim("tabulated.EHdn.outliers.overlapped.any.tsv", stringsAsFactors = F)
rare.expansion <- rare.expansion[rare.expansion$chr %in% paste0("chr", 1:22), ]
fam.data <- read.delim("fam.data.kids.clean.tsv", stringsAsFactors = F)
fam.data <- fam.data[fam.data$Status %in% c("AffectedKid", "UnaffectedKid"), ]
outlier <- read.delim("MSSNG.samples.with.excessive.no.detected.EHdn.both.sides.tsv", stringsAsFactors = F)
fam.data <- fam.data[!fam.data$Sample.ID %in% outlier$Sample, ]

sample.rare <- unique(rare.expansion$Sample.ID)
sample.rare <- sample.rare[sample.rare %in% fam.data$Sample.ID[fam.data$Dataset == "MSSNG"]]

sample.common <- fam.data$Sample.ID[!fam.data$Sample.ID %in% sample.rare & fam.data$Dataset == "MSSNG"]

vineland <- read.delim("../phenotypeData/VABSSCOREALL.csv", sep=",", stringsAsFactors = F)
iqall <- read.delim("../phenotypeData/IQSCOREALL.csv", sep=",", stringsAsFactors = F)
owls <- read.delim("../phenotypeData/LANGOWLSCOREALL.csv", sep=",", stringsAsFactors = F)
repetitve.revise.behavior <- read.delim("../phenotypeData/RBSRSCOREALL.csv", sep=",", stringsAsFactors = F)
repetitve.behavior <- read.delim("../phenotypeData/RBSSCOREALL.csv", sep=",", stringsAsFactors = F)

social.commu <- read.delim("../phenotypeData/SCQSCOREALL.csv", sep=",", stringsAsFactors = F)
social.resp <- read.delim("../phenotypeData/SRSSCOREALL.csv", sep=",", stringsAsFactors = F)

iqall$max <- NA
for(i in 1:nrow(iqall)){
  iqall$max[i] <- max(iqall[i, 2:5], na.rm = T)
}
iqall$max[which(iqall$max == 999 | iqall$max == 777 | iqall$max == -Inf | iqall$max > 200)] <- NA

vineland <- vineland[vineland$INDEXID %in% c(sample.rare, sample.common), ]
iqall <- iqall[iqall$INDEXID %in% c(sample.rare, sample.common), ]
owls <- owls[owls$INDEXID %in% c(sample.rare, sample.common), ]
repetitve.behavior <- repetitve.behavior[repetitve.behavior$INDEXID %in% c(sample.rare, sample.common), ]
social.commu <- social.commu[social.commu$INDEXID %in% c(sample.rare, sample.common), ]
social.resp <- social.resp[social.resp$INDEXID %in% c(sample.rare, sample.common), ]

vineland$label <- ifelse(vineland$INDEXID %in% sample.rare, "Rare expansion", "No rare expansion")
iqall$label <- ifelse(iqall$INDEXID %in% sample.rare, "Rare expansion", "No rare expansion")
owls$label <- ifelse(owls$INDEXID %in% sample.rare, "Rare expansion", "No rare expansion")
repetitve.behavior$label <- ifelse(repetitve.behavior$INDEXID %in% sample.rare, "Rare expansion", "No rare expansion")
social.commu$label <- ifelse(social.commu$INDEXID %in% sample.rare, "Rare expansion", "No rare expansion")
social.resp$label <- ifelse(social.resp$INDEXID %in% sample.rare, "Rare expansion", "No rare expansion")

samples <- unique(c(vineland$INDEXID, iqall$INDEXID, owls$INDEXID, repetitve.behavior$INDEXID,
                    social.commu$INDEXID, social.resp$INDEXID))
sum(sample.rare %in% samples) #485
sum(sample.common %in% samples) #1194

vineland$score <- "Vineland \nadaptive behavior\n(n=1113)"
iqall$score <- "IQ full scale\n(n=565)"
owls$score <- "OWLS test\n(n=588)"
repetitve.behavior$score <- "Repetitive behavior\n(n=155)"
social.commu$score <- "Social communication\n(n=701)"
social.resp$score <- "Social responsiveness\n(n=531)"

vineland <- na.omit(vineland[, c("max", "label", "score")])
iqall <- na.omit(iqall[, c("max", "label", "score")])
owls <- na.omit(owls[, c("max", "label", "score")])
repetitve.behavior <- na.omit(repetitve.behavior[, c("max", "label", "score")])
social.commu <- na.omit(social.commu[, c("max", "label", "score")])
social.resp <- na.omit(social.resp[, c("max", "label", "score")])

getPercentile <- function(values){
  return(sapply(lapply(values, ">", values), sum)/length(values))
}

vineland$perc <- getPercentile(vineland$max)
iqall$perc <- getPercentile(iqall$max)
owls$perc <- getPercentile(owls$max)
repetitve.behavior$perc <- getPercentile(repetitve.behavior$max)
social.commu$perc <- getPercentile(social.commu$max)
social.resp$perc <- getPercentile(social.resp$max)

dt.pheno <- rbind(owls, rbind(repetitve.behavior, rbind(social.commu, rbind(social.resp, rbind(vineland, iqall)))))
dt.pheno$score <- factor(dt.pheno$score, levels = c("IQ full scale\n(n=565)",
                                                    "Vineland \nadaptive behavior\n(n=1113)",
                                                    "OWLS test\n(n=588)",
                                                    "Repetitive behavior\n(n=155)",
                                                    "Social communication\n(n=701)",
                                                    "Social responsiveness\n(n=531)"))
pab <- wilcox.test(vineland$max[vineland$label == "Rare expansion"],
              vineland$max[vineland$label == "No rare expansion"], conf.int = T, alternative = "less")
pab$estimate/sqrt(nrow(vineland)) #-0.06
pab$conf.int[1]/sqrt(nrow(vineland)) #-Inf
pab$conf.int[2]/sqrt(nrow(vineland)) #-9.6e-07

pab.pvalue <- pab$p.value

piq <- wilcox.test(iqall$max[iqall$label == "Rare expansion"],
                   iqall$max[iqall$label == "No rare expansion"], conf.int = T, alternative = "less")
piq.pvalue <- piq$p.value
piq$estimate/sqrt(nrow(iqall)) #-0.29
piq$conf.int[1]/sqrt(nrow(iqall)) #-Inf
piq$conf.int[2]/sqrt(nrow(iqall)) #-0.13

powls <- wilcox.test(owls$max[owls$label == "Rare expansion"],
                   owls$max[owls$label == "No rare expansion"])
powls.pvalue <- powls$p.value

prep <- wilcox.test(repetitve.behavior$max[repetitve.behavior$label == "Rare expansion"],
                   repetitve.behavior$max[repetitve.behavior$label == "No rare expansion"], alternative = "less", conf.int = T)
prep.pvalue <- prep$p.value
prep$estimate/sqrt(nrow(repetitve.behavior)) #-0.64
prep$conf.int[1]/sqrt(nrow(repetitve.behavior)) #-Inf
prep$conf.int[2]/sqrt(nrow(repetitve.behavior)) #-0.24

pcom <- wilcox.test(social.commu$max[social.commu$label == "Rare expansion"],
                   social.commu$max[social.commu$label == "No rare expansion"])
pcom.pvalue <- pcom$p.value


pres <- wilcox.test(social.resp$max[social.resp$label == "Rare expansion"],
                   social.resp$max[social.resp$label == "No rare expansion"])
pres.pvalue <- pres$p.value

library(ggplot2)
library(fmsb)

ggplot(dt.pheno[dt.pheno$score %in% c(
  "IQ full scale\n(n=565)", "Vineland \nadaptive behavior\n(n=1113)" #, "Repetitive behavior\n(n=155)"
), ], aes(x = score, y = perc*100, fill = label)) + geom_violin(width = 0.4, color = NA, alpha = 0) + theme_classic() +
  ylab("Percentile") + xlab("") + 
  annotate("line", x = c(0.8, 1.2, 1.8, 2.2), y = c(105, 105, 105, 105),
           group = c("A", "A", "B", "B")) +
  annotate("text", x = c(1, 2), y = c(110, 110), label = c(
    sprintf("p=%s", round(piq.pvalue, digits = 3)),
    sprintf("p=%s", round(pab.pvalue, digits = 3))
  ), cex = 3, fontface="italic") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                                         axis.text.y = element_text(size = 12),
                                         axis.title.y = element_text(size = 12),
                      panel.border = element_rect(fill = NA),
                      legend.title = element_blank(),
                      legend.position = c(0.5, 0.9),
                      legend.direction = "horizontal",
                      legend.key.size = unit(0.3, "cm"),
                      legend.text = element_text(size = 9)) +
  geom_boxplot(position = position_dodge(width = 0.6), width = .4) +
  scale_y_continuous(breaks = c(0, 90), limits = c(0, 130)) +
  scale_fill_manual(values = c("#4393C3", "#D6604D"))

ggsave("../plots/pdf/figure3b.pdf", width = 4, height = 4)
write.table(dt.pheno[dt.pheno$score %in% c(
  "IQ full scale\n(n=565)", "Vineland \nadaptive behavior\n(n=1113)" #, "Repetitive behavior\n(n=155)"
), ], "../plots/pdf/figure3b.underly.tsv", sep="\t", row.names=F, quote=F, col.names=T)

table(dt.pheno$label[dt.pheno$score == "IQ full scale\n(n=565)"])
table(dt.pheno$label[dt.pheno$score == "Vineland \nadaptive behavior\n(n=1113)"])

