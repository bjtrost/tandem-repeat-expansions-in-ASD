## script to generate figures based on transmission test results

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(data.table)
library(GenomicRanges)
library(GenomicScores)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(magrittr)

#######################
######### figure 3c
rare.function <- read.delim("ssc.affected.99perc.element.remained_repeats_as_background.tsv", stringsAsFactors = F)
rare.function <- rare.function[rare.function$gs.query %in% c("upstream", "downstream", "intron", "UTR5", "UTR3", "exon", "splicing"), ]
rare.function <- rare.function[order(rare.function$gs.query), ]
rare.function$geneset <- factor(c("downstream", "exon", "intron", "splicing", "upstream", "3'UTR", "5'UTR"), 
                                levels = c("upstream", "5'UTR", "exon", "splicing", "intron", "3'UTR", "downstream"))

p1 <- ggplot(rare.function, aes(x = geneset, y = or, fill = pvalue < 0.05)) + 
  geom_bar(stat="identity", show.legend = F, color = "black", width = .5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), show.legend = F,
                position = position_dodge(width = .5), size = .4, width = .2) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_classic() + ylab("Odds Ratio") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(fill = NA)) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_fill_manual(values = c("#4393C3", "#D6604D")) + 
  scale_y_continuous(breaks = c(0, 1, 3)) +
  annotate("text", x = c(0.7, 1.7, 4.9, 6.9), y = c(3, 3.4, 2.2, 2), angle = 90,
           label = paste0("p==", c("5.4~x~10^-3",
                                   "5.4~x~10^-8",
                                   "2.7~x~10^-2",
                                   "4.9~x~10^-2")), parse = T,
           fontface = "italic", cex = 4)
###mssng
rare.function <- read.delim("mssng.affected.99perc.element.remained_repeats_as_background.tsv", stringsAsFactors = F)
rare.function <- rare.function[rare.function$gs.query %in% c("upstream", "downstream", "intron", "UTR5", "UTR3", "exon", "splicing"), ]
rare.function <- rare.function[order(rare.function$gs.query), ]
rare.function$geneset <- factor(c("downstream", "exon", "intron", "splicing", "upstream", "3'UTR", "5'UTR"), 
                                levels = c("upstream", "5'UTR", "exon", "splicing", "intron", "3'UTR", "downstream"))

p2 <- ggplot(rare.function, aes(x = geneset, y = or, fill = pvalue < 0.05)) + 
  geom_bar(stat="identity", show.legend = F, color = "black", width = .5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), show.legend = F,
                position = position_dodge(width = .5), size = .4, width = .2) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_classic() + ylab("Odds Ratio") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(fill = NA)) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_fill_manual(values = c("#4393C3", "#D6604D")) + 
  scale_y_continuous(breaks = c(0, 1, 3)) +
  ###mssng
  annotate("text", x = c(1.7, 2.7, 5.9), y = c(2.7, 3, 2), angle = 90,
         label = paste0("p==", c("0.01",
                                 "5.4~x~10^-4",
                                 "2.3~x~10^-2")), parse = T,
         fontface = "italic", cex = 4)

###eur
rare.function <- read.delim("eur.affected.99perc.element.remained_repeats_as_background.tsv", stringsAsFactors = F)
rare.function <- rare.function[rare.function$gs.query %in% c("upstream", "downstream", "intron", "UTR5", "UTR3", "exon", "splicing"), ]
rare.function <- rare.function[order(rare.function$gs.query), ]
rare.function$geneset <- factor(c("downstream", "exon", "intron", "splicing", "upstream", "3'UTR", "5'UTR"), 
                                levels = c("upstream", "5'UTR", "exon", "splicing", "intron", "3'UTR", "downstream"))

p3 <- ggplot(rare.function, aes(x = geneset, y = or, fill = pvalue < 0.05)) + 
  geom_bar(stat="identity", show.legend = F, color = "black", width = .5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), show.legend = F,
                position = position_dodge(width = .5), size = .4, width = .2) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_classic() + ylab("Odds Ratio") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(fill = NA)) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_fill_manual(values = c("#4393C3", "#D6604D")) + 
  scale_y_continuous(breaks = c(0, 1, 3)) +
  ###mssng
  annotate("text", x = c(0.7, 1.7, 2.9, 5.9), y = c(3, 3.1, 2.8, 1.8), angle = 90,
           label = paste0("p==", c("2.3~x~10^-4",
                                     "1.1~x~10^-10",
                                   "7.1~x~10^-7",
                                   "9.1~x~10^3")), parse = T,
           fontface = "italic", cex = 4)

plot_grid(p1, p2, p3, labels = c("a", "b", "c"), nrow = 1)


###### figure 3d
selected.columns <- c("all.aggr.glm.OR",  
                      "all.aggr.glm.upper", 
                      "all.aggr.glm.lower", 
                      "all.aggr.glm.pvalue"
                      )
global <- read.delim("ssc.affected.99perc.element.remained_repeats_as_background.tsv", stringsAsFactors = F)
global <- global[global$gs.query %in% c("all", "genic"), ]
rare.function <- read.delim("ssc.affected.99perc.element.remained_repeats_as_background.tsv", stringsAsFactors = F)
global <- rbind(global, rare.function)
global$gs.query <- factor(global$gs.query, levels = c("intergenic", "genic"))
p1 <- ggplot(global[global$gs.query %in% c("genic", "intergenic"), ], 
             aes(x = gs.query, y = or)) + 
  geom_bar(stat = "identity", color = "black", fill = "#2166AC", width = .5) + 
  geom_hline(yintercept = 1, lty = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5), size = 0.5, width = 0.2) +
  # geom_point(aes(y = minuslog10pvalue), position = position_dodge(width = 0.9), shape = 8, show.legend = F) +
  theme_classic() + ylab("Odds ratio") + xlab("SSC") + scale_y_continuous(breaks = c(0,1), limits = c(0, 1.4)) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(fill = NA)) +
annotate("text", x = c(2), y = c(1.35), label =
           paste0("p==", c("3.9~x~10^-5")), parse = T, fontface = "italic", cex = 4)

global <- read.delim("mssng.affected.99perc.element.remained_repeats_as_background.tsv", stringsAsFactors = F)
rare.function <- read.delim("mssng.affected.99perc.element.remained_repeats_as_background.tsv", stringsAsFactors = F)
global <- rbind(global, rare.function)
global$gs.query <- factor(global$gs.query, levels = c("intergenic", "genic"))

p2 <- ggplot(global[global$gs.query %in% c("genic", "intergenic"), ], aes(x = gs.query, y = or)) +
  geom_bar(stat = "identity", color = "black", fill = "#2166AC", width = .5) + 
  geom_hline(yintercept = 1, lty = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5), size = 0.5, width = 0.2) +
  theme_classic() + ylab("Odds ratio") + xlab("MSSNG") + scale_y_continuous(breaks = c(0,1), limits = c(0, 1.4)) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(fill = NA))
  # annotate("text", x = c(2), y = c(1.35), label =
  #            paste0("p==", c("4~x~10^-03")), parse = T, fontface = "italic", cex = 4)

global <- read.delim("eur.affected.99perc.element.remained_repeats_as_background.tsv", stringsAsFactors = F)
rare.function <- read.delim("eur.affected.99perc.element.remained_repeats_as_background.tsv", stringsAsFactors = F)
global <- rbind(global, rare.function)
global$gs.query <- factor(global$gs.query, levels = c("intergenic", "genic"))

p3 <- ggplot(global[global$gs.query %in% c("genic", "intergenic"), ], aes(x = gs.query, y = or)) +
  geom_bar(stat = "identity", color = "black", fill = "#2166AC", width = .5) + 
  geom_hline(yintercept = 1, lty = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5), size = 0.5, width = 0.2) +
  # geom_point(aes(y = minuslog10pvalue), position = position_dodge(width = 0.9), shape = 8, show.legend = F) +
  theme_classic() + ylab("Odds ratio") + xlab("Combined") + scale_y_continuous(breaks = c(0,1), limits = c(0, 1.4)) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(fill = NA))+
annotate("text", x = c(2), y = c(1.30), label =
           paste0("p==", c("1.6~x~10^-7")), parse = T, fontface = "italic", cex = 4)

plot_grid(p1, p2, p3, labels = c("a", "b", "c"), nrow = 1)

###### figure 3e
plot_figure11 <- function(){
  gs.name <- read.delim("geneset.name.txt", stringsAsFactors = F, header = F)
  geneset <- read.delim("ssc.affected.99perc.geneset.remained_repeats_as_background.tsv", stringsAsFactors = F)
  names(geneset)[1] <- c("geneset")
  
  geneset <- merge(geneset, gs.name, by.x = "geneset", by.y = "V2", all.x = T)
  
  geneset <- geneset[order(geneset$or, decreasing = F), ]
  geneset$V1[!geneset$geneset %in% c("PhMm_Aggr_CardvascMuscle_all", "Neurof_GoNervSysDev")] <- ""
  geneset$minuslog10pvalue <- -log10(geneset$pvalue)
  
  p1 <- ggplot(geneset, aes(x = minuslog10pvalue, y = or)) +
    geom_point(aes(color = FWER < 0.25), show.legend = F) +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper,
                      color = FWER < 0.25), size = .3, width = .1, show.legend = F) + 
    geom_hline(yintercept = 1, lty = 2) +
    annotate("line", x = c(1.2, 2.7,
                           geneset$minuslog10pvalue[geneset$geneset %in% c("PhMm_Aggr_CardvascMuscle_all", "Neurof_GoNervSysDev")]),
             y = c(2.3,   2.82,
                   geneset$or[geneset$geneset %in% c("PhMm_Aggr_CardvascMuscle_all", "Neurof_GoNervSysDev")]), 
             group = c("A", "A", "B", "B")) +
    annotate("text_repel", label = c(
      "MPO: Cardiovascular & Muscle", "GO:Nervous system dev"), 
             x = c(1,  2.5), 
             y = c(2.7,   2.82), cex = 5) + coord_cartesian(ylim = c(0, 3)) +
    scale_color_manual(values = c("black", "#D6604D")) +
    theme_classic() + ylab("Odds ratio") + xlab("-log10(p-value)") + 
    theme(panel.border = element_rect(fill = NA),
          axis.text.y = element_text(size = 14), 
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14), 
          axis.title.x = element_text(size = 14)) +
    scale_y_continuous(breaks = c(0, 1, 2, 3)) +
    scale_x_continuous(breaks = c(0, 3), limits = c(0, 4.5))
  
  ###mssng
  geneset <- read.delim("mssng.affected.99perc.geneset.remained_repeats_as_background.tsv", stringsAsFactors = F)
  names(geneset)[1] <- c("geneset")
  
  geneset <- merge(geneset, gs.name, by.x = "geneset", by.y = "V2", all.x = T)
  
  geneset <- geneset[order(geneset$or, decreasing = F), ]
  geneset$V1[!geneset$geneset %in% c("PhMm_Aggr_CardvascMuscle_all", "Neurof_GoNervSysDev")] <- ""
  geneset$minuslog10pvalue <- -log10(geneset$pvalue)
  
  p2 <- ggplot(geneset, aes(x = minuslog10pvalue, y = or)) +
    geom_point(aes(color = FWER < 0.25), show.legend = F) +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper,
                      color = FWER < 0.25), size = .3, width = .1, show.legend = F) + 
    geom_hline(yintercept = 1, lty = 2) +
    annotate("line", x = c(1.2, 3.4,
                           geneset$minuslog10pvalue[geneset$geneset %in% c("PhMm_Aggr_CardvascMuscle_all", "Neurof_GoNervSysDev")]),
             y = c(2.3,   2.85,
                   geneset$or[geneset$geneset %in% c("PhMm_Aggr_CardvascMuscle_all", "Neurof_GoNervSysDev")]), 
             group = c("A", "A", "B", "B")) +
    annotate("text_repel", label = c("GO:Nervous system dev", 
      "MPO: Cardiovascular & Muscle"), 
      x = c(1,  2.5), 
      y = c(2.7,   2.82), cex = 5) + coord_cartesian(ylim = c(0, 3)) +
    scale_color_manual(values = c("black", "#D6604D")) +
    theme_classic() + ylab("Odds ratio") + xlab("-log10(p-value)") + 
    theme(panel.border = element_rect(fill = NA),
          axis.text.y = element_text(size = 14), 
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14), 
          axis.title.x = element_text(size = 14)) +
    scale_y_continuous(breaks = c(0, 1, 2, 3)) +
    scale_x_continuous(breaks = c(0, 3), limits = c(0, 4))
  
  ###combine
  geneset <- read.delim("eur.affected.99perc.geneset.remained_repeats_as_background.tsv", stringsAsFactors = F)
  names(geneset)[1] <- c("geneset")
  
  geneset <- merge(geneset, gs.name, by.x = "geneset", by.y = "V2", all.x = T)
  
  geneset <- geneset[order(geneset$or, decreasing = F), ]
  geneset$V1[!geneset$geneset %in% c("PhMm_Aggr_CardvascMuscle_all", "Neurof_GoNervSysDev")] <- ""
  geneset$minuslog10pvalue <- -log10(geneset$pvalue)
  
  p3 <- ggplot(geneset, aes(x = minuslog10pvalue, y = or)) +
    geom_point(aes(color = FWER < 0.25), show.legend = F) +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper,
                      color = FWER < 0.25), size = .3, width = .1, show.legend = F) + 
    geom_hline(yintercept = 1, lty = 2) +
    annotate("line", x = c(2.7, 5.2,
                           geneset$minuslog10pvalue[geneset$geneset %in% c("PhMm_Aggr_CardvascMuscle_all", "Neurof_GoNervSysDev")]),
             y = c(2,   2.6,
                   geneset$or[geneset$geneset %in% c("PhMm_Aggr_CardvascMuscle_all", "Neurof_GoNervSysDev")]), 
             group = c("B", "A", "B", "A")) +
    annotate("text_repel", label = c(
      "MPO: Cardiovascular & Muscle", "GO:Nervous system dev"), 
      x = c(1,  4), 
      y = c(2.3,   2.82), cex = 5) + coord_cartesian(ylim = c(0, 3)) +
    scale_color_manual(values = c("black", "#D6604D")) +
    theme_classic() + ylab("Odds ratio") + xlab("-log10(p-value)") + 
    theme(panel.border = element_rect(fill = NA),
          axis.text.y = element_text(size = 14), 
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14), 
          axis.title.x = element_text(size = 14)) +
    scale_y_continuous(breaks = c(0, 1, 2, 3)) +
    scale_x_continuous(breaks = c(0, 3), limits = c(0, 7))
  plot_grid(p1, p2, p3, labels = c("a", "b", "c"), nrow = 3)
}

plot_figure11()