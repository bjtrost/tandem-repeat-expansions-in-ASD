# This script is to obtain families with transmitted allele in both affected and unaffected kids.
# It is to show that the expansion is already large in the parent and that it is further expanded in affected kid.
# To generate figure 2a in the paper.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list=ls())
library(data.table)
library(ggplot2)
library(ggExtra)
library(ggforce)
load("gs.asd.2019.RData")

annotation <- fread("merged.expansion.coordinates_forAnno.txt.annovar.out_rev38.1mFmt.tsv", stringsAsFactors = F)
annotation <- data.frame(annotation)
annotation$key <- paste(annotation$chr, annotation$start, annotation$end, sep="#")
annotation$typeseq_priority <- factor(annotation$typeseq_priority, 
                                      levels = c("exonic", "splicing", "exonic;splicing", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_exonic;ncRNA_splicing", 
                                                 "UTR5", "UTR3", "intronic", "ncRNA_intronic", "upstream", "downstream", "intergenic"))
annotation <- annotation[order(annotation$typeseq_priority), ]
annotation <- annotation[!duplicated(annotation$key), ]

ehdn.m <- read.delim("merged.any.expansions.autosomal.sex.chr.nobadsamples.tsv", stringsAsFactors = F) #chr, start, end, motif, v1, v2, genotype

annotation <- merge(annotation, ehdn.m, by.x = "key", by.y = "repeatID", all.x = T)
ehdn <- fread("output_regions.min2.1000G+SSC+MSSNG.txt") #chr, start, end, motif, v1, v2, genotype
ehdn <- as.data.frame(ehdn)
ehdn$repeatID <- paste(ehdn$V1, ehdn$V2, ehdn$V3, ehdn$V4, sep="#")

############

# annotation.tmp <- annotation #[!is.na(annotation$entrez_id) & !is.na(annotation$allRepeat), ]
# annotation.tmp <- annotation.tmp$allRepeat[annotation.tmp$entrez_id %in% unlist(gsASD[c("Neurof_GoNeuronBody")])]
# annotation.tmp <- strsplit(paste(annotation.tmp, collapse = ";"), ";")[[1]]


ssc.errors1 <- read.delim('ssc_ME.problems.txt',
                          header = F, stringsAsFactors = F)
ssc.errors2 <- read.delim('ssc_gender.problems.txt',
                          header = F, sep=" ", stringsAsFactors = F)
ssc.error <- c(ssc.errors1$V2[-27], ssc.errors2$V7, ssc.errors2$V11)
ssc.error <- unique(ssc.error[!ssc.error %in% c("", "0")])
all.error <- c(ssc.error, "2-1547-003", "3-0051-000", "3-0008-000", "7-0064-003", "7-0075-003", "3-0395-000")

brett.samples <- read.delim("manifest.1000G+SSC+MSSNG.txt", stringsAsFactors = F, header = F)
brett.samples <- brett.samples[!brett.samples$V1 %in% all.error, ]

tmp.outlier <- read.delim("MSSNG.samples.with.excessive.no.detected.EHdn.both.sides.tsv", stringsAsFactors = F)

ehdn.result.m <- read.delim("outlier.EHdn.MSSNG.SSC.1000G.Feb2020.no3SDoutlier.tsv", stringsAsFactors = F)
ehdn.result.m$repeatID <- paste(ehdn.result.m$chr, ehdn.result.m$start, ehdn.result.m$end, ehdn.result.m$motif, sep="#")
ehdn.result.m <- ehdn.result.m[ehdn.result.m$repeatID %in% ehdn$repeatID, ]

fam.data.tmp <- read.delim("Data_Merged_MSSNG.SSC.1KG.OCt162019.tsv", stringsAsFactors = F)
fam.data.tmp <- fam.data.tmp[fam.data.tmp$Sample.ID %in% brett.samples$V1, ]
fam.data.tmp <- merge(fam.data.tmp, readdepth, by.x = "Sample.ID", by.y = "Sample", all.x = T)

fam.data.tmp <- fam.data.tmp[!fam.data.tmp$Sample.ID %in% tmp.outlier$Sample, ]
fam.1000g <- fam.data.tmp[fam.data.tmp$Dataset == "1000G", ]
fam.data <- fam.data.tmp[which(!is.na(fam.data.tmp$Sample.ID)), ]
fam.data <- fam.data[fam.data$Status != "", ]

### limit to one affected kid per family
fam.data <- fam.data[which((!duplicated(paste(fam.data$Family.ID, fam.data$Status)) & fam.data$Status == "AffectedKid") |
                             fam.data$Status != "AffectedKid"), ]
fam.data <- fam.data[fam.data$DNASOURCE == "Other cell", ]
fam.data <- fam.data[fam.data$Library.type == "PCR-free" & fam.data$Platform == "Illumina HiSeqX", ]
fam.data <- fam.data[fam.data$Family.ID %in% names(which(table(fam.data$Family.ID) == 4)), ]
fam.data <- fam.data[fam.data$Dataset == "SSC" & fam.data$ETH_TAG == "EUR", ]

proband.fam <- unique(fam.data$Family.ID[fam.data$Status == "AffectedKid"])
both.parent.fam <- names(which(table(fam.data$Family.ID[fam.data$Status %in% c("UnaffectedParent", "AffectedParent")]) > 0))
sib.fam <- unique(fam.data$Family.ID[fam.data$Status == "UnaffectedKid"])

fam.data <- fam.data[fam.data$Family.ID %in% intersect(sib.fam, intersect(proband.fam, both.parent.fam)), ]

ehdn.result.m <- ehdn.result.m[ehdn.result.m$chr %in% paste0("chr", 1:22), ]

parent.geno <- c()
proband.geno <- c()
sib.geno <- c()
parent.id <- c()
proband <- c()

annotation$repeatID <- paste(annotation$key, annotation$motif, sep="#")
ehdn.result.m <- ehdn.result.m[ehdn.result.m$repeatID %in% annotation$repeatID[
  # annotation$typeseq_priority %in% c("exonic", "exonic;splicing", "splicing", "upstream", 
  #                                    "UTR3", "UTR5", "downstream", "intronic")
  annotation$entrez_id %in% c(gsASD$PhMm_Aggr_CardvascMuscle_all, gsASD$Neurof_GoNervSysDev)
], ]
genes <- c()
for(i in 1:nrow(ehdn.result.m)){
  tmp <- ehdn.result.m[i, ]
  ehdn.tmp <- ehdn[ehdn$repeatID == tmp$repeatID, ]
  ehdn.tmp.geno <- unique(strsplit(paste(ehdn.tmp$V7, collapse = ","), ",")[[1]])
  ehdn.tmp.geno <- data.frame(do.call(rbind, sapply(ehdn.tmp.geno, strsplit, ":")), row.names = NULL, stringsAsFactors = F)
  names(ehdn.tmp.geno) <- c("Sample", "EHdn")
  ehdn.tmp.geno <- ehdn.tmp.geno[ehdn.tmp.geno$Sample %in% fam.data$Sample.ID, ]
  if(nrow(ehdn.tmp.geno) > 0){
          no.geno <- nrow(fam.data) - nrow(ehdn.tmp.geno)
          ehdn.tmp.geno$EHdn <- as.numeric(as.character(ehdn.tmp.geno$EHdn))  
          ehdn.tmp.geno$EHdn <- (sapply(lapply(ehdn.tmp.geno$EHdn, ">", ehdn.tmp.geno$EHdn), sum) + no.geno)/(nrow(fam.data) - 1)
            
          samples <- fam.data[grep(gsub(";", "|", tmp$outliers), fam.data$Sample.ID), ]
          samples <- samples[samples$Status %in% c("AffectedKid"), ]
          
          fam.check <- c()
          if(nrow(samples) > 0)
          for(j in 1:nrow(samples)){
            if(fam.data$Family.ID[fam.data$Sample.ID == samples$Sample.ID[j]] %in% fam.check){
              message("skip")
            }else{
            # proband.tmp <- ehdn.tmp.geno[which(ehdn.tmp.geno$Sample == samples$Sample.ID[j]), ]\
                proband.tmp <- ehdn.tmp.geno[which(ehdn.tmp.geno$Sample %in% 
                                                            fam.data$Sample.ID[fam.data$Status == "AffectedKid" &
                                                                                 fam.data$Family.ID == samples$Family.ID]), ]
                parent <- ehdn.tmp.geno[which(ehdn.tmp.geno$Sample %in% c(samples$Mother.ID[j], samples$Father.ID[j])), ]
                parent <- parent[order(parent$EHdn, decreasing = T), ][1, ]
                sib <- ehdn.tmp.geno[which(ehdn.tmp.geno$Sample %in% 
                                             fam.data$Sample.ID[fam.data$Status == "UnaffectedKid" &
                                                                  fam.data$Family.ID == samples$Family.ID[j]]), ]
                if(nrow(proband.tmp) > 0 & nrow(parent) > 0 & nrow(sib) > 0){
                  if(sum(annotation$repeatID == tmp$repeatID) > 0)
                    if(!is.na(annotation$gene_symbol[annotation$repeatID == tmp$repeatID]) &
                       annotation$gene_symbol[annotation$repeatID == tmp$repeatID]== "DMPK")
                      message(i)
                  proband.geno <- c(proband.geno, max(proband.tmp$EHdn, na.rm = T))
                  proband <- c(proband, samples$Sample.ID[j])
                  if(nrow(parent) == 0){
                    #parent.geno <- c(parent.geno, 0.5)
                  }else{
                    parent.geno <- c(parent.geno, max(parent$EHdn, na.rm = T))
                    parent.id <- c(parent.id, parent$Sample)
                  }
                  
                  if(nrow(sib) == 0){
                    sib.geno <- c(sib.geno, 0.5)
                  }else{
                    sib.geno <- c(sib.geno, max(sib$EHdn))
                  }
                }
                
                fam.check <- c(fam.check, fam.data$Family.ID[fam.data$Sample.ID == samples$Sample.ID[j]])
            }
          }
  }
}

dt <- data.frame("Sample.ID" = proband, proband.geno, parent.geno, sib.geno, genes,  stringsAsFactors = F)

table(fam.data$Sex[fam.data$Sample.ID %in% parent.id])

dt <- merge(dt, fam.data, by = "Sample.ID")
dt$EHdn.diff <- dt$proband.geno - dt$parent.geno

dt.proband.ssc <- dt[dt$Dataset == "SSC", c("Sample.ID", "proband.geno")]
dt.parent.ssc <- dt[dt$Dataset == "SSC", c("Sample.ID", "parent.geno")]
dt.sib.ssc <- dt[dt$Dataset == "SSC", c("Sample.ID", "sib.geno")]

dt.proband.ssc$Relation <- "Proband"
dt.proband.ssc$Inclusion <- "SSC"
dt.parent.ssc$Relation <- "Parent"
dt.parent.ssc$Inclusion <- "SSC"
dt.sib.ssc$Relation <- "Sibs"
dt.sib.ssc$Inclusion <- "SSC"

dt.parent <- dt.parent.ssc
dt.proband <- dt.proband.ssc
dt.sib <- dt.sib.ssc

names(dt.parent)[2] <- "EHdn"
names(dt.proband)[2] <- "EHdn"
names(dt.sib)[2] <- "EHdn"

dt.plot <- rbind(dt.parent, dt.proband, dt.sib)

dt.plot <- merge(dt.plot, fam.data[, c("Sample.ID", "Sex")], by = "Sample.ID", all.x = T)

table(dt.plot$Sex[dt.plot$Relation == "Parent"])
b <- binom.test(184, 184+32)

dt.plot.ssc <- dt.plot[dt.plot$Inclusion == "SSC", ]
dt.plot.ssc$Relation <- gsub("Sibs", "Unaffected sibling", dt.plot.ssc$Relation)
names(dt.plot.ssc)[2] <- c("EHdn percentile")
ggplot(dt.plot.ssc, aes(x = Relation, y = `EHdn percentile` * 100, fill = Relation)) +
  geom_boxplot(outlier.shape = 1, width = 0.5) + ylab ("Percentile") + theme_classic() + xlab("") +
  theme(panel.border = element_rect(fill = NA), axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        legend.title = element_blank(), legend.direction = "horizontal", legend.position = "top",
        legend.key.size = unit(0.7, "cm"),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_y_continuous(breaks = c(50, 100)) +
  facet_zoom(ylim = c(99, 100), zoom.size = 1) +
  scale_fill_manual(values = c("#4393C3", "#D6604D", "#24AA34"))
