# This script is to test the association between expansions and other rare loss-of-function variants.
# It is used to generate supplementary figure 3 in the paper.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(cowplot)

annotation <- fread("merged.expansion.coordinates_forAnno.txt.annovar.out_rev38.1mFmt.tsv", stringsAsFactors = F)
annotation <- data.frame(annotation)
annotation$key <- paste(annotation$chr, annotation$start, annotation$end, sep="#")
annotation$typeseq_priority <- factor(annotation$typeseq_priority, 
                                      levels = c("exonic", "splicing", "exonic;splicing", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_exonic;ncRNA_splicing", 
                                                 "UTR5", "UTR3", "intronic", "ncRNA_intronic", "upstream", "downstream", "intergenic"))
annotation <- annotation[order(annotation$typeseq_priority), ]
annotation <- annotation[!duplicated(annotation$key), ]
annotation <- annotation[annotation$typeseq_priority != "intergenic", ]

load("gs.asd.2019.RData")
annotation <- annotation[annotation$entrez_id %in% c(gsASD$Neurof_GoNervSysDev, gsASD$PhMm_Aggr_CardvascMuscle_all), ]

del.out <- data.frame()
lof.out <- data.frame()
for(dataset in c("SSC", "MSSNG")){
  #### denovo genic deletions
  fam.1000g <- read.delim("Data_Merged_MSSNG.SSC.1KG.OCt162019.tsv", stringsAsFactors = F)
  aff.fam <- fam.1000g$Family.ID[fam.1000g$Status == "AffectedKid"]
  mpx.fam <- names(which(table(aff.fam) > 1))
  fam.1000g <- fam.1000g[fam.1000g$Dataset == "1000G", ]
  
  fam.data <- read.delim("fam.data.kids.clean.tsv", stringsAsFactors = F)
  fam.data <- fam.data[fam.data$Status == "AffectedKid" & fam.data$Dataset == dataset, ]
  fam.data$family <- "SPX"
  fam.data$family[fam.data$Family.ID %in% mpx.fam] <- "MPX"
  fam.data <- fam.data[fam.data$ETH_TAG == "EUR", ]

  ehdn.rare <- read.delim("merged.any.outlier.expansions.autosomal.sex.chr.nobadsamples.tsv", stringsAsFactors = F)
  ehdn.rare <- ehdn.rare[ehdn.rare$repeatID %in% annotation$key, ]
  ehdn.count <- sapply(ehdn.rare$outliers, strsplit, ";")
  ehdn.rare <- ehdn.rare[sapply(sapply(ehdn.count, "%in%", fam.1000g$Sample.ID), sum) < 2.5, ]
  ehdn.count <- paste(ehdn.rare$outliers, collapse = ";")
  ehdn.count <- strsplit(ehdn.count, ";")[[1]]
  ehdn.count <- data.frame(table(ehdn.count))
  names(ehdn.count)[2] <- "RareExpansions"
  
  denovo.cnv.samples <- read.delim(sprintf("QC.%s.txt", dataset), stringsAsFactors = F, header = F)
  denovo.cnv.samples <- denovo.cnv.samples[denovo.cnv.samples$V2 == "ok", ]
  fam.data <- fam.data[fam.data$Sample.ID %in% denovo.cnv.samples$V1, ]
  
  denovo.cnv <- read.delim(sprintf("%s_ERDS+_CNVs.variants_HQR.samples_passQC.p_denovo.tsv", dataset), stringsAsFactors = F)
  denovo.cnv <- denovo.cnv[denovo.cnv$SVTYPE == "DEL" & denovo.cnv$gene_symbol != "", ]
  
  fam.data$denovo_del <- F
  fam.data$denovo_del[fam.data$Sample.ID %in% denovo.cnv$sample] <- T
  
  fam.data <- merge(fam.data, ehdn.count, by.x = "Sample.ID", by.y = "ehdn.count", all.x = T)
  fam.data$RareExpansions[is.na(fam.data$RareExpansions)] <- 0
  
  denovo.expansion <- sum(fam.data$denovo_del & fam.data$RareExpansions > 0)
  denovo.noexpansion <- sum(fam.data$denovo_del & fam.data$RareExpansions == 0)
  nodenovo.expansion <- sum(!fam.data$denovo_del & fam.data$RareExpansions > 0)
  nodenovo.noexpansion <- sum(!fam.data$denovo_del & fam.data$RareExpansions == 0)
    
  del.test <- fisher.test(data.frame("X1" = c(denovo.expansion, nodenovo.expansion), "X2" = c(denovo.noexpansion, nodenovo.noexpansion)))
  del.out <- rbind(del.out, 
                   data.frame(
                     "Dataset" = dataset,
                     "OR" = del.test$estimate,
                     "lower" = del.test$conf.int[1],
                     "upper" = del.test$conf.int[2],
                     "p" = del.test$p.value
                   ))

  ssc.map <- read.delim("nygc_sfari_id_sample_map.data", stringsAsFactors = F, sep=",")
  all.ssc <- readLines("ssc.ids")
  all.mssng <- readLines("mssng_ilmn.ids")
  
  ssc.map <- ssc.map[ssc.map$SFARI.ID %in% all.ssc, ]
  lof.ssc <- read.delim("aat6576_Table-S2_denovo_mutations_ssc.hg38_withAnnotations.lof.tsv", stringsAsFactors = F)
  lof.ssc <- merge(lof.ssc, ssc.map, by.x = "SampleID", by.y = "SFARI.ID", all.x = T)
  lof.mssng <- read.delim("mssng_db6_ilmn_filtered_dng_20190929.lof.txt", stringsAsFactors = F)
  
  if(dataset == "MSSNG"){
    lof <- lof.mssng$X.sample_id # 
    all <- all.mssng # 
  }else{
    lof <- lof.ssc$Sample.ID #lof.mssng$X.sample_id # 
    all <- ssc.map$Sample.ID #all.mssng # 
  }

  
  fam.1000g <- read.delim("Data_Merged_MSSNG.SSC.1KG.OCt162019.tsv", stringsAsFactors = F)
  aff.fam <- fam.1000g$Family.ID[fam.1000g$Status == "AffectedKid"]
  mpx.fam <- names(which(table(aff.fam) > 1))
  
  fam.1000g <- fam.1000g[fam.1000g$Dataset == "1000G", ]
  
  fam.data <- read.delim("fam.data.kids.clean.tsv", stringsAsFactors = F)
  fam.data <- fam.data[fam.data$Status == "AffectedKid" & fam.data$Dataset == dataset, ]
  fam.data$family <- "SPX"
  fam.data$family[fam.data$Family.ID %in% mpx.fam] <- "MPX"
  #fam.data <- fam.data[fam.data$ETH_TAG == "EUR", ]
  
  # fam.data <- fam.data[fam.data$family == "MPX", ]
  
  ehdn.rare <- read.delim("merged.any.outlier.expansions.autosomal.sex.chr.nobadsamples.tsv", stringsAsFactors = F)
  ehdn.rare <- ehdn.rare[ehdn.rare$repeatID %in% annotation$key, ]
  ehdn.count <- sapply(ehdn.rare$outliers, strsplit, ";")
  ehdn.rare <- ehdn.rare[sapply(sapply(ehdn.count, "%in%", fam.1000g$Sample.ID), sum) < 2.5, ]
  ehdn.count <- paste(ehdn.rare$outliers, collapse = ";")
  ehdn.count <- strsplit(ehdn.count, ";")[[1]]
  ehdn.count <- data.frame(table(ehdn.count))
  names(ehdn.count)[2] <- "RareExpansions"
  
  fam.data <- fam.data[fam.data$Sample.ID %in% all, ]
  
  fam.data <- merge(fam.data, ehdn.count, by.x = "Sample.ID", by.y = "ehdn.count", all.x = T)
  fam.data$RareExpansions[is.na(fam.data$RareExpansions)] <- 0
  fam.data$lof <- fam.data$Sample.ID %in% lof
  
  print(nrow(fam.data))
  
  denovo.expansion <- sum(fam.data$lof & fam.data$RareExpansions > 0)
  denovo.noexpansion <- sum(fam.data$lof & fam.data$RareExpansions == 0)
  nodenovo.expansion <- sum(!fam.data$lof & fam.data$RareExpansions > 0)
  nodenovo.noexpansion <- sum(!fam.data$lof & fam.data$RareExpansions == 0)
  
  del.test <- fisher.test(data.frame("X1" = c(denovo.expansion, nodenovo.expansion), "X2" = c(denovo.noexpansion, nodenovo.noexpansion)))

  lof.out <- rbind(lof.out, 
                   data.frame(
                     "Dataset" = dataset,
                     "OR" = del.test$estimate,
                     "lower" = del.test$conf.int[1],
                     "upper" = del.test$conf.int[2],
                     "p" = del.test$p.value
                   ))
  
 
}

del <- ggplot(del.out, aes(x = Dataset, y = OR)) + geom_point() + geom_errorbar(aes(ymin = lower, ymax = upper), width = .1) + 
  geom_hline(yintercept = 1, lty = 2) + ylab("Odds ratio of having rare expansions") + theme_classic() +
  ggtitle("De novo deletion")
lof <- ggplot(lof.out, aes(x = Dataset, y = OR)) + geom_point() + geom_errorbar(aes(ymin = lower, ymax = upper), width = .1) + 
  geom_hline(yintercept = 1, lty = 2) + ylab("Odds ratio of having rare expansions") + theme_classic() +
  ggtitle("De novo LoFs")
plot_grid(del, lof, nrow = 1)

