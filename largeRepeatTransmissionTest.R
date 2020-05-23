## This script is to perform transmission test for large tandem repeats 
## by comparing the rates of transmission of certain element to other large tandem repeats

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(cowplot)

source("requiredFunctions.R")
load("gs.asd.2019.RData")
gsASD <- gsASD[-33]
annotation <- fread("merged.expansion.coordinates_forAnno.txt.annovar.out_rev38.1mFmt.tsv", stringsAsFactors = F)
annotation <- data.frame(annotation)
annotation <- annotation[annotation$chr %in% paste0("chr", 1:22), ]
annotation$key <- paste(annotation$chr, annotation$start, annotation$end, sep="#")
annotation.intergenic <- annotation[annotation$typeseq_priority ==  "intergenic", ]
annotation.all <- annotation
annotation <- annotation[annotation$typeseq_priority %in% c("exonic", "splicing", "exonic;splicing", 
                                                            "UTR5", "UTR3", "intronic", "upstream", "downstream"), ]
ehdn <- fread("output_regions.min2.1000G+SSC+MSSNG.txt") #chr, start, end, motif, v1, v2, genotype
ehdn$key <- paste(ehdn$V1, ehdn$V2, ehdn$V3, ehdn$V4, sep="#")

annotation$typeseq_priority <- factor(annotation$typeseq_priority, 
                                      levels = c("exonic", "splicing", "exonic;splicing", 
                                                 "UTR5", "UTR3", "intronic", "upstream", "downstream"))
annotation <- annotation[order(annotation$typeseq_priority), ]
annotation <- annotation[!duplicated(annotation$key), ]

all.error <- readLines("samplewitherror.txt")
brett.samples <- read.delim("manifest.1000G+SSC+MSSNG.txt", stringsAsFactors = F, header = F)

outlier <- read.delim("MSSNG.samples.with.excessive.no.detected.EHdn.both.sides.ASD.tsv", stringsAsFactors = F)
ehdn.result.m <- read.delim("merged.any.outlier.expansions.autosomal.sex.chr.nobadsamples.tsv", stringsAsFactors = F)

fam.data.tmp <- read.delim(".Data_Merged_MSSNG.SSC.1KG.tsv", stringsAsFactors = F)
fam.data.tmp <- fam.data.tmp[fam.data.tmp$Sample.ID %in% brett.samples$V1, ]
fam.data <- fam.data.tmp[!fam.data.tmp$Sample.ID %in% all.error, ]
fam.data <- fam.data[!fam.data.tmp$Sample.ID %in% outlier$Sample, ]
fam.1000g <- fam.data.tmp[fam.data.tmp$Dataset == "1000G", ]
fam.data <- fam.data[which(!is.na(fam.data$Sample.ID)), ]
fam.data <- fam.data[fam.data$Status != "", ]
### limit to one affected kid per family
fam.data <- fam.data[which((!duplicated(paste(fam.data$Family.ID, fam.data$Status)) & fam.data$Status == "AffectedKid") |
                             fam.data$Status != "AffectedKid"), ]

fam.data <- fam.data[fam.data$DNASOURCE == "Other cell", ]
fam.data <- fam.data[fam.data$Library.type == "PCR-free" & fam.data$Platform == "Illumina HiSeqX", ]
fam.data <- fam.data[fam.data$ETH_TAG == "EUR" & fam.data$Dataset == "SSC", ]

dt.out <- data.frame()

### this needs to be run three times for SSC, MSSNG and SSC+MSSNG
### please change condition in line 50 to get different set of samples
for(i in 1:nrow(ehdn.result.m)){
  this.locus <- ehdn.result.m[i, ]
  
  outliers <- strsplit(this.locus$outliers, ";")[[1]]
  outliers <- outliers[outliers %in% fam.data$Sample.ID]
  
  if(length(outliers) > 0){
    this.ehdn <- ehdn[ehdn$key %in% strsplit(this.locus$allRepeat, ";")[[1]], ]
    
    tmp <- paste(this.ehdn$V7, collapse = ",")
    tmp <- gsub(":|,", "#", tmp)
    tmp <- unlist(strsplit(tmp, "#")[[1]])
    tmp <- data.frame(matrix(tmp, ncol = 2, byrow = T), stringsAsFactors = F)
    tmp[, 2] <- as.numeric(tmp[, 2])
    tmp <- aggregate(X2 ~ X1, tmp, max)
    
    ref <- as.numeric(names(which.max(table(round(tmp$X2, digits = 1)))))
    
    tmp <- tmp[tmp$X1 %in% c(fam.data$Sample.ID, fam.1000g$Sample.ID), ]
    tmp <- tmp[order(tmp$X2, decreasing = T), ]  
    tmp$order <- 1:nrow(tmp)
    tmp <- tmp[tmp$order < (nrow(fam.data) + nrow(fam.1000g))*0.01 & tmp$X2 > ref*2 & tmp$X2 > 2, ]
    
    this.fam <- fam.data[fam.data$Sample.ID %in% tmp$X1, ]
    
    parent <- this.fam$Family.ID[this.fam$Relation %in% c("father", "mother")]
    affected <- sum(fam.data$Sample.ID %in% this.fam$Sample.ID &
                           fam.data$Status == "AffectedKid" & fam.data$Family.ID %in% parent)
    unaffected <- sum(fam.data$Sample.ID %in% this.fam$Sample.ID &
                           fam.data$Status == "UnaffectedKid" & fam.data$Family.ID %in% parent)
    
    possible.affected <- sum(fam.data$Status == "AffectedKid" & fam.data$Family.ID %in% parent)
    possible.unaffected <- sum(fam.data$Status == "UnaffectedKid" & fam.data$Family.ID %in% parent)
    
    dt.out <- rbind(dt.out, data.frame(
      this.locus, affected, unaffected, possible.affected, possible.unaffected #, male, female, possible.male, possible.female
    ))
  }
}

write.table(dt.out, "transmitted.ssc.expansions.99perc.tsv", sep="\t", row.names=F, quote=F, col.names=T)


for(cohort in c("ssc", "mssng", "eur")){
  dt.out <- read.delim(sprintf("transmitted.%s.expansions.99perc.tsv", cohort), stringsAsFactors = F)
  ### gene-set
  test.out <- data.frame()
  for(igs in 1:length(gsASD)){

      gs <- names(gsASD)[igs]
      tmp.dt <- dt.out[dt.out$repeatID %in% annotation$key[
        annotation$entrez_id %in% gsASD[[gs]]], ]
      
      affected <- sum(tmp.dt$affected)
      affected.denom <- sum(tmp.dt$possible.affected) - affected
  
      tmp.dt <- dt.out[!dt.out$repeatID %in% annotation$key[
        annotation$entrez_id %in% gsASD[[gs]]], ]
      
      unaffected <-  sum(tmp.dt$affected)
      unaffected.denom <- sum(tmp.dt$possible.affected) - unaffected
      test <- fisher.test(data.frame("V1" = c(affected, unaffected), "V2" = c(affected.denom, unaffected.denom)))

      test.out <- rbind(test.out, data.frame("gs.query" = gs, #"gs.target" = gs.tar, 
                                             "or" = test$estimate, "lower" = test$conf.int[1], "upper" = test$conf.int[2], "pvalue" = test$p.value,
                                             "transmitted.case" = affected, "transmitted.ctrl" = unaffected,
                                             "case.denom" = affected.denom, "ctrl.denom" = unaffected.denom, "nrepeat" = nrow(tmp.dt)))
  }

  test.out$FWER <- p.adjust(test.out$pvalue, method = "bonferroni")
  write.table(test.out[order(test.out$pvalue), ], 
              sprintf("%s.affected.99perc.geneset.remained_repeats_as_background.tsv",cohort),
              sep="\t", row.names=F, quote=F, col.names=T)
  
  dt.out <- read.delim(sprintf("transmitted.%s.expansions.99perc.tsv", cohort), stringsAsFactors = F)
  ### gene-set
  test.out <- data.frame()
  
  functions <- list()
  functions[["intergenic"]] <- annotation.intergenic$key
  functions[["genic"]] <- annotation$key[annotation$typeseq_priority %in% c("exonic", "exonic;splicing", "splicing", "upstream", 
                                                                               "UTR3", "UTR5", "downstream", "intronic")]
  functions[["upstream"]] <- annotation$key[annotation$typeseq_priority == "upstream"]
  functions[["downstream"]] <- annotation$key[annotation$typeseq_priority == "downstream"]
  functions[["intron"]] <- annotation$key[annotation$typeseq_priority == "intronic"]
  functions[["UTR5"]] <- annotation$key[annotation$typeseq_priority == "UTR5"]
  functions[["UTR3"]] <- annotation$key[annotation$typeseq_priority == "UTR3"]
  functions[["exon"]] <- annotation$key[annotation$typeseq_priority %in% c("exonic")]
  functions[["splicing"]] <- annotation$key[(annotation$typeseq_priority %in% c("splicing", "exonic;splicing"))]

  for(igs in 1:length(functions)){
    
    gs <- names(functions)[igs]
    
    tmp.dt <- dt.out[dt.out$repeatID %in% functions[[igs]], ]
    
    affected <- sum(tmp.dt$affected)
    affected.denom <- sum(tmp.dt$possible.affected) - affected
    
    tmp.dt <- dt.out[!dt.out$repeatID %in% functions[[igs]], ]
    unaffected <-  sum(tmp.dt$affected)
    unaffected.denom <- sum(tmp.dt$possible.affected) - unaffected
    test <- fisher.test(data.frame("V1" = c(affected, unaffected), "V2" = c(affected.denom, unaffected.denom)))
    
    test.out <- rbind(test.out, data.frame("gs.query" = gs, #"gs.target" = gs.tar, 
                                           "or" = test$estimate, "lower" = test$conf.int[1], "upper" = test$conf.int[2], "pvalue" = test$p.value,
                                           "transmitted.case" = affected, "transmitted.ctrl" = unaffected,
                                           "case.denom" = affected.denom, "ctrl.denom" = unaffected.denom, "nrepeat" = nrow(tmp.dt)))
  }
  
  test.out$FWER <- p.adjust(test.out$pvalue, method = "bonferroni")
  write.table(test.out[order(test.out$pvalue), ], 
              sprintf("%s.affected.99perc.element.remained_repeats_as_background.tsv",cohort),
              sep="\t", row.names=F, quote=F, col.names=T)
  
      ###unaffected
      if(cohort == "ssc"){
      dt.out <- read.delim(sprintf("transmitted.%s.expansions.99perc.tsv", cohort), stringsAsFactors = F)
      test.out <- data.frame()
      for(igs in 1:length(gsASD)){

        gs <- names(gsASD)[igs]
        tmp.dt <- dt.out[dt.out$repeatID %in% annotation$key[
          annotation$entrez_id %in% gsASD[[gs]]], ]

        affected <- sum(tmp.dt$unaffected)
        affected.denom <- sum(tmp.dt$possible.unaffected) - affected
        
        tmp.dt <- dt.out[!dt.out$repeatID %in% annotation$key[
          annotation$entrez_id %in% gsASD[[gs]]], ]
        unaffected <-  sum(tmp.dt$unaffected)
        unaffected.denom <- sum(tmp.dt$possible.unaffected) - unaffected
        test <- fisher.test(data.frame("V1" = c(affected, unaffected), "V2" = c(affected.denom, unaffected.denom)))

        test.out <- rbind(test.out, data.frame("gs.query" = gs, #"gs.target" = gs.tar,
                                               "or" = test$estimate, "lower" = test$conf.int[1], "upper" = test$conf.int[2], "pvalue" = test$p.value,
                                               "transmitted.case" = affected, "transmitted.ctrl" = unaffected,
                                               "case.denom" = affected.denom, "ctrl.denom" = unaffected.denom, "nrepeat" = nrow(tmp.dt)))
      }

      test.out$FWER <- p.adjust(test.out$pvalue, method = "bonferroni")
      write.table(test.out[order(test.out$pvalue), ], 
                  sprintf("%s.unaffected.99perc.geneset.remained_repeats_as_background.tsv",cohort),
                  sep="\t", row.names=F, quote=F, col.names=T)

      dt.out <- read.delim(sprintf("transmitted.%s.expansions.99perc.tsv", cohort), stringsAsFactors = F)
      ### gene-set
      test.out <- data.frame()

      functions <- list()
      functions[["intergenic"]] <- annotation.intergenic$key
      functions[["genic"]] <- annotation$key[annotation$typeseq_priority %in% c("exonic", "exonic;splicing", "splicing", "upstream",
                                                                                "UTR3", "UTR5", "downstream", "intronic")]
      functions[["upstream"]] <- annotation$key[annotation$typeseq_priority == "upstream"]
      functions[["downstream"]] <- annotation$key[annotation$typeseq_priority == "downstream"]
      functions[["intron"]] <- annotation$key[annotation$typeseq_priority == "intronic"]
      functions[["UTR5"]] <- annotation$key[annotation$typeseq_priority == "UTR5"]
      functions[["UTR3"]] <- annotation$key[annotation$typeseq_priority == "UTR3"]
      functions[["exon"]] <- annotation$key[annotation$typeseq_priority %in% c("exonic")]
      functions[["splicing"]] <- annotation$key[(annotation$typeseq_priority %in% c("splicing", "exonic;splicing"))]

      for(igs in 1:length(functions)){

        gs <- names(functions)[igs]

        tmp.dt <- dt.out[dt.out$repeatID %in% functions[[igs]], ]

        affected <- sum(tmp.dt$unaffected)
        affected.denom <- sum(tmp.dt$possible.unaffected) - affected
        
        tmp.dt <- dt.out[!dt.out$repeatID %in% functions[[igs]], ]
        unaffected <-  sum(tmp.dt$unaffected)
        unaffected.denom <- sum(tmp.dt$possible.unaffected) - unaffected
        test <- fisher.test(data.frame("V1" = c(affected, unaffected), "V2" = c(affected.denom, unaffected.denom)))

        test.out <- rbind(test.out, data.frame("gs.query" = gs, #"gs.target" = gs.tar,
                                               "or" = test$estimate, "lower" = test$conf.int[1], "upper" = test$conf.int[2], "pvalue" = test$p.value,
                                               "transmitted.case" = affected, "transmitted.ctrl" = unaffected,
                                               "case.denom" = affected.denom, "ctrl.denom" = unaffected.denom, "nrepeat" = nrow(tmp.dt)))
      }

      test.out$FWER <- p.adjust(test.out$pvalue, method = "bonferroni")
      write.table(test.out[order(test.out$pvalue), ], 
                  sprintf("%s.unaffected.99perc.element.remained_repeats_as_background.tsv",cohort),
                  sep="\t", row.names=F, quote=F, col.names=T)
      }
}  