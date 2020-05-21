### Analysis for x-linked expansions, separately for male and female

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(data.table)
library(GenomicRanges)
library(pwr)
library(rms)

source("requiredFunctions.R")

annotation <- fread("merged.expansion.coordinates_forAnno.txt.annovar.out_rev38.1mFmt.tsv", stringsAsFactors = F)
annotation <- data.frame(annotation)
annotation$key <- paste(annotation$chr, annotation$start, annotation$end, sep="#")
annotation$typeseq_priority <- factor(annotation$typeseq_priority, 
                                      levels = c("exonic", "splicing", "exonic;splicing", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_exonic;ncRNA_splicing", 
                                                 "UTR5", "UTR3", 
                                                 "intronic", "ncRNA_intronic", "upstream", "downstream", "intergenic"))
annotation <- annotation[order(annotation$typeseq_priority), ]
annotation <- annotation[!duplicated(annotation$key), ]


load("gs.asd.2019.RData")
gsASD <- gsASD[-33]
for(sex in c("male", "female")){
      
      ehdn <- fread("output_regions.min2.1000G+SSC+MSSNG.txt") #chr, start, end, motif, v1, v2, genotype
      trf <- fread("genotype_table.1000G+MSSNG+SSC.txt") #sample, relation, affection, dataset, region1
      
      ssc.errors1 <- read.delim('ssc_ME.problems.txt',
                                header = F, stringsAsFactors = F)
      ssc.errors2 <- read.delim('ssc_gender.problems.txt',
                                header = F, sep=" ", stringsAsFactors = F)
      ssc.error <- c(ssc.errors1$V2[-27], ssc.errors2$V7, ssc.errors2$V11)
      ssc.error <- unique(ssc.error[!ssc.error %in% c("", "0")])
      all.error <- c(ssc.error, "2-1547-003", "3-0051-000", "3-0008-000", "7-0064-003", "7-0075-003", "3-0395-000")
      
      brett.samples <- read.delim("manifest.1000G+SSC+MSSNG.txt", stringsAsFactors = F, header = F)
      brett.samples <- brett.samples[!brett.samples$V1 %in% all.error, ]
      
      ehdn <- as.data.frame(ehdn)
      trf <- as.data.frame(trf)
      ehdn$repeatID <- paste(ehdn$V1, ehdn$V2, ehdn$V3, ehdn$V4, sep="#")
      names(trf)[-c(1:4)] <- gsub(":|-", "#", names(trf)[-c(1:4)])
      
      tmp.outlier <- read.delim("MSSNG.samples.with.excessive.no.detected.EHdn.both.sides.tsv", stringsAsFactors = F)
      
      fam.data.tmp <- read.delim("Data_Merged_MSSNG.SSC.1KG.tsv", stringsAsFactors = F)
      fam.data.tmp <- fam.data.tmp[fam.data.tmp$Sample.ID %in% brett.samples$V1, ]

      fam.data.tmp <- fam.data.tmp[!fam.data.tmp$Sample.ID %in% tmp.outlier$Sample, ]
      fam.1000g <- fam.data.tmp[fam.data.tmp$Dataset == "1000G", ]

      fam.data <- fam.data.tmp[which(!is.na(fam.data.tmp$Sample.ID)), ]
      fam.data <- fam.data[fam.data$Status != "", ]
      ### limit to one affected kid per family
      fam.data <- fam.data[which((!duplicated(paste(fam.data$Family.ID, fam.data$Status)) & fam.data$Status == "AffectedKid") |
                                   fam.data$Status != "AffectedKid"), ]
      fam.data <- fam.data[fam.data$DNASOURCE == "Other cell", ]
      fam.data <- fam.data[fam.data$Library.type == "PCR-free" & fam.data$Platform == "Illumina HiSeqX", ]
      fam.data <- fam.data[fam.data$Status %in% c("AffectedKid", "UnaffectedKid"), ]
      fam.data <- fam.data[fam.data$Sex == sex, ]
      
      ehdn.result <- read.delim(sprintf("outlier.EHdn.chrx.%s.only.tsv", sex), stringsAsFactors = F)
      ehdn.result$repeatID <- paste(ehdn.result$chr, ehdn.result$start, ehdn.result$end, ehdn.result$motif, sep="#")
      
      map <- read.delim("map.TRF.EHdn.0.66.tsv", stringsAsFactors = F)
      ehdn.result <- merge(ehdn.result, map, by.x = "repeatID", by.y = "ehdn.id", all.x = T)
      ehdn.result$freq_1000g <- sapply(ehdn.result$outliers, get1000Freq, fam.1000g)
      ehdn.result <- ehdn.result[ehdn.result$freq_1000g <= 0.001, ]
      
      ehdn.all <- merge(ehdn[, -c(5:7, 9:12)], map, by.x = "repeatID", by.y = "ehdn.id", all.x = T)
      names(ehdn.all) <- c("repeatID", "chr", "start", "end", "motif", "trf.id", "reciprocalIdentity", "length.diff")
      ehdn.all <- ehdn.all[ehdn.all$chr == "chrX", ]
      ehdn.m <- mergeLoci(ehdn.all)
      ehdn.m$uniqueMotif <- sapply(sapply(ehdn.m$motif, strsplit, ";"), length)
      ehdn.result.m <- mapOutlierToMergeLoci(ehdn.result, ehdn.m)
      ehdn.m$hasOutlier <- ehdn.m$repeatID %in% ehdn.result.m$repeatID
      
      ehdn.m$allRepeat <- ehdn.m$repeatID
      ehdn.m$repeatID <- paste(ehdn.m$chr, ehdn.m$start, ehdn.m$end, sep="#")
      ehdn.result.m$allRepeat <- ehdn.result.m$repeatID
      ehdn.result.m$repeatID <- paste(ehdn.result.m$chr, ehdn.result.m$start, ehdn.result.m$end, sep="#")
      
      start <- Sys.time()
      tmp <- do.call(rbind, lapply(1:nrow(ehdn.result.m), getOutlierMergeData, ehdn.result.m, fam.data, fam.1000g))
      end <- Sys.time()
      difftime(end, start)

      names(tmp)[1] <- "Sample"
      tmp.f <- tmp[tmp$Dataset != "1000G", ]
      
      tmp.f <- tmp.f[tmp.f$chr == "chrX", ]
      annotation.f <- annotation$key[!is.na(annotation$gene_symbol) & annotation$key %in% tmp.f$repeatID &
                                       annotation$typeseq_priority %in% c("exonic", "splicing", "exonic;splicing", 
                                                                          "UTR5", "UTR3", 
                                                                          "intronic", "upstream", "downstream")]
      annotation.intergenic <- annotation$key[annotation$typeseq_priority == "intergenic"& annotation$key %in% tmp.f$repeatID]
      
      #### non coding loci
      ehdn.genotype.data <- strsplit(paste(ehdn.result.m$outliers[ehdn.result.m$chr %in% "chrX" & 
                                                                    ehdn.result.m$repeatID %in% annotation.intergenic], 
                                           collapse = ";"), ";")[[1]]
      
      ehdn.genotype.data <- data.frame(table(ehdn.genotype.data),
                                          stringsAsFactors = F)
      names(ehdn.genotype.data) <- c("Sample", "EHdnGenotypeNCLoci")
      
      
      intergenic <- merge(ehdn.genotype.data, fam.data, by.x = "Sample", by.y = "Sample.ID", all.y = T)
      intergenic$EHdnGenotypeNCLoci[is.na(intergenic$EHdnGenotypeNCLoci)] <- 0
      intergenic <- intergenic[which(intergenic$Status %in% c("AffectedKid", "UnaffectedKid")), ]
      intergenic$Status <- factor(intergenic$Status, levels = c("UnaffectedKid", "AffectedKid"))
      
      tmp.tmp <- tmp.f
      
      
      cov = c("EHdnGenotypeNCLoci")
      tmp.f <- tmp.tmp[which(tmp.tmp$repeatID %in% annotation.f), ]
      tmp.agr.res <- getGenesetTable(unique(tmp.f$repeatID[which(tmp.f$repeatID %in% annotation.f)]), 
                                     tmp.f, fam.data[fam.data$Dataset == "SSC" & fam.data$ETH_TAG == "EUR", ], 
                                     ehdn.genotype.data, cov, "global")
      
      tmp.agr.res <- merge(tmp.agr.res, fam.data, by.x = "Sample", by.y = "Sample.ID", all.x = T)
      tmp.agr.res <- merge(tmp.agr.res, ehdn.genotype.data, by = "Sample", all.x =T)
      tmp.agr.res$EHdnGenotypeNCLoci[is.na(tmp.agr.res$EHdnGenotypeNCLoci)] <- 0
      tmp.agr.res$Status <- factor(tmp.agr.res$Status, levels = c("UnaffectedKid", "AffectedKid"))
      lm.ref <- glm("Status ~ EHdnGenotypeNCLoci", tmp.agr.res, family = binomial(link = "logit"))
      lm.add <- glm("Status ~ EHdnGenotypeNCLoci + global", tmp.agr.res, family = binomial(link = "logit"))
      ano <- anova(lm.ref, lm.add, test = "Chisq")
      
      or <- exp(lm.add$coefficients)["global"]
      pvalue <- ano$`Pr(>Chi)`[2]

      lower <- confint(lm.add)["global", 1]
      upper <- confint(lm.add)["global", 2]
      
      n <- sum(fam.data$Dataset == "SSC" & fam.data$ETH_TAG == "EUR")
      u <- 1 # number of coefficients in full model excepts intercept and covariates
      w <- 1 # number of covariates
      
      pR1 = 1 - lm.ref$deviance / lm.ref$null.deviance 
      pR2 = 1 - lm.add$deviance / lm.add$null.deviance # R2 based on deviance
      
      f2 <- (pR2 - pR1) / (1 - pR2)#based on https://www.statmethods.net/stats/power.html two models
      
      # current.power <- pwr.f2.test(u = u, v = n - w - u - 1, f2 = f2, sig.level = pvalue)$power
      # new.power <- pwr.f2.test(u = u, v = 815 - w - u - 1, f2 = f2, sig.level = pvalue)$power
      new.n <- pwr.f2.test(u = u, power = 0.8, f2 = f2, sig.level = 0.05)$v
      new.n <- new.n + w + u + 1 
      
      tmp.agr.res <- data.frame("global", or, exp(lower), exp(upper), pvalue, new.n)
      write.table(tmp.agr.res, sprintf("global.any.motif.genic.%s.tsv", sex), sep="\t", row.names=F, quote = F, col.names=T)
      
      ###########################################
      ########### gene set analysis #############
      cov = c("EHdnGenotypeNCLoci")
      
      gene.info <- annotation[!is.na(annotation$entrez_id), ]
      load("gs.asd.2019.RData")
      gsASD <- gsASD[-33]
      
      tmp.gs <- data.frame()
      type <- "genicCorrectedByNC"
      for(gs in names(gsASD)[!names(gsASD) %in% c("PhMm_Aggr_SkeCranioLimbs_all", "FMR1_Targets_Ascano")]){
        print(gs)
        flush.console()
        ann.tmp <- gene.info$key[gene.info$entrez_id %in% as.numeric(gsASD[[gs]])]
        
        ann.tmp <- ann.tmp[grep("chrX", ann.tmp)]
        if(length(ann.tmp) > 0){
          if(sum(tmp.f$repeatID %in% ann.tmp & tmp.f$Sample %in% fam.data$Sample.ID[fam.data$Dataset == "SSC" & fam.data$ETH_TAG == "EUR"]) == 0){
            
          }else{
            tmp.gs.res <- testByAggregate(unique(tmp.f$repeatID[tmp.f$repeatID %in% ann.tmp]), 
                                          tmp.f[tmp.f$repeatID %in% ann.tmp, ],fam.data[fam.data$Dataset == "SSC" & fam.data$ETH_TAG == "EUR", ],
                                          ehdn.genotype.data, cov, standardization = F)
            tmp.gs.res$geneset = gs
            tmp.gs.res$type = type
          
            tmp.gs <- rbind(tmp.gs, tmp.gs.res)
          }
        }
      }
      
      tmp.gs$aggr.all.bhfdr <- 1
      tmp.gs$aggr.all.fwer <- 1
      tmp.gs$aggr.all.bhfdr <- p.adjust(tmp.gs$all.aggr.glm.pvalue, method = "BH")
      tmp.gs$aggr.all.fwer <- p.adjust(tmp.gs$all.aggr.glm.pvalue, method = "bonferroni")

      write.table(tmp.gs, sprintf("geneset.any.motif.%s.tsv", sex),
                  sep = "\t", row.names=F, quote=F, col.names=T)

}

