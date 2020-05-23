### This is the script to perform main analysis shown in the paper.
### It includes global burden test, gene-set burden test and genic element burden test.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(cowplot)

source("requiredFunctions.R")

load("gs.asd.2019.RData")

annotation <- fread("merged.expansion.coordinates_forAnno.txt.annovar.out_rev38.1mFmt.tsv", stringsAsFactors = F)
annotation <- data.frame(annotation)
annotation$key <- paste(annotation$chr, annotation$start, annotation$end, sep="#")
annotation$typeseq_priority <- factor(annotation$typeseq_priority, 
                                      levels = c("exonic", "splicing", "exonic;splicing", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_exonic;ncRNA_splicing", 
                                                 "UTR5", "UTR3", "intronic", "ncRNA_intronic", "upstream", "downstream", "intergenic"))
annotation <- annotation[order(annotation$typeseq_priority), ]
annotation <- annotation[!duplicated(annotation$key), ]

ehdn <- fread("output_regions.min2.1000G+SSC+MSSNG.txt") #chr, start, end, motif, v1, v2, genotype
trf <- fread("genotype_table.1000G+MSSNG+SSC.txt") #sample, relation, affection, dataset, region1
tmp.count <- read.delim("call_counts.txt")

tmp.1000g <- tmp.count[tmp.count$Dataset.platform.library.type %in% c("1000G/Illumina NovaSeq/PCR-free"), ]
tmp.count <- tmp.count[tmp.count$Dataset.platform.library.type %in% c("SSC/Illumina HiSeq X/PCR-free", "MSSNG/Illumina HiSeq X/PCR-free"), ]
mean.1000g <- mean(tmp.1000g$EHDN.call.count)
sd.1000g <- sd(tmp.1000g$EHDN.call.count)
tmp.1000g.outlier <- tmp.1000g[tmp.1000g$EHDN.call.count > mean.1000g+3*sd.1000g |
                                 tmp.1000g$EHDN.call.count < mean.1000g-3*sd.1000g, ]

mean <- mean(tmp.count$EHDN.call.count)
sd <- sd(tmp.count$EHDN.call.count)
tmp.outlier <- tmp.count[tmp.count$EHDN.call.count > mean+3*sd |
                           tmp.count$EHDN.call.count < mean-3*sd, ] #74 outliers

tmp.outlier <- rbind(tmp.outlier, tmp.count[grep("AU2218", tmp.count$Sample), ])

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

fam.setting <- c()
count.complete <- 0
singletons <- 0
mssng.fam <- c()
mssng.sin <- c()
for(fam in unique(fam.data$Family.ID[fam.data$Dataset == "MSSNG"])){
  tmp <- fam.data[fam.data$Family.ID == fam & fam.data$Dataset == "MSSNG", ]
  #tmp$Status <- tmp$Status.x
  if(sum(tmp$Status %in% c("AffectedParent", "UnaffectedParent")) == 2 & sum(tmp$Status == "AffectedKid") > 0 ){
    count.complete <- count.complete + 1
    mssng.fam <- c(mssng.fam, fam.data$Sample.ID[which(fam.data$Family.ID == fam & fam.data$Status == "AffectedKid" & fam.data$Dataset == "MSSNG")])
  }
  else if(sum(tmp$Status == "AffectedKid") > 0 ){
    singletons <- singletons + 1 #nrow(tmp)
    mssng.sin <- c(mssng.sin, fam.data$Sample.ID[which(fam.data$Family.ID == fam & fam.data$Status == "AffectedKid" & fam.data$Dataset == "MSSNG")])
  }
  tmp <- paste(sort(tmp$Status), collapse = ";")

  fam.setting  <- c(fam.setting, tmp)
}

fam.data <- fam.data[fam.data$Status %in% c("AffectedKid", "UnaffectedKid"), ]
fam.data <- fam.data[fam.data$Sample.ID %in% tmp.count$Sample, ]
fam.data <- merge(fam.data, tmp.count[, c("Sample", "EHDN.call.count")], by.x = "Sample.ID", by.y = "Sample", all.x = T)


map <- read.delim("map.TRF.EHdn.0.66.tsv", stringsAsFactors = F)

clean.ehdn.m <- read.delim( "merged.EHdn.genotype.samples.without.outliers.tsv", stringsAsFactors = F)

ehdn.result <- read.delim("outlier.EHdn.MSSNG.SSC.1000G.tsv", stringsAsFactors = F)

ehdn.result$repeatID <- paste(ehdn.result$chr, ehdn.result$start, ehdn.result$end, ehdn.result$motif, sep="#")

ehdn.result <- merge(ehdn.result, map, by.x = "repeatID", by.y = "ehdn.id", all.x = T)
ehdn.result$freq_1000g <- sapply(ehdn.result$outliers, get1000Freq, fam.1000g)
ehdn.result <- ehdn.result[ehdn.result$freq_1000g <= 0.001, ]

ehdn.all <- merge(ehdn[, -c(5:7, 9:12)], map, by.x = "repeatID", by.y = "ehdn.id", all.x = T)
names(ehdn.all) <- c("repeatID", "chr", "start", "end", "motif", "trf.id", "reciprocalIdentity", "length.diff")
ehdn.all <- ehdn.all[ehdn.all$chr %in% paste0("chr", c(1:22, "X", "Y")), ]
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
write.table(tmp, "tabulated.EHdn.outliers.overlapped.any.tsv", sep="\t", row.names=F, quote=F, col.names=T)
#################################################################################################
# 
tmp <- read.delim("tabulated.EHdn.outliers.overlapped.any.tsv", stringsAsFactors = F)
names(tmp)[1] <- "Sample"

tmp.f <- tmp[tmp$Dataset != "1000G", ]
tmp.f <- tmp.f[tmp.f$chr %in% paste0("chr", c(1:22, "X")), ]
# # 
start <- Sys.time()
tmp.res <- do.call(rbind, lapply(unique(tmp.f$repeatID), testByLocus, tmp.f, fam.data))
end <- Sys.time()
difftime(end, start)

write.table(tmp.res, "tabulated.merge.any.EHdn.outliers.with.tests.tsv", sep="\t", row.names=F, quote=F, col.names=T)

####
tmp.res <- read.delim("tabulated.merge.any.EHdn.outliers.with.tests.tsv", stringsAsFactors = F)
annotation.tmp <- annotation

tmp.res$freq_affected <- tmp.res$AGGR.AffectedKid/tmp.res$AGGR.AffectedKid.denominator
tmp.res$freq_unaffected <- tmp.res$AGGR.UnaffectedKid/tmp.res$AGGR.UnaffectedKid.denominator
tmp.res$freq_1000g <- tmp.res$REF1000G/nrow(fam.1000g)

tmp.annotated <- merge(tmp.res, annotation, by.x = "repeatID", by.y = "key", all.x = T)
tmp.annotated$insigset <- tmp.annotated$entrez_id %in% c(gsASD$PhMm_Aggr_CardvascMuscle_all,
                                                         gsASD$Neurof_GoNervSysDev)

tmp.sup <- tmp.annotated[tmp.annotated$freq_affected > pmax(tmp.annotated$freq_unaffected, tmp.annotated$freq_1000g, na.rm = T) &
                           tmp.annotated$insigset & tmp.annotated$chr.x %in% paste0("chr", c(1:22)), c(1:12, 19:21, 32, 37)]
tmp.sup$mssng.fam <- 0
tmp.sup$mssng.sin <- 0

for(i in 1:nrow(tmp.sup)){
  samples <- strsplit(as.character(tmp.sup$samples[i]), ";")[[1]]
  tmp.sup$mssng.fam[i] <- sum(samples %in% mssng.fam)
  tmp.sup$mssng.sin[i] <- sum(samples %in% mssng.sin)
}

tmp.sup <- tmp.sup[, c(1, 3:5, 2, 6, 18:19, 8:10, 13:15, 16:17)]
names(tmp.sup) <- c("repeatID", "chr", "start", "end", "motif", "TRF_loci", "affected_MSSNG_family", "affected_MSSNG_singleton", "affected_SSC", "1000_genomes",
                    "affected_outliers", "freq_affected", "freq_unaffected", "freq_1000_genomes", "gene_symbol", "omim_phenotype")

tmp.table2 <- tmp.annotated[tmp.annotated$freq_affected > 0.001 & tmp.annotated$freq_unaffected < 0.001 &
                              tmp.annotated$freq_1000g < 0.001 &
                              tmp.annotated$insigset, ]

write.table(tmp.sup, "supplementary.table.tsv", sep="\t", row.names=F, quote=F, col.names=T)
write.table(tmp.table2, "top.loci.tsv", sep="\t", row.names=F, quote=F, col.names=T)
####
tmp <- read.delim("tabulated.EHdn.outliers.overlapped.any.tsv", stringsAsFactors = F)
names(tmp)[1] <- "Sample"

tmp.f <- tmp[tmp$Dataset != "1000G", ]
tmp.f <- tmp.f[tmp.f$chr %in% paste0("chr", c(1:22)), ]

annotation.f <- annotation$key[!is.na(annotation$gene_symbol) & annotation$key %in% tmp.f$repeatID &
                                 annotation$typeseq_priority %in% c("exonic", "exonic;splicing", "upstream", "downstream",
                                                                    "UTR3", "UTR5", "splicing", "intronic")]

annotation.intergenic <- annotation$key[annotation$typeseq_priority == "intergenic" & annotation$key %in% tmp.f$repeatID]

#### non coding loci
ehdn.genotype.data <- strsplit(paste(ehdn.result.m$outliers[ehdn.result.m$chr %in% paste0("chr", 1:22) & 
                                                              ehdn.result.m$repeatID %in% annotation.intergenic], 
                                     collapse = ";"), ";")[[1]]

ehdn.genotype.data <- data.frame(table(ehdn.genotype.data),
                                    stringsAsFactors = F)
names(ehdn.genotype.data) <- c("Sample", "EHdnGenotypeNCLoci")

ehdn.all.data <- strsplit(paste(ehdn.result.m$outliers[ehdn.result.m$chr %in% paste0("chr", 1:22)],
                                collapse = ";"), ";")[[1]]
ehdn.all.data <- data.frame(table(ehdn.all.data), stringsAsFactors = F)
names(ehdn.all.data) <- c("Sample", "EHdnGenotypeNCLoci")

### genic
ehdn.genic <- strsplit(paste(ehdn.result.m$outliers[ehdn.result.m$chr %in% paste0("chr", 1:22) & 
                                                      ehdn.result.m$repeatID %in% annotation.f], 
                             collapse = ";"), ";")[[1]]

ehdn.genic <- data.frame(table(ehdn.genic),
                         stringsAsFactors = F)
names(ehdn.genic) <- c("Sample", "EHdnGenotypeNCLoci")


merge.data <- merge(ehdn.genotype.data, ehdn.genic, by = "Sample", all = T)
names(merge.data) <- c("Sample", "intergenic", "genic")
merge.data <- merge(merge.data, fam.data, by.x = "Sample", by.y = "Sample.ID", all.y = T)
merge.data[is.na(merge.data)] <- 0
merge.data$Status <- factor(merge.data$Status, levels = c("UnaffectedKid", "AffectedKid"))

mssng.rare <- sum((merge.data$intergenic[merge.data$Dataset == "MSSNG"] + merge.data$genic[merge.data$Dataset == "MSSNG"]) > 0)
ssc.rare <- sum((merge.data$intergenic[merge.data$Dataset == "SSC"] + merge.data$genic[merge.data$Dataset == "SSC"]) > 0)
mssng.all <- sum(merge.data$Dataset == "MSSNG")
ssc.all <- sum(merge.data$Dataset == "SSC")

chisq.test(data.frame("X1" = c(mssng.rare, ssc.rare),
                      "X2" = c(mssng.all, ssc.all)))
fisher.test(data.frame("X1" = c(mssng.rare, ssc.rare),
                      "X2" = c(mssng.all, ssc.all)))

s <- c()
for(i in c("intergenic", "V1_AMR", "V2_EUR", "V3_EAS", "V4_SAS", "V5_AFR")){
  lm <- glm(sprintf("Status ~ %s", i), merge.data, family = binomial(link = "logit"))
  s <- c(s, summary(lm)$coefficients[2, "Pr(>|z|)"])
}

s

lm <- glm(Status ~ genic, merge.data[merge.data$Dataset == "SSC", ], family = binomial(link = "logit"))
summary(lm)
exp(lm$coefficients["genic"])

ref <- glm(Status ~ 1, merge.data[merge.data$Dataset == "SSC", ], family = binomial(link = "logit"))
add <- glm(Status ~ 1 + genic, merge.data[merge.data$Dataset == "SSC", ], family = binomial(link = "logit"))
anova(ref, add, test = "Chisq")
summary(lm)
exp(add$coefficients["genic"])


lm <- glm(Status ~ V2_EUR + V4_SAS + V5_AFR, merge.data, family = binomial(link = "logit"))
summary(lm)

lm.ref <- glm(factor(Dataset) ~  V4_SAS + V5_AFR + V2_EUR , merge.data,
              family = binomial(link = "logit"))
lm.add <- glm(factor(Dataset) ~ V4_SAS + V5_AFR + V2_EUR + intergenic + genic, merge.data,
              family = binomial(link = "logit"))
anova(lm.ref, lm.add, test = "Chisq") #p=0.29
exp(lm.add$coefficients)

ggplot(merge.data[merge.data$Status == "AffectedKid", ], aes(x = genic, fill = Dataset)) + 
  geom_density(alpha = .5) + xlab("# rare tandem repeats")# + annotate("text", x = 1.5, y = 1.5, label = "p=0.24")
ggsave("ehdn.rare.repeat.count.mssng.ssc.png", width = 5, height = 5)

count <- data.frame()
mssng.sum <- 0
ssc.sum <- 0
for(i in 0:18){
  mssng.i <- sum(merge.data$intergenic[merge.data$Dataset == "MSSNG"] + merge.data$genic[merge.data$Dataset == "MSSNG"] == i)  
  ssc.i <- sum(merge.data$intergenic[merge.data$Dataset == "SSC"] + merge.data$genic[merge.data$Dataset == "SSC"] == i)
  
  mssng.sum <- mssng.sum + mssng.i
  ssc.sum <- ssc.sum + ssc.i
  
  count <- rbind(count, data.frame(
    "rare_expansion_count" = i,
    "cohort" = "MSSNG",
    "samples" = mssng.i, stringsAsFactors = F
  ))
  
  count <- rbind(count, data.frame(
    "rare_expansion_count" = i,
    "cohort" = "SSC",
    "samples" = ssc.i, stringsAsFactors = F
  ))
  
}

count$samples[count$cohort == "MSSNG"] <- count$samples[count$cohort == "MSSNG"]/mssng.sum
count$samples[count$cohort == "SSC"] <- count$samples[count$cohort == "SSC"]/ssc.sum

ggplot(count, aes(x = rare_expansion_count, y = samples, fill = cohort)) + geom_bar(stat = "identity", position = position_dodge()) +
  xlim(-0.5, 5) + theme_classic() + ylab("Proportion") + xlab("Rare expansion count") +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(fill = NA), legend.title = element_blank())
ggsave("ext.fig12.pdf", width = 6, height = 5)
tmp.tmp <- tmp.f

### global burden test
fam.data$dummy <- 1
fam.data$FAMILY_TYPE[fam.data$Dataset == "SSC"] <- "SPX"
cov = c("EHdnGenotypeNCLoci", "Sex")

tmp.f <- tmp.tmp[which(tmp.tmp$repeatID %in% annotation.f), ]
tmp.agr.res <- testByAggregate(unique(tmp.f$repeatID[which(tmp.f$repeatID %in% annotation.f)]), 
                               tmp.f, 
                               fam.data[fam.data$Dataset == "SSC" & fam.data$ETH_TAG == "EUR", ], 
                               ehdn.genotype.data, cov, standardization = F)
tmp.agr.res$all.aggr.glm.OR
tmp.agr.res$all.aggr.glm.pvalue

write.table(tmp.agr.res, "global.any.motif.genic.tsv", sep="\t", row.names=F, quote = F, col.names=T)

########### gene set analysis #############
cov = c("EHdnGenotypeNCLoci", "Sex")

gene.info <- annotation[!is.na(annotation$entrez_id) & 
                          annotation$typeseq_priority %in% c("exonic", "exonic;splicing", "splicing", "upstream", 
                                                             "UTR3", "UTR5", "downstream", "intronic"), ]
  
load("gs.asd.2019.RData")
gsASD <- gsASD[-c(33)]
tmp.gs <- data.frame()
type <- "genicCorrectedByNC"

for(gs in names(gsASD)){
  print(gs)
  flush.console()
  ann.tmp <- gene.info$key[gene.info$entrez_id %in% as.numeric(gsASD[[gs]])]
  
  if(length(ann.tmp) > 0){
    tmp.gs.res <- tryCatch(testByAggregate(unique(tmp.f$repeatID[tmp.f$repeatID %in% ann.tmp]), 
                                  tmp.f[tmp.f$repeatID %in% ann.tmp, ],
                                  fam.data[fam.data$Dataset == "SSC" & fam.data$ETH_TAG == "EUR", ], ehdn.genotype.data, cov, standardization = F),
                           error = function(e){
                             message(sprintf("Error %s", gs))
                           })
    
    if(!is.null(tmp.gs.res)){
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
write.table(tmp.gs, "geneset.any.motif.tsv",
            sep = "\t", row.names=F, quote=F, col.names=T)

######### Gene element test
##############################################
functions <- list()

functions[["upstream"]] <- annotation$key[annotation$typeseq_priority == "upstream" & annotation$chr %in% tmp.f$chr]
functions[["downstream"]] <- annotation$key[annotation$typeseq_priority == "downstream" & annotation$chr %in% tmp.f$chr]
functions[["intron"]] <- annotation$key[annotation$typeseq_priority == "intronic" & annotation$chr %in% tmp.f$chr]
functions[["UTR5"]] <- annotation$key[annotation$typeseq_priority == "UTR5" & annotation$chr %in% tmp.f$chr]
functions[["UTR3"]] <- annotation$key[annotation$typeseq_priority == "UTR3" & annotation$chr %in% tmp.f$chr]
functions[["exon"]] <- annotation$key[annotation$typeseq_priority %in% c("exonic") & annotation$chr %in% tmp.f$chr]
functions[["splicing"]] <- annotation$key[(annotation$typeseq_priority %in% c("splicing", "exonic;splicing")) & annotation$chr %in% tmp.f$chr]
functions[["intergenic"]] <- annotation$key[annotation$typeseq_priority %in% c("intergenic") & annotation$chr %in% tmp.f$chr]

sapply(functions, length)

cov = c("EHdnGenotypeNCLoci", "Sex")

tmp.f <- tmp.tmp
tmp.fs <- data.frame()
for(fs in names(functions)){
  print(fs)
  flush.console()
  
  tmp.fs.res <- testByAggregate(unique(tmp.f$repeatID[tmp.f$repeatID %in% functions[[fs]]]), 
                                tmp.f[tmp.f$repeatID %in% functions[[fs]], ], 
                                fam.data[fam.data$Dataset == "SSC", ], 
                                ehdn.genotype.data, cov, standardization = F)
  tmp.fs.res$geneset = fs
  
  tmp.fs <- rbind(tmp.fs, tmp.fs.res)
}

tmp.fs$aggr.all.bhfdr <- p.adjust(tmp.fs$all.aggr.glm.pvalue, method = "BH")
tmp.fs$aggr.all.fwer <- p.adjust(tmp.fs$all.aggr.glm.pvalue, method = "bonferroni")
tmp.fs

write.table(tmp.fs, "functional.any.motif.intergenic.correction.tsv",
            sep = "\t", row.names=F, quote=F, col.names=T)



##### test additional ASD-related geneset
library(GSBurden)
cov = c("EHdnGenotypeNCLoci", "Sex")

gene.info <- annotation[!is.na(annotation$entrez_id) & 
                          annotation$typeseq_priority %in% c("exonic", "exonic;splicing", "splicing", "upstream", 
                                                             "UTR3", "UTR5", "downstream", "intronic"), ]


add.gs <- read.delim("HRG_iRIGs_SCZ-DEPICT-SFARI-ASD_GWAS-ASD_102-OTHER_ASD.clean.txt", stringsAsFactors = F)
add.gs <- add.gs[, c("Entrez_gene_id", "gene_set")]

tmp.gs <- data.frame()
type <- "genicCorrectedByNC"

for(gs in names(add.gs)[!names(add.gs) %in% c("ASD_GWAS")]){
  print(gs)
  flush.console()
  ann.tmp <- gene.info$key[gene.info$entrez_id %in% as.numeric(add.gs[[gs]]) & 
                             gene.info$entrez_id %in% c(gsASD$Neurof_GoNervSysDev, gsASD$PhMm_Aggr_CardvascMuscle_all)]
  
  if(length(ann.tmp) > 0){
    tmp.gs.res <- tryCatch(testByAggregate(unique(tmp.f$repeatID[tmp.f$repeatID %in% ann.tmp]), 
                                           tmp.f[tmp.f$repeatID %in% ann.tmp, ],
                                           fam.data[fam.data$Dataset == "SSC" & fam.data$ETH_TAG == "EUR", ], ehdn.genotype.data, cov, standardization = F),
                           error = function(e){
                             message(sprintf("Error %s", gs))
                           })
    
    if(!is.null(tmp.gs.res)){
      tmp.gs.res$geneset = gs
      tmp.gs.res$type = type
      
      tmp.gs <- rbind(tmp.gs, tmp.gs.res)
    }
  }
}

tmp.gs$geneset.size <- sapply(add.gs[!names(add.gs) %in% c("ASD_GWAS")], length)
tmp.gs$aggr.all.bhfdr <- 1
tmp.gs$aggr.all.fwer <- 1
tmp.gs$aggr.all.bhfdr <- p.adjust(tmp.gs$all.aggr.glm.pvalue, method = "BH")
tmp.gs$aggr.all.fwer <- p.adjust(tmp.gs$all.aggr.glm.pvalue, method = "bonferroni")
write.table(tmp.gs, "additional.geneset.any.motif.tsv",
            sep = "\t", row.names=F, quote=F, col.names=T)

