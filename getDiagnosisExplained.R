# script to get diagnosis explained by expansions

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
annotation <- read.delim("merged.expansion.coordinates_forAnno.txt.annovar.out_rev38.1mFmt.tsv", stringsAsFactors = F)
annotation$key <- paste(annotation$chr, annotation$start, annotation$end, sep="#")
annotation$typeseq_priority <- factor(annotation$typeseq_priority, 
                                      levels = c("exonic", "splicing", "exonic;splicing", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_exonic;ncRNA_splicing", 
                                                 "UTR5", "UTR3", 
                                                 "intronic", "ncRNA_intronic", "upstream", "downstream", "intergenic"))
annotation <- annotation[order(annotation$typeseq_priority), ]
annotation <- annotation[!duplicated(annotation$key), ]
annotation.intergenic <- annotation$key[annotation$typeseq_priority == "intergenic"]
all.error <- readLines("error.with.ssc.hiseq2500.txt")

brett.samples <- read.delim("manifest.1000G+SSC+MSSNG.txt", stringsAsFactors = F, header = F)
brett.samples <- brett.samples[!brett.samples$V1 %in% all.error, ]

tmp.outlier <- read.delim("MSSNG.samples.with.excessive.no.detected.EHdn.both.sides.tsv", stringsAsFactors = F)

fam.data.tmp <- read.delim("Data_Merged_MSSNG.SSC.1KG.tsv", stringsAsFactors = F)
fam.data.tmp <- fam.data.tmp[fam.data.tmp$Sample.ID %in% brett.samples$V1, ]

fam.data.tmp <- fam.data.tmp[!fam.data.tmp$Sample.ID %in% tmp.outlier$Sample, ]
fam.1000g <- fam.data.tmp[fam.data.tmp$Dataset == "1000G", ]
fam.data <- fam.data.tmp[which(!is.na(fam.data.tmp$Sample.ID)), ]
fam.data <- fam.data[fam.data$Status != "", ]
fam.data <- fam.data[which(fam.data$ETH_TAG == "EUR" & fam.data$Dataset == "SSC"), ]

### limit to one affected kid per family
fam.data <- fam.data[which((!duplicated(paste(fam.data$Family.ID, fam.data$Status)) & fam.data$Status == "AffectedKid") |
                             fam.data$Status != "AffectedKid"), ]
fam.data <- fam.data[fam.data$DNASOURCE == "Other cell", ]
fam.data <- fam.data[fam.data$Library.type == "PCR-free" & fam.data$Platform == "Illumina HiSeqX", ]
fam.data <- fam.data[fam.data$Status %in% c("AffectedKid", "UnaffectedKid"), ]

ehdn.result.m <- read.delim("merged.any.outlier.expansions.autosomal.sex.chr.nobadsamples.tsv", stringsAsFactors = F)
ehdn.genotype.data <- strsplit(paste(ehdn.result.m$outliers[ehdn.result.m$chr %in% paste0("chr", 1:22) & 
                                                              ehdn.result.m$repeatID %in% annotation.intergenic], collapse = ";"), ";")[[1]]
ehdn.genotype.data <- data.frame(table(ehdn.genotype.data),
                                 stringsAsFactors = F)
names(ehdn.genotype.data) <- c("Sample", "EHdnGenotypeNCLoci")

for(type in c("all rare", "genic rare", "geneset rare")){
  dt <- read.delim("Merged.Any.Expansion.with.outliers.EHdn.ssc.eur.tsv", stringsAsFactors = F)
  dt <- dt[dt$chr.x %in% paste0("chr", c(1:22)), ]
  if(type == "genic rare")
    dt <- dt[!is.na(dt$gene_symbol), ]
  if(type == "geneset rare"){
    dt <- dt[dt$insigset, ]
    sup.table <- read.delim("supplementary.table.with.cutoff.tsv", stringsAsFactors = F)
    dt <- dt[dt$repeatID %in% sup.table$repeatID, ]
  }
  samples <- strsplit(paste(dt$samples, collapse = ";"), ";")[[1]]
  samples <- data.frame(table(samples))
  
  dt.lm <- merge(fam.data, samples, by.x = "Sample.ID", by.y = "samples", all.x = T)
  
  dt.lm$Freq[is.na(dt.lm$Freq)] <- 0
  affected.rates <- sum(dt.lm$Status == "AffectedKid" & dt.lm$Freq > 0, na.rm = T)/
    sum(dt.lm$Status == "AffectedKid", na.rm = T)
  
  unaffected.rates <- sum(dt.lm$Status == "UnaffectedKid" & dt.lm$Freq > 0, na.rm = T)/
    sum(dt.lm$Status == "UnaffectedKid", na.rm = T)
  explained <- (affected.rates - unaffected.rates)
  test <- wilcox.test(dt.lm$Freq[dt.lm$Status == "AffectedKid"],dt.lm$Freq[dt.lm$Status == "UnaffectedKid"], alternative = "greater")
  
  print(sprintf("%s, affected = %s, affected count = %s, unaffected = %s, contribution = %s, pvalue=%s", 
                type, sum(dt.lm$Status == "AffectedKid" & dt.lm$Freq > 0), round(affected.rates, digits = 3), round(unaffected.rates, digits = 3), 
                round(explained, digits = 4), test$p.value))
}


### EUR SSC april update
# [1] "all rare, affected = 418, affected count = 0.231, unaffected = 0.207, contribution = 0.0233, pvalue=0.0501155516232843"
# [1] "genic rare, affected = 249, affected count = 0.137, unaffected = 0.115, contribution = 0.0223, pvalue=0.0275504111660168"
# [1] "geneset rare, affected = 47, affected count = 0.026, unaffected = 0.007, contribution = 0.0192, pvalue=1.28612423235789e-05"



dt <- read.delim("Merged.Any.Expansion.with.outliers.EHdn.ssc.eur.tsv", stringsAsFactors = F)
dt <- dt[dt$chr.x %in% paste0("chr", 1:22), ]
dt <- dt[dt$insigset, ]
dt <- dt[dt$repeatID %in% sup.table$repeatID, ]

dt <- dt[!is.na(dt$gene_symbol), ]
dt <- dt[!duplicated(dt$repeatID), ]
samples <- strsplit(paste(dt$samples, collapse = ";"), ";")[[1]]
samples <- data.frame(table(samples))

dt.lm <- merge(fam.data, samples, by.x = "Sample.ID", by.y = "samples", all.x = T)
dt.lm <- merge(dt.lm, ehdn.genotype.data, by.x = "Sample.ID", by.y = "Sample", all.x = T)

dt.lm$Freq[is.na(dt.lm$Freq)] <- 0
dt.lm$EHdnGenotypeNCLoci[is.na(dt.lm$EHdnGenotypeNCLoci)] <- 0

dt.lm$Status <- factor(dt.lm$Status, levels = c("UnaffectedKid", "AffectedKid"))
dt.lm <- dt.lm[dt.lm$Dataset == "SSC", ]
sum(dt.lm$Status == "AffectedKid" & dt.lm$Freq > 0)/
  sum(dt.lm$Status == "AffectedKid")
# [1] 0.02593819 

fisher.test(data.frame(c(sum(dt.lm$Status == "AffectedKid" & dt.lm$Freq > 0), 
                       sum(dt.lm$Status == "UnaffectedKid" & dt.lm$Freq > 0)),
            c(sum(dt.lm$Status == "AffectedKid" & dt.lm$Freq == 0),
                       sum(dt.lm$Status == "UnaffectedKid" & dt.lm$Freq == 0))))
 
#   p-value = 1.945e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.949417 8.746200
# sample estimates:
#   odds ratio 
# 3.926296 
