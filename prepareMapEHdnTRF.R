# This script is to map EHdn tandem repeats to TRF repeats based on coordinate overlap and motif matching.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
queryString <- function(str1, str2){
  ori.str2 <- str2
  str2 <- paste0(str2, str2)
  
  length.str1 <- nchar(str1)
  length.str2 <- nchar(str2)
  revcom.str1 <- as.character(reverseComplement(DNAString(str1)))
  
  length.diff <- ifelse(length.str2 > length.str1, length.str2 - length.str1, 1)
  
  identity <- c()
  for(i in 1:length.diff){
    tmp <- substr(str2, i, length.str1+i-1)
    match <- sum(strsplit(tmp, "")[[1]] == strsplit(str1, "")[[1]])
    match.recom <- sum(strsplit(tmp, "")[[1]] == strsplit(revcom.str1, "")[[1]])
    
    match <- min(max(match, match.recom), nchar(ori.str2))
    tmp.identity <- match/length.str1
    
    identity <- c(identity, tmp.identity)
  }
  return(max(identity))
}


library(Biostrings)
library(data.table)
library(GenomicRanges)

ehdn <- fread("output_regions.min2.1000G+SSC+MSSNG.txt")
trf <- fread("UCSC_simple_repeats_hg38.txt")
ehdn <- as.data.frame(ehdn)[, c(1:6)]
ehdn$repeatID <- paste(ehdn$V1, ehdn$V2, ehdn$V3, ehdn$V4, sep="#")
names(ehdn) <- c("chr", "start", "end", "motif", "var1", "var2", "repeatID")

trf <- as.data.frame(trf)
trf$repeatID <- paste(trf$V1, trf$V2, trf$V3, trf$V4, sep="#")

ehdn.g <- GRanges(ehdn$chr, IRanges(ehdn$start, ehdn$end), "*")
trf.g <- GRanges(trf$V1, IRanges(trf$V2, trf$V3), "*")

olap <- data.frame(findOverlaps(trf.g, ehdn.g))
olap$trf.motif <- trf$V4[olap$queryHits]
olap$ehdn.motif <- ehdn$motif[olap$subjectHits]
olap$trf.id <- trf$repeatID[olap$queryHits]
olap$ehdn.id <- ehdn$repeatID[olap$subjectHits]

for(i in 1:nrow(olap)){
  trf.identity <- queryString(olap$trf.motif[i], olap$ehdn.motif[i])
  ehdn.identity <- queryString(olap$ehdn.motif[i], olap$trf.motif[i])
  
  olap$reciprocalIdentity[i] <- min(trf.identity, ehdn.identity)
}

olap$length.diff <- nchar(olap$trf.motif) - nchar(olap$ehdn.motif)
write.table(olap[, c(5:6, 7:8)], "map.all.tsv", sep="\t", row.names=F, quote=F, col.names=T)

olap <- read.delim("map.all.tsv", stringsAsFactors = F)
olap.f <- olap[order(olap$reciprocalIdentity, decreasing = T), ]
olap.f <- olap.f[olap.f$reciprocalIdentity > 0.66, ]
olap.f <- aggregate(trf.id ~ ehdn.id + reciprocalIdentity + length.diff, olap.f, paste, collapse = ";")
olap.f <- unique(olap.f)
olap.f <- olap.f[order(olap.f$reciprocalIdentity, decreasing = T), ]

olap.f <- olap.f[!duplicated(olap.f$ehdn.id), ]

write.table(olap.f[, c(4, 1, 2, 3)], "map.TRF.EHdn.0.66.tsv", sep="\t", row.names=F, quote=F, col.names=T)
