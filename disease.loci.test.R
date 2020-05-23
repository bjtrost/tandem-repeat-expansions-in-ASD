### Fisher's Exact test for disease loci burden

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

fam.data <- read.delim("fam.data.kids.clean.tsv", stringsAsFactors = F)
fam.data <- fam.data[fam.data$Dataset == "SSC" & fam.data$ETH_TAG == "EUR", ]
geno <- read.delim("potentially_pathogenic_STRs.txt", stringsAsFactors = F)

dt.out <- data.frame()

for(gene in unique(geno$Gene)){
  tmp <- geno[geno$Gene == gene, ]
  sex <- c("female", "male")
  if(length(grep("chrX", tmp$Locus[1])) > 0){
    sex <- c("male")
  }
  
  case.n <- length(unique(fam.data$Family.ID[which(fam.data$Status == "AffectedKid" & fam.data$Sex %in% sex)]))
  ctrl.n <- length(unique(fam.data$Family.ID[which(fam.data$Status == "UnaffectedKid" & fam.data$Sex %in% sex)]))
  
  case <- length(unique(fam.data$Family.ID[which(fam.data$Sample.ID[fam.data$Status == "AffectedKid" & fam.data$Sex %in% sex] %in% tmp$Individual)]))
  ctrl <- length(unique(fam.data$Family.ID[which(fam.data$Sample.ID[fam.data$Status == "UnaffectedKid" & fam.data$Sex %in% sex] %in% tmp$Individual)]))
  
  fisher <- fisher.test(data.frame("obs" = c(case, ctrl), "unobs" = c(case.n - case, ctrl.n - ctrl)))
  
  dt.out <- rbind(dt.out,
                  data.frame("Locus" = tmp$Locus[1],
                             "Gene" = gene,
                             "N.cases" = case,
                             "N.ctrls" = ctrl,
                             "Total.cases" = case.n,
                             "Total.ctrls" = ctrl.n,
                             "OR" = round(fisher$estimate, digits = 2),
                             "OR.upper" = round(fisher$conf.int[2], digits = 2),
                             "OR.lower" = round(fisher$conf.int[1], digits = 2),
                             "pvalue" = round(fisher$p.value, digits = 4)))
}

write.table(dt.out[order(dt.out$OR, decreasing =  T), ], "Fisher.OR.disease.loci.09032020.tsv", sep="\t", row.names=F,
            col.names = T, quote=F)
