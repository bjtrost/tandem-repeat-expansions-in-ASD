### library required to perform statistical analysis in other scripts.

get1000Freq <- function(outlier, fam.1000g){
  outlier <- strsplit(outlier, ";")[[1]]
  onek.count <- sum(unique(outlier) %in% fam.1000g$Sample.ID)
  return(onek.count/nrow(fam.1000g))
}

getOutlierData <- function(i, ehdn.result, ehdn, trf, fam.data, fam.1000g){
  ehdn.rec <- ehdn.result[i, ]
  ehdn.data <- ehdn[ehdn$repeatID == ehdn.rec$repeatID, ]
  
  outliers <- strsplit(ehdn.rec$outliers, ";")[[1]]
  
  ehdn.data <- unlist(strsplit(gsub(":|,", "#", ehdn.data[, 7]), "#")[[1]]) #V7 is genotype data from EHDN
  ehdn.data <- data.frame(matrix(ehdn.data, ncol = 2, byrow = T), stringsAsFactors = F)
  ehdn.data[, 2] <- as.numeric(ehdn.data[, 2])
  names(ehdn.data) <- c("Sample", "EHDN")
  
  ehdn.data$ehdn.genotypeable <- sum(!is.na(ehdn.data$EHDN))
  ehdn.data$ehdn.mode <- names(table(ehdn.data$EHDN)[which.max(table(ehdn.data$EHDN))])
  trf.data <- trf[, c(1,2, which(names(trf) == ehdn.rec$trf.id))]
  trf.data$longAllele <- NA
  trf.mode <- NA
  rsquared <- NA
  if(!is.na(ehdn.rec$trf.id)){
    trf.data$longAllele <- sapply(lapply(sapply(trf.data[, 3], strsplit, "/"), as.numeric), max, na.rm = T)
    ehdn.tmp <- merge(ehdn.data, trf.data, by = "Sample", all = T)
    ehdn.no.na <- na.omit(ehdn.tmp)
    rsquared <- cor(ehdn.no.na$EHDN, ehdn.no.na$longAllele)
    trf.mode <- names(table(trf.data$longAllele)[which.max(table(trf.data$longAllele))])
  }
  
  ehdn.data$trf.genotypeable <- sum(!is.na(trf.data$longAllele))
  ehdn.data$both.genotypeable <- sum(ehdn.data$Sample %in% trf.data$Sample[!is.na(trf.data$longAllele)])
  ehdn.data$trf.mode <- ifelse(!is.na("trf.mode"), trf.mode, NA)
  ehdn.data$rsquared <- ifelse(!is.na("rsquared"), rsquared, NA)
  
  ehdn.data <- merge(ehdn.data, trf.data[, c("Sample", "longAllele")], by = "Sample", all = T)
  
  ehdn.data <- ehdn.data[ehdn.data$Sample %in% outliers, ]
  
  ehdn.data <- merge(ehdn.data, fam.data[, c("Sample.ID", "Family.ID", "Sex", "Relation", "Affection", "Dataset", "Status", "DNASOURCE",
                                             "Platform", "FAMILY_TYPE")], by.x = "Sample", 
                     by.y = "Sample.ID", all.x = T)
  
  tmp <- do.call(rbind, lapply(1:nrow(ehdn.data), getInheritance, ehdn.data, fam.data))
  ehdn.data[, names(tmp)] <- tmp
  ehdn.data$MSSNG <- sum(ehdn.data$Dataset == "MSSNG", na.rm = T)
  ehdn.data$SSC <- sum(ehdn.data$Dataset == "SSC", na.rm = T)
  ehdn.data$REF1000G <- sum(outliers %in% fam.1000g$Sample.ID, na.rm = T)

  ehdn.data[, c("repeatID", "motif", "chr", "start", "end", "trf.id")] <- ehdn.rec[, c("repeatID", "motif", "chr", "start", "end", "trf.id")]
  return(ehdn.data)
}

getOutlierMergeData <- function(i, ehdn.result.m, fam.data, fam.1000g){
  ehdn.rec <- ehdn.result.m[i, ]
  outliers <- strsplit(ehdn.rec$outliers, ";")[[1]]
  
  
  ehdn.data <-  fam.data[fam.data$Sample.ID %in% outliers, c("Sample.ID", "Family.ID", "Sex", "Relation", "Affection", "Dataset", "Status", "DNASOURCE",
                                                             "Platform", "FAMILY_TYPE")]
  if(nrow(ehdn.data) == 0){
    return(NULL)
  }else{
    tmp <- do.call(rbind, lapply(1:nrow(ehdn.data), getInheritance, ehdn.data, fam.data))
    ehdn.data[, names(tmp)] <- tmp
    ehdn.data$MSSNG <- sum(ehdn.data$Dataset == "MSSNG", na.rm = T)
    ehdn.data$SSC <- sum(ehdn.data$Dataset == "SSC", na.rm = T)
    ehdn.data$REF1000G <- sum(outliers %in% fam.1000g$Sample.ID, na.rm = T)

    ehdn.data[, c("repeatID", "motif", "chr", "start", "end", "trf.id")] <- ehdn.rec[, c("repeatID", "motif", "chr", "start", "end", "trf.id")]
  }
  return(ehdn.data)
}

mapOutlierToMergeLoci <- function(outlier.loci, merge.loci){
  tmp <- unique(as.numeric(unlist(sapply(outlier.loci$repeatID, grep, merge.loci$repeatID))))
  merge.loci <- merge.loci[tmp, ]
  
  for(i in 1:nrow(merge.loci)){
    repeatIDs <- strsplit(merge.loci$repeatID[i], ";")[[1]]
    outliers <- paste(unique(unlist(sapply(outlier.loci$outliers[outlier.loci$repeatID %in% repeatIDs], strsplit, ";"))), collapse = ";")
    merge.loci$outliers[i] <- outliers
  }
  
  return(merge.loci)
}

mergeLoci <- function(loci, mergeAny = T){
  
  cnv.tmp <- loci
  
  firstIt <- T
  olap <- data.frame()
  while(nrow(olap) != 0 | firstIt){
    firstIt <- F
    cnv.query <- GRanges(cnv.tmp$chr, IRanges(cnv.tmp$start, cnv.tmp$end), "*")
    
    olap <- data.frame(findOverlaps(cnv.query, cnv.query))
    olap <- olap[which(olap$queryHits < olap$subjectHits), ]
    olap <- olap[!duplicated(olap$queryHits), ]
    olap <- olap[!olap$queryHits %in% olap$subjectHits, ]
    if(!mergeAny){
      olap <- olap[nchar(cnv.tmp$motif[olap$queryHits]) == nchar(cnv.tmp$motif[olap$subjectHits]), ]
    }
    
    message(sprintf("%s found overlapped, from %s loci", nrow(olap), nrow(cnv.tmp)))
    flush.console()
    
    if(nrow(olap) > 0){
      ids <- union(olap$queryHits, olap$subjectHits)
      merge.dt <- data.frame()
      for(i in 1:nrow(olap)){
        olap.th <- olap[i, ]
        
        olap.th[, c("qstart", "qend")] <- cnv.tmp[olap.th$queryHits, c("start", "end")]
        olap.th[, c("sstart", "send")] <- cnv.tmp[olap.th$subjectHits, c("start", "end")]
        
        
        repeatID <- paste(union(strsplit(cnv.tmp$repeatID[olap.th$queryHits], ";")[[1]], 
                                 strsplit(cnv.tmp$repeatID[olap.th$subjectHits], ";")[[1]]), collapse=";")
        motif <- paste(union(strsplit(cnv.tmp$motif[olap.th$queryHits], ";")[[1]], 
                              strsplit(cnv.tmp$motif[olap.th$subjectHits], ";")[[1]]), collapse = ";")
        trf.id <- paste(na.omit(union(strsplit(cnv.tmp$trf.id[olap.th$queryHits], ";")[[1]], 
                               strsplit(cnv.tmp$trf.id[olap.th$subjectHits], ";")[[1]])), collapse = ";")
        reciprocalIdentity <- max(cnv.tmp$reciprocalIdentity[olap.th$queryHits], cnv.tmp$reciprocalIdentity[olap.th$subjectHits], na.rm =T)
        length.diff <- 0
        
        merge.dt <- rbind(merge.dt, data.frame(repeatID, "chr" = cnv.tmp$chr[olap.th$queryHits],
                                      "start" = pmin(olap.th$qstart, olap.th$sstart), "end" = pmax(olap.th$qend, olap.th$send),
                                      motif, trf.id, reciprocalIdentity, length.diff))
      }
    
      
      cnv.tmp <- cnv.tmp[-(ids), ]
      cnv.tmp <- rbind(cnv.tmp, merge.dt)
    }
    
    gc();
  }
  
  return(cnv.tmp)
}


getInheritance <- function(i, dt, fam.data){
  rec <- dt[i, ]
  denovo <- F
  transmitted.case.count <- NA
  nontransmitted.case.count <- NA
  transmitted.ctrl.count <- NA
  nontransmitted.ctrl.count <- NA
  if(rec$Status %in% c("AffectedKid", "UnaffectedKid")){
    if(sum(dt$Family.ID == rec$Family.ID & dt$Status %in% c("AffectedParent", "UnaffectedParent"), na.rm = T) == 0)
      denovo <- T
  }else{
    transmitted.case.count <- sum(dt$Family.ID == rec$Family.ID & dt$Status == "AffectedKid", na.rm = T)
    nontransmitted.case.count <- sum(fam.data$Family.ID == rec$Family.ID &
                                       fam.data$Status == "AffectedKid", na.rm = T) - transmitted.case.count
    
    transmitted.ctrl.count <- sum(dt$Family.ID == rec$Family.ID & dt$Status == "UnaffectedKid", na.rm = T)
    nontransmitted.ctrl.count <- sum(fam.data$Family.ID == rec$Family.ID &
                                       fam.data$Status == "UnaffectedKid", na.rm = T) - transmitted.ctrl.count
  }
  
  return(data.frame(denovo, transmitted.case.count, nontransmitted.case.count, transmitted.ctrl.count, nontransmitted.ctrl.count, stringsAsFactors = F))
}

testTransmission <- function(transmitted.dt, label, label.list, test = "chisq"){
  dt.out <- data.frame("OR" = 1, "Pvalue" = 1)
  dt.out[, label.list] <- 0
  
  if(label == "Affection"){
    case.numerator <- sum(transmitted.dt$transmitted.case.count)
    ctrl.numerator <- sum(transmitted.dt$transmitted.ctrl.count)
    case.denominator <- sum(transmitted.dt$nontransmitted.case.count) + case.numerator
    ctrl.denominator <- sum(transmitted.dt$nontransmitted.ctrl.count) + ctrl.numerator
  }else if(label == "FAMILY_TYPE"){
    case.numerator <- sum(transmitted.dt$transmitted.case.count[transmitted.dt[, label] == label.list[1]])
    ctrl.numerator <- sum(transmitted.dt$transmitted.case.count[transmitted.dt[, label] == label.list[2]])
    case.denominator <- sum(transmitted.dt$nontransmitted.case.count[transmitted.dt[, label] == label.list[1]]) + case.numerator
    ctrl.denominator <- sum(transmitted.dt$nontransmitted.case.count[transmitted.dt[, label] == label.list[2]]) + ctrl.numerator
  }
  
  dt.out[, paste0(label.list[1], ".denominator")] <- case.denominator
  dt.out[, paste0(label.list[2], ".denominator")] <- ctrl.denominator
  dt.out[, paste0(label.list[1], ".transmission.rate")] <- 0
  dt.out[, paste0(label.list[2], ".transmission.rate")] <- 0

  if(nrow(transmitted.dt) > 0){
    
    if(test == "chisq"){
      ftest <- chisq.test(data.frame(c(case.numerator, ctrl.numerator), c(case.denominator-case.numerator, ctrl.denominator-ctrl.numerator)))
      dt.out$OR = ftest$statistic
    }else{
      ftest <- fisher.test(data.frame(c(case.numerator, ctrl.numerator), c(case.denominator-case.numerator, ctrl.denominator-ctrl.numerator)))
      dt.out$OR = ftest$estimate
    }
    dt.out$Pvalue = ftest$p.value
    dt.out[, label.list[1]] <- case.numerator
    dt.out[, label.list[2]] <- ctrl.numerator
    dt.out[, paste0(label.list[1], ".transmission.rate")] <- round(case.numerator/case.denominator, digits = 3)
    dt.out[, paste0(label.list[2], ".transmission.rate")] <- round(ctrl.numerator/ctrl.denominator, digits = 3)
  }
  
  return(dt.out)
}

testDenovo <- function(denovo.dt, label, label.list, case.denominator, ctrl.denominator, test = "chisq"){
  dt.out <- data.frame("OR" = 1, "Pvalue" = 1)
  dt.out[, label.list] <- 0
  dt.out[, paste0(label.list[1], ".denominator")] <- case.denominator
  dt.out[, paste0(label.list[2], ".denominator")] <- ctrl.denominator
  
  if(nrow(denovo.dt) > 0){
    case <- length(unique(denovo.dt$Family.ID[denovo.dt[, label] == label.list[1]]))
    ctrl <- length(unique(denovo.dt$Family.ID[denovo.dt[, label] == label.list[2]]))
    if(test == "chisq"){
      ftest <- chisq.test(data.frame(c(case, ctrl), c(case.denominator-case, ctrl.denominator-ctrl)))
      dt.out$OR = ftest$statistic
    }else{
      ftest <- fisher.test(data.frame(c(case, ctrl), c(case.denominator-case, ctrl.denominator-ctrl)))
      dt.out$OR = ftest$estimate
    }
    dt.out$Pvalue = ftest$p.value
    dt.out[, label.list[1]] <- case
    dt.out[, label.list[2]] <- ctrl
  }
  
  return(dt.out)
}

testCaseCtrl <- function(caseCtrl.dt, label, label.list, case.denominator, ctrl.denominator, test = "chisq", fam.data){
  dt.out <- data.frame("OR" = 1, "Pvalue" = 1)
  dt.out[, label.list] <- 0
    
  dt.out[, paste0(label.list[1], ".denominator")] <- case.denominator
  dt.out[, paste0(label.list[2], ".denominator")] <- ctrl.denominator
  
  if(nrow(caseCtrl.dt) > 0){
    case <- length(unique(caseCtrl.dt$Family.ID[caseCtrl.dt[, label] == label.list[1] & caseCtrl.dt$Sample %in% fam.data$Sample.ID]))
    ctrl <- length(unique(caseCtrl.dt$Family.ID[caseCtrl.dt[, label] == label.list[2] & caseCtrl.dt$Sample %in% fam.data$Sample.ID]))
    # if(test == "chisq"){
    #   ftest <- chisq.test(data.frame(c(case, ctrl), c(case.denominator-case, ctrl.denominator-ctrl)))
    #   dt.out$OR = ftest$statistic
    # }else{
    #   ftest <- fisher.test(data.frame(c(case, ctrl), c(case.denominator-case, ctrl.denominator-ctrl)))
    #   dt.out$OR = ftest$estimate
    # }
    # 
    # dt.out$Pvalue = ftest$p.value
    dt.out[, label.list[1]] <- case
    dt.out[, label.list[2]] <- ctrl
  }
  return(dt.out[, -c(1,2)])
}

getCountTable <- function(dt, label, label.list, fam.data, ehdn.genotype.data, cov){
  if(nrow(dt) == 0){
    return(NA)
  }else{
    fam.tmp <- fam.data[which(fam.data$Dataset %in% unique(dt$Dataset) & 
                                fam.data[, label] %in% label.list & fam.data$Status %in% dt$Status), ]
    dt <- dt[dt$Sample %in% fam.tmp$Sample.ID, ]
    
    if(sum(names(fam.tmp) %in% cov) > 0)
      cov.dt <- unique(fam.tmp[, c("Sample.ID", names(fam.tmp)[names(fam.tmp) %in% cov])])

    table.tmp <- aggregate(repeatID ~ Sample, dt, length)
    names(table.tmp)[2] <- "ExpansionCount" 
    
    table.tmp <- merge(table.tmp, fam.tmp[, c("Sample.ID", label)], by.x = "Sample", by.y = "Sample.ID", all = T)
    table.tmp[is.na(table.tmp)] <- 0
    table.tmp[, label] <- ifelse(table.tmp[, label] == label.list[1], 1, 0)
    table.tmp[, label] <- factor(table.tmp[, label])
    if(sum(names(fam.tmp) %in% cov) > 0)
      table.tmp <- merge(table.tmp, cov.dt, by.x = "Sample", by.y = "Sample.ID", all = T)
    table.tmp <- merge(table.tmp, ehdn.genotype.data, by = "Sample", all.x = T)
    table.tmp[is.na(table.tmp)] <- 0
    return(table.tmp)
  }
}

testGLM <- function(dt, label, cov, feature, standardization = T){
  if(is.null(nrow(dt))){
    return(list("OR" = NA, "upper" = NA, "lower" = NA,  "pvalue" = NA))
  }else if (length(unique(na.omit(dt[, label]))) < 2){
    return(list("OR" = NA, "upper" = NA, "lower" = NA,  "pvalue" = NA))
  }else{
    for(cov.each in cov){
      if(length(unique(dt[, cov.each])) < 2){
        cov <- cov[-which(cov == cov.each)]
      }
    }
    
    if(standardization){
      dt[, feature] <- scale(dt[, feature])
    }
    
    message(sprintf("Testing...%s samples (%s cases, %s ctrl)", nrow(dt), sum(dt[, label] == 1), sum(dt[, label] == 0)))
    
    if(length(cov) == 0)
      ref <- paste0(label, " ~ 1")
    else
      ref <- paste0(label, " ~ ", paste(cov, collapse = " + "))
    add <- paste0(ref, " + ", feature)
    ref.lm <- glm(ref, dt, family = binomial(link = "logit"))
    add.lm <- glm(add, dt, family = binomial(link = "logit"))
    ano <- anova(ref.lm, add.lm, test = "Chisq")
    
    or <- exp(add.lm$coefficients[feature])
    pvalue <- ano$`Pr(>Chi)`[2]
    
    conf <- confint(add.lm)
    low <- exp(conf[feature, 1])
    up <- exp(conf[feature, 2])
    
    return(list("OR" = or, "upper" = up, "lower" = low, "pvalue" = pvalue))
  }
}

testByLocus <- function(repeatID, tmp.f, fam.data){
  # mssng.fam.spx.denom <- length(unique(fam.data$Family.ID[which(fam.data$FAMILY_TYPE %in% c("SPX") & 
  #                                                                 fam.data$Status == "AffectedKid")]))
  # mssng.fam.mpx.denom <- length(unique(fam.data$Family.ID[which(fam.data$FAMILY_TYPE %in% c("MPX") & 
  #                                                                 fam.data$Status == "AffectedKid")]))
  # ssc.fam.case.denom <- length(unique(fam.data$Family.ID[which(fam.data$Dataset == "SSC" & 
  #                                                                fam.data$Status == "AffectedKid")]))
  # ssc.fam.ctrl.denom <- length(unique(fam.data$Family.ID[which(fam.data$Dataset == "SSC" & 
  #                                                                fam.data$Status == "UnaffectedKid")]))
  # 
  # mssng.spx.denom <- length(unique(fam.data$Sample.ID[which(fam.data$FAMILY_TYPE %in% c("SPX") & 
  #                                                             fam.data$Status == "AffectedKid")]))
  # mssng.mpx.denom <- length(unique(fam.data$Sample.ID[which(fam.data$FAMILY_TYPE %in% c("MPX") & 
  #                                                             fam.data$Status == "AffectedKid")]))
  # ssc.case.denom <- length(unique(fam.data$Sample.ID[which(fam.data$Dataset == "SSC" & 
  #                                                            fam.data$Status == "AffectedKid")]))
  # ssc.ctrl.denom <- length(unique(fam.data$Sample.ID[which(fam.data$Dataset == "SSC" & 
  #                                                            fam.data$Status == "UnaffectedKid")]))
  case.denom <- length(unique(fam.data$Family.ID[which((fam.data$FAMILY_TYPE %in% c("MPX", "SPX") | fam.data$Dataset == "SSC") & 
                                                         fam.data$Status == "AffectedKid")]))
  ctrl.denom <- length(unique(fam.data$Family.ID[which((fam.data$FAMILY_TYPE %in% c("MPX", "SPX") | fam.data$Dataset == "SSC") & 
                                                         fam.data$Status == "UnaffectedKid")]))
  
  
  caseCtrl.dt <- tmp.f[which(tmp.f$repeatID == repeatID & 
                               tmp.f$Status %in% c("AffectedKid", "UnaffectedKid")), ]
  
  # denovo.mssng.dt <- tmp.f[which(tmp.f$repeatID == repeatID & tmp.f$FAMILY_TYPE %in% c("SPX", "MPX") & 
  #                                  tmp.f$Status %in% c("AffectedKid", "UnaffectedKid") & tmp.f$denovo == T), ]
  # denovo.ssc.dt <- tmp.f[which(tmp.f$repeatID == repeatID & tmp.f$Dataset == "SSC" & 
  #                                tmp.f$Status %in% c("AffectedKid", "UnaffectedKid") & tmp.f$denovo == T), ]
  # denovo.all.dt <- tmp.f[which(tmp.f$repeatID == repeatID & (tmp.f$FAMILY_TYPE %in% c("MPX", "SPX") | tmp.f$Dataset == "SSC") & 
  #                                tmp.f$Status %in% c("AffectedKid", "UnaffectedKid") & tmp.f$denovo == T), ]
  
  
  # inherit.mssng.dt <- tmp.f[which(tmp.f$repeatID == repeatID & tmp.f$FAMILY_TYPE %in% c("SPX", "MPX") & 
  #                                  tmp.f$Status %in% c("AffectedKid", "UnaffectedKid") & tmp.f$denovo == F), ]
  # inherit.ssc.dt <- tmp.f[which(tmp.f$repeatID == repeatID & tmp.f$Dataset == "SSC" & 
  #                                tmp.f$Status %in% c("AffectedKid", "UnaffectedKid") & tmp.f$denovo == F), ]
  # inherit.all.dt <- tmp.f[which(tmp.f$repeatID == repeatID &(tmp.f$FAMILY_TYPE %in% c("MPX", "SPX") | tmp.f$Dataset == "SSC") & 
  #                                 tmp.f$Status %in% c("AffectedKid", "UnaffectedKid") & tmp.f$denovo == F), ]
  # 
  # transmitted.mssng.dt <- tmp.f[which(tmp.f$repeatID == repeatID & tmp.f$FAMILY_TYPE %in% c("SPX", "MPX") & 
  #                                       tmp.f$Status %in% c("AffectedParent", "UnaffectedParent")), ]
  # transmitted.ssc.dt <- tmp.f[which(tmp.f$repeatID == repeatID & tmp.f$Dataset == "SSC" & 
  #                                     tmp.f$Status %in% c("AffectedParent", "UnaffectedParent")), ]
  
  all.result <- testCaseCtrl(caseCtrl.dt, "Status", c("AffectedKid", "UnaffectedKid"), case.denom, ctrl.denom, test = "fisher", fam.data)
  # all.mssng.result <- testCaseCtrl(caseCtrl.dt[caseCtrl.dt$Dataset == "MSSNG", ], "FAMILY_TYPE", c("MPX", "SPX"), mssng.fam.spx.denom, mssng.fam.mpx.denom, test = "fisher")
  # all.ssc.result <- testCaseCtrl(caseCtrl.dt[caseCtrl.dt$Dataset == "SSC", ], "Affection", c("2", "1"), ssc.fam.case.denom, ssc.fam.ctrl.denom, test = "fisher")
  
  # names(all.mssng.result) <- paste0("MSSNG.", names(all.mssng.result))
  # names(all.ssc.result) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("SSC.", names(all.ssc.result))))
  names(all.result) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("AGGR.", names(all.result))))
  
  # mssng.denovo <- testDenovo(denovo.mssng.dt, "FAMILY_TYPE", c("MPX", "SPX"), mssng.fam.spx.denom, mssng.fam.mpx.denom, test = "fisher")
  # ssc.denovo <- testDenovo(denovo.ssc.dt, "Affection", c("2", "1"), ssc.fam.case.denom, ssc.fam.ctrl.denom, test = "fisher")
  # all.denovo <- testDenovo(denovo.all.dt, "Affection", c("2", "1"), case.denom, ctrl.denom, test = "fisher")
  
  # mssng.inherit <- testDenovo(inherit.mssng.dt, "FAMILY_TYPE", c("MPX", "SPX"), mssng.fam.spx.denom, mssng.fam.mpx.denom, test = "fisher")
  # ssc.inherit <- testDenovo(inherit.ssc.dt, "Affection", c("2", "1"), ssc.fam.case.denom, ssc.fam.ctrl.denom, test = "fisher")
  # all.inherit <- testDenovo(inherit.all.dt, "Affection", c("2", "1"), case.denom, ctrl.denom, test = "fisher")
  
  # mssng.transmission <- testTransmission(transmitted.mssng.dt, "FAMILY_TYPE", c("MPX", "SPX"), test = "fisher")
  # ssc.transmission <- testTransmission(transmitted.ssc.dt, "Affection", c("2", "1"), test = "fisher")
  
  # names(mssng.denovo) <- paste0("MSSNG.DeNovo.", names(mssng.denovo))
  # names(ssc.denovo) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("SSC.DeNovo.", names(ssc.denovo))))
  # names(all.denovo) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("AGGR.DeNovo.", names(all.denovo))))

  # names(mssng.inherit) <- paste0("MSSNG.Inherit.", names(mssng.inherit))
  # names(ssc.inherit) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("SSC.Inherit.", names(ssc.inherit))))
  # names(all.inherit) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("AGGR.Inherit.", names(all.inherit))))
  
  # names(mssng.transmission) <- paste0("MSSNG.Transmitted.", names(mssng.transmission))
  # names(ssc.transmission) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("SSC.Transmitted.", names(ssc.transmission))))
  
  samples <- paste(tmp.f$Sample[which(tmp.f$repeatID == repeatID)], collapse=";")
  ehdn.estimates <- paste(tmp.f$EHDN[which(tmp.f$repeatID == repeatID)], collapse=";")
  trf.estimates <- paste(tmp.f$longAllele[which(tmp.f$repeatID == repeatID)], collapse=";")
  status <- paste(tmp.f$Status[which(tmp.f$repeatID == repeatID)], collapse=";")
  relation <- paste(tmp.f$Relation[which(tmp.f$repeatID == repeatID)], collapse=";")
  
  prefer.columns <- c("repeatID", "motif", "chr", "start", "end", "trf.id", "MSSNG", "SSC", "REF1000G", "ehdn.genotypeable",
    "ehdn.mode",	"trf.genotypeable",	"both.genotypeable",	"trf.mode",	"rsquared")
  prefer.columns <- prefer.columns[prefer.columns %in% names(tmp.f)]
  return(data.frame(tmp.f[which(tmp.f$repeatID == repeatID), ][1, prefer.columns], samples, ehdn.estimates, trf.estimates, status, relation,
                    all.result))
  # return(data.frame(tmp.f[which(tmp.f$repeatID == repeatID), ][1, prefer.columns], samples, ehdn.estimates, trf.estimates, status, relation,
  #                   all.mssng.result, all.ssc.result, all.result, mssng.denovo, ssc.denovo, all.denovo, mssng.inherit, ssc.inherit, all.inherit, mssng.transmission, ssc.transmission))
}

testByAggregate <- function(repeatIDs, tmp.f, fam.data, ehdn.genotype.data, cov, standardization = T){
  
  # mssng.fam.spx.denom <- length(unique(fam.data$Family.ID[which(fam.data$FAMILY_TYPE %in% c("SPX") & 
  #                                                                 fam.data$Status == "AffectedKid")]))
  # mssng.fam.mpx.denom <- length(unique(fam.data$Family.ID[which(fam.data$FAMILY_TYPE %in% c("MPX") & 
  #                                                                 fam.data$Status == "AffectedKid")]))
  # ssc.fam.case.denom <- length(unique(fam.data$Family.ID[which(fam.data$Dataset == "SSC" & 
  #                                                                fam.data$Status == "AffectedKid")]))
  # ssc.fam.ctrl.denom <- length(unique(fam.data$Family.ID[which(fam.data$Dataset == "SSC" & 
  #                                                                fam.data$Status == "UnaffectedKid")]))
  # 
  # mssng.spx.denom <- length(unique(fam.data$Sample.ID[which(fam.data$FAMILY_TYPE %in% c("SPX") & 
  #                                                             fam.data$Status == "AffectedKid")]))
  # mssng.mpx.denom <- length(unique(fam.data$Sample.ID[which(fam.data$FAMILY_TYPE %in% c("MPX") & 
  #                                                             fam.data$Status == "AffectedKid")]))
  # ssc.case.denom <- length(unique(fam.data$Sample.ID[which(fam.data$Dataset == "SSC" & 
  #                                                            fam.data$Status == "AffectedKid")]))
  # ssc.ctrl.denom <- length(unique(fam.data$Sample.ID[which(fam.data$Dataset == "SSC" & 
  #                                                            fam.data$Status == "UnaffectedKid")]))
  case.denom <- length(unique(fam.data$Family.ID[which((fam.data$FAMILY_TYPE %in% c("MPX", "SPX") | fam.data$Dataset == "SSC") & 
                                                         fam.data$Status == "AffectedKid")]))
  ctrl.denom <- length(unique(fam.data$Family.ID[which((fam.data$FAMILY_TYPE %in% c("MPX", "SPX") | fam.data$Dataset == "SSC") & 
                                                         fam.data$Status == "UnaffectedKid")]))
  
  caseCtrl.dt <- tmp.f[which(tmp.f$repeatID %in% repeatIDs & 
                               tmp.f$Status %in% c("AffectedKid", "UnaffectedKid")), ]
  
  # denovo.mssng.dt <- tmp.f[which(tmp.f$repeatID %in% repeatIDs & tmp.f$FAMILY_TYPE %in% c("MPX", "SPX") & 
  #                                  tmp.f$Status %in% c("AffectedKid") & tmp.f$denovo == T), ]
  # denovo.ssc.dt <- tmp.f[which(tmp.f$repeatID %in% repeatIDs & tmp.f$Dataset == "SSC" & 
  #                                tmp.f$Status %in% c("AffectedKid", "UnaffectedKid") & tmp.f$denovo == T), ]
  # denovo.all.dt <- tmp.f[which((tmp.f$FAMILY_TYPE %in% c("MPX", "SPX") | tmp.f$Dataset == "SSC") & 
  #                                tmp.f$Status %in% c("AffectedKid", "UnaffectedKid") & tmp.f$denovo == T), ]
  
  # inherit.mssng.dt <- tmp.f[which(tmp.f$repeatID %in% repeatIDs & tmp.f$FAMILY_TYPE %in% c("MPX", "SPX") & 
  #                                   tmp.f$Status %in% c("AffectedKid", "UnaffectedKid") & tmp.f$denovo == F), ]
  # inherit.ssc.dt <- tmp.f[which(tmp.f$repeatID %in% repeatIDs & tmp.f$Dataset == "SSC" & 
  #                                 tmp.f$Status %in% c("AffectedKid", "UnaffectedKid") & tmp.f$denovo == F), ]
  # inherit.all.dt <- tmp.f[which((tmp.f$FAMILY_TYPE %in% c("MPX", "SPX") | tmp.f$Dataset == "SSC") & 
  #                                 tmp.f$Status %in% c("AffectedKid", "UnaffectedKid") & tmp.f$denovo == F), ]

  # denovo.mssng.table <- getCountTable(denovo.mssng.dt, "FAMILY_TYPE", c("MPX", "SPX"), fam.data, ehdn.genotype.data, cov)
  # denovo.ssc.table <- getCountTable(denovo.ssc.dt, "Affection", c("2", "1"), fam.data, ehdn.genotype.data, cov)
  # denovo.all.table <- getCountTable(denovo.all.dt, "Affection", c("2", "1"), fam.data, ehdn.genotype.data, cov)
  
  feature <- "ExpansionCount"
  
  # sm.mssng.denovo <- testGLM(denovo.mssng.table, "FAMILY_TYPE", cov, feature)
  # sm.ssc.denovo <- testGLM(denovo.ssc.table, "Affection", cov, feature)
  # sm.all.denovo <- testGLM(denovo.all.table, "Affection", cov, feature)
  # 
  # inherit.mssng.table <- getCountTable(inherit.mssng.dt, "FAMILY_TYPE", c("MPX", "SPX"), fam.data, ehdn.genotype.data, cov)
  # inherit.ssc.table <- getCountTable(inherit.ssc.dt, "Affection", c("2", "1"), fam.data, ehdn.genotype.data, cov)
  # inherit.all.table <- getCountTable(inherit.all.dt, "Affection", c("2", "1"), fam.data, ehdn.genotype.data, cov)
  
  # sm.mssng.inherit <- testGLM(inherit.mssng.table, "FAMILY_TYPE", cov, feature)
  # sm.ssc.inherit <- testGLM(inherit.ssc.table, "Affection", cov, feature)
  # sm.all.inherit <- testGLM(inherit.all.table, "Affection", cov, feature)
  # 
  # transmitted.mssng.dt <- tmp.f[which(tmp.f$repeatID %in% repeatIDs & tmp.f$FAMILY_TYPE %in% c("MPX", "SPX") & 
  #                                       tmp.f$Status %in% c("AffectedParent", "UnaffectedParent")), ]
  # transmitted.ssc.dt <- tmp.f[which(tmp.f$repeatID %in% repeatIDs & tmp.f$Dataset == "SSC" & 
  #                                     tmp.f$Status %in% c("AffectedParent", "UnaffectedParent")), ]
  # transmitted.all.dt <- tmp.f[which((tmp.f$FAMILY_TYPE %in% c("MPX", "SPX") | tmp.f$Dataset == "SSC") & 
  #                                     tmp.f$Status %in% c("AffectedParent", "UnaffectedParent")), ]
  
  all.result <- testCaseCtrl(caseCtrl.dt[caseCtrl.dt$Dataset %in% c("SSC", "MSSNG"), ], "Status", c("AffectedKid",
                                                                                                    "UnaffectedKid"),
                             case.denom, ctrl.denom, fam.data = fam.data)
  
  # all.mssng.result <- testCaseCtrl(caseCtrl.dt[caseCtrl.dt$Dataset == "MSSNG", ], "FAMILY_TYPE", c("MPX", "SPX"), mssng.fam.mpx.denom, mssng.fam.spx.denom)
  # all.ssc.result <- testCaseCtrl(caseCtrl.dt[caseCtrl.dt$Dataset == "SSC", ], "Affection", c("2", "1"), ssc.fam.case.denom, ssc.fam.ctrl.denom)
  
  # mssng.table <- getCountTable(caseCtrl.dt[caseCtrl.dt$Dataset == "MSSNG", ], "FAMILY_TYPE", c("MPX", "SPX"), fam.data, ehdn.genotype.data, cov)
  # ssc.table <- getCountTable(caseCtrl.dt[caseCtrl.dt$Dataset == "SSC", ], "Affection", c("2", "1"), fam.data, ehdn.genotype.data, cov)
  all.table <- getCountTable(caseCtrl.dt[(caseCtrl.dt$Dataset %in% c("SSC", "MSSNG")), ], "Status", 
                             c("AffectedKid", "UnaffectedKid"), 
                             fam.data, ehdn.genotype.data, cov)
  
  # sm.mssng.all <- testGLM(mssng.table, "FAMILY_TYPE", cov, feature, standardization)
  # sm.ssc.all <- testGLM(ssc.table, "Affection", cov, feature, standardization)
  sm.all.all <- testGLM(all.table, "Status", cov, feature, standardization)
  
  # names(all.mssng.result) <- paste0("MSSNG.", names(all.mssng.result))
  # names(all.ssc.result) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("SSC.", names(all.ssc.result))))
  # names(all.result) <- paste0("Aggr.", names(all.result))
  
  # mssng.denovo <- testDenovo(denovo.mssng.dt, "FAMILY_TYPE", c("MPX", "SPX"), mssng.fam.mpx.denom, mssng.fam.spx.denom)
  # ssc.denovo <- testDenovo(denovo.ssc.dt, "Affection", c("2", "1"), ssc.fam.case.denom, ssc.fam.ctrl.denom)
  # all.denovo <- testDenovo(denovo.all.dt, "Affection", c("2", "1"), case.denom, ctrl.denom)
  # 
  # mssng.inherit <- testDenovo(inherit.mssng.dt, "FAMILY_TYPE", c("MPX", "SPX"), mssng.fam.mpx.denom, mssng.fam.spx.denom)
  # ssc.inherit <- testDenovo(inherit.ssc.dt, "Affection", c("2", "1"), ssc.fam.case.denom, ssc.fam.ctrl.denom)
  # all.inherit <- testDenovo(inherit.all.dt, "Affection", c("2", "1"), case.denom, ctrl.denom)
  # 
  # mssng.ssc.denovo <- testDenovo(tmp.f[which(tmp.f$repeatID %in% repeatIDs & (tmp.f$FAMILY_TYPE == "SPX" |
  #                                                                                tmp.f$Dataset == "SSC") & 
  #                                              tmp.f$Status %in% c("AffectedKid") & tmp.f$denovo == T), ], 
  #                                "Dataset", c("MSSNG", "SSC"),mssng.fam.spx.denom, ssc.fam.case.denom, test = "fisher")
  # mssng.ssc.inherit <- testDenovo(tmp.f[which(tmp.f$repeatID %in% repeatIDs & (tmp.f$FAMILY_TYPE == "SPX" |
  #                                                                                 tmp.f$Dataset == "SSC") & 
  #                                               tmp.f$Status %in% c("AffectedKid") & tmp.f$denovo == F), ], 
  #                                 "Dataset", c("MSSNG", "SSC"),mssng.fam.spx.denom, ssc.fam.case.denom, test = "fisher")
  # 
  # denovo.mssng.ssc.table <- getCountTable(tmp.f[which(tmp.f$repeatID %in% repeatIDs & (tmp.f$FAMILY_TYPE == "SPX" |
  #                                                                                    tmp.f$Dataset == "SSC") & 
  #                                                   tmp.f$Status %in% c("AffectedKid") & tmp.f$denovo == T), ], 
  #                                     "Dataset", c("MSSNG", "SSC"), fam.data, ehdn.genotype.data, cov)
  # inherit.mssng.ssc.table <- getCountTable(tmp.f[which(tmp.f$repeatID %in% repeatIDs & (tmp.f$FAMILY_TYPE == "SPX" |
  #                                                                                    tmp.f$Dataset == "SSC") & 
  #                                                   tmp.f$Status %in% c("AffectedKid") & tmp.f$denovo == F), ], 
  #                                     "Dataset", c("MSSNG", "SSC"), fam.data, ehdn.genotype.data, cov)
  # 
  # sm.mssng.ssc.denovo <- testGLM(denovo.mssng.ssc.table, "Dataset", cov, feature, standardization)
  # sm.mssng.ssc.inherit <- testGLM(inherit.mssng.ssc.table, "Dataset", cov, feature, standardization)
  # 
  # 
  # mssng.transmission <- testTransmission(transmitted.mssng.dt, "FAMILY_TYPE", c("MPX", "SPX"))
  # ssc.transmission <- testTransmission(transmitted.ssc.dt, "Affection", c("2", "1"))
  # all.transmission <- testTransmission(transmitted.all.dt, "Affection", c("2", "1"))
  # 
  # names(mssng.ssc.denovo) <- paste0("MSSNG.SSC.Denovo.", names(mssng.ssc.denovo))
  # names(mssng.ssc.inherit) <- paste0("MSSNG.SSC.Inherit.", names(mssng.ssc.inherit))
  # 
  # names(mssng.denovo) <- paste0("MSSNG.DeNovo.", names(mssng.denovo))
  # names(ssc.denovo) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("SSC.DeNovo.", names(ssc.denovo))))
  # names(all.denovo) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("Aggr.DeNovo.", names(all.denovo))))
  # 
  # names(mssng.inherit) <- paste0("MSSNG.Inherit.", names(mssng.inherit))
  # names(ssc.inherit) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("SSC.Inherit.", names(ssc.inherit))))
  # names(all.inherit) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("Aggr.Inherit.", names(all.inherit))))
  # 
  # names(mssng.transmission) <- paste0("MSSNG.Transmitted.", names(mssng.transmission))
  # names(ssc.transmission) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("SSC.Transmitted.", names(ssc.transmission))))
  # names(all.transmission) <- gsub("2", "Case", gsub("1", "Ctrl", paste0("Aggr.Transmitted.", names(all.transmission))))
  
  # return(data.frame(all.mssng.result, all.ssc.result, all.result, mssng.ssc.denovo, mssng.ssc.inherit,
  #                   mssng.denovo, ssc.denovo, all.denovo,
  #                   mssng.inherit, ssc.inherit, all.inherit,
  #                   mssng.transmission, ssc.transmission, all.transmission,
  #                   mssng.ssc.glm.denovo.OR = sm.mssng.ssc.denovo$OR, mssng.ssc.glm.denovo.pvalue = sm.mssng.ssc.denovo$pvalue, 
  #                   mssng.ssc.glm.inherit.OR = sm.mssng.ssc.inherit$OR, mssng.ssc.glm.inherit.pvalue = sm.mssng.ssc.inherit$pvalue, 
  #                   all.mssng.glm.OR = sm.mssng.all$OR, all.mssng.glm.pvalue = sm.mssng.all$pvalue,
  #                   all.ssc.glm.OR = sm.ssc.all$OR, all.ssc.glm.pvalue = sm.ssc.all$pvalue,
  #                   all.aggr.glm.OR = sm.all.all$OR, all.aggr.glm.upper = sm.all.all$upper, all.aggr.glm.lower = sm.all.all$lower, all.aggr.glm.pvalue = sm.all.all$pvalue,
  #                   mssng.denovo.glm.OR = sm.mssng.denovo$OR, mssng.denovo.glm.pvalue = sm.mssng.denovo$pvalue,
  #                   ssc.denovo.glm.OR = sm.ssc.denovo$OR, ssc.denovo.glm.pvalue = sm.ssc.denovo$pvalue,
  #                   aggr.denovo.glm.OR = sm.all.denovo$OR, aggr.denovo.glm.upper = sm.all.denovo$upper, aggr.denovo.glm.lower = sm.all.denovo$lower, aggr.denovo.glm.pvalue = sm.all.denovo$pvalue,
  #                   mssng.inherit.glm.OR = sm.mssng.inherit$OR, mssng.inherit.glm.pvalue = sm.mssng.inherit$pvalue,
  #                   ssc.inherit.glm.OR = sm.ssc.inherit$OR, ssc.inherit.glm.pvalue = sm.ssc.inherit$pvalue,
  #                   aggr.inherit.glm.OR = sm.all.inherit$OR, aggr.inherit.glm.upper = sm.all.inherit$upper, aggr.inherit.glm.lower = sm.all.inherit$lower, aggr.inherit.glm.pvalue = sm.all.inherit$pvalue
  #                   ))
  
  return(data.frame(all.result,
                    all.aggr.glm.OR = sm.all.all$OR, all.aggr.glm.upper = sm.all.all$upper, all.aggr.glm.lower = sm.all.all$lower, all.aggr.glm.pvalue = sm.all.all$pvalue
  ))
}

getGenesetTable <- function(repeatIDs, tmp.f, fam.data, ehdn.genotype.data, cov, feature){
    
    caseCtrl.dt <- tmp.f[which(tmp.f$repeatID %in% repeatIDs & 
                                 tmp.f$Status %in% c("AffectedKid", "UnaffectedKid")), ]
    
    
    all.table <- getCountTable(caseCtrl.dt[(caseCtrl.dt$Dataset %in% c("SSC", "MSSNG")), ], "Status", 
                               c("AffectedKid", "UnaffectedKid"), 
                               fam.data, ehdn.genotype.data, cov)
    names(all.table)[2] <- feature
    return(all.table[, -c(3:4)])
}
