rm(list = ls())
#### Expansion detection by DBSCAN
#### to run
#### Rscript DBSCAN.EHdn.parallel.R --infile output_regions.min2.1000G+SSC+MSSNG+QASD.txt 
#### --outpath output/ --samplelist samples.txt

#### read arguments from command line
args = commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
  stop(call. = T, "Require more argument(s)")
}

paramNames <- grep("--", args)
paramValues <- paramNames + 1

if(length(grep("--", args[paramValues])) > 0){
  stop(call. = T, "Failed in reading arguments, '--' found in argument value")
}

if(length(args) != length(c(paramNames, paramValues))){
  message(sprintf("Ignore argument(s):%s", paste(args[-union(paramNames, paramValues)], collapse=",")))
}

params <- list()
paramNames <- gsub("--", "", args[paramNames])
for(i in 1:length(paramNames)){
  params[paramNames[i]] <- args[paramValues[i]]
}

library(dbscan)
library(ggplot2)
library(data.table)

set.seed(2000)

outpath <- params$outpath
if(length(grep("\\/$", outpath)) == 0){
  outpath <- paste0(outpath, "/")
}

if(!dir.exists(outpath)){
  dir.create(outpath, recursive = T)
}

infile <- params$infile
message(sprintf("reading %s##", infile))
dt <- fread(infile)
dt <- data.frame(dt)
dt$varid <- paste(dt$V1, dt$V2, dt$V3, sep="_")

### get count per sample
sample.count <- paste(dt$V7, collapse = ",")
sample.count <- gsub(":", ",", sample.count)
sample.count <- strsplit(sample.count, ",")[[1]] 
sample.count <- matrix(sample.count, ncol = 2, byrow = T)
sample.count <- data.frame(table(sample.count[, 1]), stringsAsFactors = F)

if("samplelist" %in% names(params)){
  samples <- readLines(params$samplelist)
  sample.count <- sample.count[sample.count$Var1 %in% samples, ]
}else{
  samples <- sample.count$Var1
}

sd <- sd(sample.count$Freq)
mean <- mean(sample.count$Freq)
outlier.count <- sum(sample.count$Freq > mean+3*sd | sample.count$Freq < mean-3*sd)

if(outlier.count > 0){
  message(sprintf("Remove %s outliers", outlier.count))
  samples <- samples[!samples %in% sample.count$Var[which(sample.count$Freq > mean+3*sd | sample.count$Freq < mean-3*sd)]]
}

message(sprintf("Samples included in the analysis are in %sclean.sample.txt", outpath))
writeLines(as.character(sample.count$Var1), sprintf("%sclean.sample.txt", outpath))

norm.dist <- round(rnorm(length(samples), mean = 1, sd = 0.25), digits = 1)
norm.dist <- sort(norm.dist)

library(parallel)
library(doSNOW)

time.start <- Sys.time()
coreNumber <- 12
cl <- makeCluster(coreNumber-1)
registerDoSNOW(cl)

dt.out <- foreach (i = 1:nrow(dt), .combine=rbind) %dopar% {
  tmp <- gsub(":|,", "#", dt[i, 7])
  tmp <- unlist(strsplit(tmp, "#")[[1]])
  tmp <- data.frame(matrix(tmp, ncol = 2, byrow = T), stringsAsFactors = F)
  tmp[, 2] <- as.numeric(tmp[, 2])
  tmp <- tmp[tmp$X1 %in% samples, ]
  if(nrow(tmp) == 0){
    dt.here <- data.frame("repeatID" = "", "motif" = "", "outliers" = "", "size" = "", "ref" = "",
                          "chr" = "", "start" = "", "end" = "")[-1, ]
  }else{
    if(nrow(tmp) < length(norm.dist)){
      dt.tmp <- norm.dist[-c(1:nrow(tmp))]
      sampleIDs <- c(paste0("Sim", 1:length(dt.tmp)), tmp$X1)
      dt.tmp <- c(dt.tmp, tmp$X2)
    }else{
      sampleIDs <- tmp$X1
      dt.tmp <- tmp$X2
    }
    
    ref <- as.numeric(names(which.max(table(dt.tmp))))
    range <-  max(2*ref, quantile(dt.tmp, 0.95, na.rm = T) - quantile(dt.tmp, 0.05, na.rm = T))
    
    scan <- dbscan::dbscan(matrix(dt.tmp), eps = range, minPts = ceiling(log2(length(samples))))
  
    if(length(unique(scan$cluster)) == 1 | sum(scan$cluster == 0) == 0){
      cutoff <- Inf
    }else{
      cutoff <- max(dt.tmp[scan$cluster != 0])
      cutoff <- ifelse(cutoff < 2, 2, cutoff)
    }
    
    if(sum(dt.tmp > cutoff) > 0){
      sim.samples <- grep("Sim", sampleIDs)
      if(length(sim.samples) > 0){
        dt.tmp <- dt.tmp[-sim.samples]
        sampleIDs <- sampleIDs[-sim.samples]
      }
      
      outliers <- paste(sampleIDs[dt.tmp > cutoff], collapse=";")
      size <- paste(dt.tmp[dt.tmp > cutoff], collapse=";")
      dt.here <- data.frame("repeatID" = dt$varid[i], "motif" = dt$V4[i], outliers, size, ref,
                            "chr" = dt$V1[i], "start" = dt$V2[i], "end" = dt$V3[i])
    }else{
      dt.here <- data.frame("repeatID" = "", "motif" = "", "outliers" = "", "size" = "", "ref" = "",
                            "chr" = "", "start" = "", "end" = "")[-1, ]
    }
  }

  dt.here
}


stopCluster(cl) 
write.table(dt.out, sprintf("%sEHdn.expansions.%s.tsv", outpath, Sys.Date()), sep="\t", row.names=F, col.names=T, quote=F)
time.end <- Sys.time()
print(difftime(time.end, time.start))
