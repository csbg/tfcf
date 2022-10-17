#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#===============================================================================
# Aim:
# a) Make table with target sgRNAs 
# b) Make a table with "noise" = sgRNAs that shouldn't be there. 
# Then calculate: 
# 1) Read counts/sample 
#===============================================================================
## PACKAGES.
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(lazyeval)
library(dplyr)
library(tibble)
#===============================================================================
# How to run me
#Rscript --vanilla 01a_MECC_sgRNA_pre_analysis_v2.r \ 
#/PATH/TO/00_MECC_sgRNA_v1.1.sh/OUTPUT > runID.Rout.txt

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {
  stop("Too many arguments supplied (input file).n", call.=FALSE)
}
  
## GLOBAL VARIABLES.
## TO BE Modified
PROJECT_DIR <- args[1]
#PROJECT_DIR<-file.path("/home/X/data/2021/CRISPR/allProcessed/merge4-492")

guidesFile <- "/home/jmendietaes/data/2021/CRISPR/finalGuides.txt"

# If the folder structure is ok this shouldnt change
Counts<-file.path(PROJECT_DIR,"idxstats")
RSession<-file.path(PROJECT_DIR,"RSession")
runName <- basename(PROJECT_DIR)

## Relevant info about the guides
# List of all allowed guides
allGuideCodes <- c("MGLibA", "R2.Br", "As", "B.Br", "TF1.Br", "Br", "R1.Br", "B", 
"TF2.Br", "R2Br", "BBr", "TF1Br", "R1Br", "LibP1", "LibP2", "OP1", "OP2", "OP3",
"OP4", "OP1234")
# Make synonyms for each time they use a different name for tha same lib
guideSynonims <- vector(mode="list", length=length(allGuideCodes) + 5)
names(guideSynonims) <- c("TF1", "TF2", "R1", "R2", "BBr", allGuideCodes)
guideSynonims[[1]] <- "TF1Br"; guideSynonims[[2]] <- "TF2Br"
guideSynonims[[3]] <- "R1Br"; guideSynonims[[4]] <- "R2Br"
guideSynonims[[5]] <- "BBr"

for (i in 1:length(allGuideCodes)) {
  guideSynonims[[i+5]] <- allGuideCodes[i]
}

#===============================================================================
######## FUNCTIONS. (the first is just an scaffold example from a known analysis)
#' @param x Character vector of length 1, f
#' @return Character vector of length 1, 
#' @export

multimerge <- function (mylist) {
  ## mimics a recursive merge or full outer join
  unames <- unique(unlist(lapply(mylist, rownames)))
  n <- length(unames)
  out <- lapply(mylist, function(df) {
    tmp <- matrix(nr = n, nc = ncol(df), dimnames = list(unames,colnames(df)))
    tmp[rownames(df), ] <- as.matrix(df)
    rm(df); gc()
    return(tmp)
  })
  stopifnot( all( sapply(out, function(x) identical(rownames(x), unames)) ) )
  bigout <- do.call(cbind, out)
  colnames(bigout) <- paste(rep(names(mylist), sapply(mylist, ncol)), unlist(sapply(mylist, colnames)), sep = "_")
  return(bigout)
}


#===============================================================================
### MAIN.
################################################################################
## Step . data loading
################################################################################

# Set the directory to where your "extracted.fastq.bam.cnt.txt.final.txt" are located
setwd (PROJECT_DIR)

## Read each of the files in a different table. See example below:
# first we get all the files in folder and their associated libraries
countreports <- dir(file.path(Counts), pattern = "*.idxstats")
allLibs <- unique(sapply(strsplit(countreports, "_"), function(x) x[3]))
# sometimes we have files with more than one library separated by '-'
allLibs <- unique(unlist(strsplit(allLibs, "-")))
# all must be in guide or synonim list
if (! sum(allLibs %in% c(allGuideCodes, names(guideSynonims))) == length(allLibs)) {
  print('Some filetered guide libraries are not in the reference list')
  print(allLibs)
  # EXIT with error
  q("no", 1, FALSE)
}

for (glib_ in allLibs) {
  # get the real ID in case he file had a synonim
  glib = guideSynonims[[glib_]]
  noDotLib <- gsub("\\.", "", glib, perl=TRUE)
  print("")
  print(paste0('---------------------- ', noDotLib, ' ----------------------'))
  
  # Load the list that relates guides to each library
  #guidesInfo <- read.table(guidesFile, stringsAsFactors = FALSE, header = TRUE) 
  #guidesInfo <- guidesInfo[,c('ID', 'Library')]
  
  #libPos <- grepl(glib, guidesInfo[,2])
  #libGuides <- guidesInfo[libPos,1]
  
  # get the list of files for this library
  gcountreports <- countreports[grepl(paste0("_", glib_, "_"), countreports)]
  # now that we have two files per ID (mapped to lib and unmapped to lib) we 
  # have to substract the base name of them
  gcountreports <- unique(gsub("\\.unMap.idxstats", '', gcountreports, perl=TRUE))
  gcountreports <- unique(gsub("\\.idxstats", '', gcountreports, perl=TRUE))
  
  ## get % of guide extraction
  statsFile <- paste0(PROJECT_DIR, '/RSession/', runName, "_sequencingStats.tsv")
  statsTable <- read.table(statsFile, sep='\t', header = TRUE)
  statsTable <- statsTable[c(1, 4)]
  # remove the percentaje symbol
  statsTable$Percentaje<-gsub("%","",as.character(statsTable$Percentaje))
  # kee only the IDs from this library
  statsTable <- statsTable[grepl(paste0("_", glib_, "_"), statsTable$SampleID), ]
  statsTable2 <- data.frame(statsTable)
  rownames(statsTable) <- statsTable$SampleID
  statsTable <- transform(statsTable, Percentaje = as.numeric(Percentaje))

  # plot guide presence percentajes
  pdf(paste0(RSession, "/", runName, "_", noDotLib, "_stats.pdf"), width=11, height=8.5)
  
  par(mar=c(15, 4.1, 4.1, 2.1), cex=1)
  # Basic barplot
  p<-ggplot(data=statsTable, aes(x=SampleID, y=Percentaje, label = Percentaje)) +
    geom_bar(stat="identity") + 
    ggtitle(paste0("Percentaje of recovered guide patterns (even if no valid guide inside): ", noDotLib)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylim(0, 100) +
    geom_text(size = 3, position = position_stack(vjust = 0.5))
  print(p)
  
 
  
  ## Load info of guide smapped to main library
  # read mapped files in memory
  for (reporte in gcountreports){
    reporte <- paste0(reporte, ".idxstats")
    t <- read.table(file.path(Counts, reporte), stringsAsFactors = FALSE,
                    col.names= c("sgRNA","length","counts", "unmapped"), 
                    row.names = 1) 
    row.names(t)[nrow(t)] <- "unmapped"
    t[nrow(t),"counts"] <- t[length(rownames(t)), "unmapped"]
    t <- t[,"counts", drop = FALSE]
    assign(reporte, t, envir = .GlobalEnv)
    rm(t)
  }
  countable_list <- Filter(function(x) is(x, "data.frame"), mget(ls()))
  # keep only the data frames that contain the idxstats info
  countable_list <- countable_list[grepl("idxstats", names(countable_list))]
  
  
  # Use the multimerge function to merge them together into a single table.
  table<- multimerge(countable_list)
  colnames(table) <- gsub ("_S.*.idxstats_counts","", colnames(table))
  # This step is to remove the dot from old library naming
  colnames(table) <- gsub ("\\.","", colnames(table))
  
  # turn values to numbers
  table <- data.frame(apply(table, 2, function(x) as.numeric(x)), row.names = rownames(table))
  colnames(table) <- gsub("\\.","", colnames(table))
  # delete lists form environment
  rm(list = ls(pattern = ".idxstats"))
  # Replace NAs for zeros
  table[is.na(table)] <- 0 
  # Write table with alignment info to all Libraries
  dir.create(RSession, showWarnings = FALSE)
  write.table(data.frame("ID"=rownames(table),table), 
                      paste0(RSession, "/", runName, "_", 
                      noDotLib, "_crisprTable_raw_MapLib.tsv"), 
              row.names=FALSE, quote=FALSE, sep='\t')
  

  ## Load unMapped that were aligned against all guides
  # read mapped files in memory
  for (reporte in gcountreports){
    reporte <- paste0(reporte, ".unMap.idxstats")
    t <- read.table(file.path(Counts, reporte), stringsAsFactors = FALSE,
                    col.names= c("sgRNA","length","counts", "unmapped"), 
                    row.names = 1) 
    row.names(t)[nrow(t)] <- "unmapped"
    t[nrow(t),"counts"] <- t[length(rownames(t)), "unmapped"]
    t <- t[,"counts", drop = FALSE]
    assign(reporte, t, envir = .GlobalEnv)
    rm(t)
  }
  countable_list <- Filter(function(x) is(x, "data.frame"), mget(ls()))
  # keep only the data frames that contain the idxstats info
  countable_list <- countable_list[grepl("idxstats", names(countable_list))]
  
  # Use the multimerge function to merge them together into a single table.
  nonValidTable<- multimerge(countable_list)
  colnames(nonValidTable) <- gsub ("_S.*.unMap.idxstats_counts","", colnames(nonValidTable))
  # This step is to remove the dot from old library naming
  colnames(nonValidTable) <- gsub ("\\.","", colnames(nonValidTable))
  
  # turn values to numbers
  nonValidTable <- data.frame(apply(nonValidTable, 2, function(x) as.numeric(x)), row.names = rownames(nonValidTable))
  colnames(nonValidTable) <- gsub ("\\.","", colnames(nonValidTable))
  # delete lists form environment
  rm(list = ls(pattern = ".idxstats"))
  # Replace NAs for zeros
  nonValidTable[is.na(nonValidTable)] <- 0 
  # Write table with alignment info to all Libraries
  dir.create(RSession, showWarnings = FALSE)
  write.table(data.frame("ID"=rownames(nonValidTable),nonValidTable), 
                        paste0(RSession, "/", runName, "_", 
                        noDotLib, "_crisprTable_raw_unMapLib.tsv"), 
              row.names=FALSE, quote=FALSE, sep='\t')
  
  
  ## Plot
  # when there is only one table it doesnt work well
  if (ncol(table) == 1) {
    stackedPlotTabl <- read.csv(text=paste(table[nrow(table),], 
                                           nonValidTable[nrow(nonValidTable),], 
                                           colnames(table), sep=','),  header=FALSE)
    rownames(stackedPlotTabl) <- colnames(nonValidTable)
    colnames(stackedPlotTabl) <- c("Unmapped to valid", "Unmapped to non-valid", "Ids")
  } else {
    tb1 <- as.data.frame(table[nrow(table),])
    tb2 <- as.data.frame(nonValidTable[nrow(nonValidTable),])
    
    stackedPlotTabl <- rbind(tb1, tb2)
    rownames(stackedPlotTabl) <- c("Unmapped to valid", "Unmapped to non-valid")
    stackedPlotTabl <- as.data.frame(t(stackedPlotTabl))
    stackedPlotTabl$Ids <- rownames(stackedPlotTabl)
    
  }
  
  stackedPlotTabl <- gather(stackedPlotTabl, reads, Unmapped, 1:2)
  stackedPlotTabl$reads <- factor(stackedPlotTabl$reads,
                                  levels = c("Unmapped to valid", "Unmapped to non-valid"))

  cols <- c("black", "darkgrey")
  par(mar=c(14,5,1,1))
  p <- ggplot(stackedPlotTabl, aes(fill=reads, y=Unmapped, x=Ids)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0("Reads unmapped to valid, and then to unvalid library: ", noDotLib)) +
    ylab("N read") +
    scale_fill_manual(values= cols)
    
  print(p)
  
  # new stacked plot
  if (length(colnames(table)) == 1) {
    stackedPlotTabl <- read.csv(text=paste(colnames(table),
                                           "Mapped_to_valid",
                                           paste0(sum(table[1:nrow(table)-1,]), '\n',
                                                  colnames(table)), 
                                           "nonValid_to_all",
                                           paste0(sum(nonValidTable[1:nrow(nonValidTable)-1,]),
                                                  "\n", colnames(table)),
                                           "Remaining_unmap",
                                           nonValidTable[nrow(nonValidTable),], 
                                           sep=','),  header=FALSE)
    colnames(stackedPlotTabl) <- c("Ids", "reads", "Unmapped")
    
    stackedPlotTabl[, "mapPercent"] <- ''
    for (c in colnames(table)){
      pos <- stackedPlotTabl[,1] == c & stackedPlotTabl[,2] == "Mapped_to_valid"
      m1 <- stackedPlotTabl[pos, 3]
      pos <- stackedPlotTabl[,1] == c & stackedPlotTabl[,2] == "nonValid_to_all"
      m2 <- stackedPlotTabl[pos, 3]
      pos <- stackedPlotTabl[,1] == c & stackedPlotTabl[,2] == "Remaining_unmap"
      m3 <- stackedPlotTabl[pos, 3]
      mapPercent <- round((m1 / (m1 + m2 + m3)) * 100, 2)
      pos <- stackedPlotTabl[,1] == c & stackedPlotTabl[,2] == "Mapped_to_valid"
      stackedPlotTabl[pos, "mapPercent"] <- mapPercent
      
    }
  } else {
    mapedGuide <- t(as.data.frame(apply(table[1:nrow(table)-1,], 2, sum)))
    rownames(mapedGuide) <- c("Mapped_to_valid")
    mapedElse <- t(as.data.frame(apply(nonValidTable[1:nrow(nonValidTable)-1,], 2, sum)))
    rownames(mapedElse) <- c("nonValid_to_all")
    remaining <- as.data.frame(nonValidTable[nrow(nonValidTable),])
    rownames(remaining) <- c("Remaining_unmap")
    
    stackedPlotTabl <- rbind(mapedGuide, mapedElse, remaining)
    #rownames(stackedPlotTabl) <- c("Mapped_to_valid", "nonValid_to_all", "Remaining_unmap")
    stackedPlotTabl <- as.data.frame(t(stackedPlotTabl))
    stackedPlotTabl$Ids <- rownames(stackedPlotTabl)
    stackedPlotTabl <- gather(stackedPlotTabl, reads, Unmapped, 1:3)
    
    stackedPlotTabl[, "mapPercent"] <- ''
    for (c in colnames(mapedGuide)){
      m1 <- mapedGuide[,c]
      m2 <- mapedElse[,c]
      m3 <- remaining[,c]
      mapPercent <- round((m1 / (m1 + m2 + m3)) * 100, 2)
      pos <- stackedPlotTabl[,1] == c & stackedPlotTabl[,2] == "Mapped_to_valid"
      stackedPlotTabl[pos, "mapPercent"] <- mapPercent
      
    }
  }
  
  
  par(mar=c(14,5,1,1))
  p <- ggplot(stackedPlotTabl, aes(fill=reads, y=Unmapped, x=Ids, label = mapPercent)) + 
    geom_bar( stat="identity") +
    scale_fill_manual(values= c(cols, "darkred")) +
    geom_text(size = 3, position = position_stack(vjust = 0.5), col='white') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0("N of reads mapping to valid, then rest, and then none: ", glib)) +
    ylab("N read")
  print(p)
  
  # we now remove the unmaped reads list
  if (length(ncol(table)) == 1) {
    validTable <- data.frame(table[-nrow(table),])
    rownames(validTable) <- rownames(table)[-nrow(table)]
    colnames(validTable) <- colnames(table)
    
    nonValidTable_ <- data.frame(nonValidTable[-nrow(nonValidTable),])
    rownames(nonValidTable_) <- rownames(nonValidTable)[-nrow(nonValidTable)]
    colnames(nonValidTable_) <- colnames(nonValidTable)
    nonValidTable <- nonValidTable_
    
  } else {
  validTable<-table[-nrow(table),]
  nonValidTable<-nonValidTable[-nrow(nonValidTable),]
  }
  
  # Write the filtered elements for the valid
  write.table(data.frame("ID"=rownames(validTable),validTable), 
              paste0(RSession, "/", runName, "_", glib, 
                     "_crisprTable_noUnmap_MapLib.tsv"), 
              row.names=FALSE, quote=FALSE, sep='\t')
 
  # remove all zeros
  if (ncol(table) == 1) {
    removed1 <- validTable[,1] == 0
    validrows <- rownames(validTable)
    validrows <- validrows[validTable[,1] > 0]
    validcols <- colnames(validTable)
    validTable<- data.frame(validTable[validTable[,1] > 0,])
    rownames(validTable) <- validrows
    colnames(validTable) <- validcols
    
    removed2 <- nonValidTable[,1] == 0
    validrows <- rownames(nonValidTable)
    validrows <- validrows[nonValidTable[,1] > 0]
    validcols <- colnames(nonValidTable)
    nonValidTable<- data.frame(nonValidTable[nonValidTable[,1] > 0,])
    rownames(nonValidTable) <- validrows
    colnames(nonValidTable) <- validcols
    
  } else {
    removed1 <- apply(validTable,1,max)==0
    validTable<-validTable[apply(validTable,1,max)>0,]
    
    removed2 <- apply(nonValidTable,1,max)==0
    nonValidTable<-nonValidTable[apply(nonValidTable,1,max)>0,]
    
  }
  
  print(paste0(sum(removed1), " libraries removed in valid library map for having zero counts in all samples"))
  print(paste0(sum(removed2), " libraries removed for having zero counts in all samples after mapping unmmaped reads to NON-valid libraries"))
  
  # Write the filtered elements without zeros for the non-valid
  write.table(data.frame("ID"=rownames(nonValidTable),nonValidTable), 
              paste0(RSession, "/", runName, "_", glib, 
                     "_crisprTable_noUnmap_unMapLib_nonZero.tsv"), 
              row.names=FALSE, quote=FALSE, sep='\t')
  
  # Remove guides with less than 1000 counts in all the samples
  #removed <- rownames(validTable[!apply(validTable,1,max)>1000,])
  #print(paste0(length(removed), ' guides removed for having less than 1000 counts in all samples'))
  #print(removed)
  #validTable<-validTable[apply(validTable,1,max)>1000,]
  #plot(density(unlist(log2(validTable+1))))
  
  
  #pdf(paste0(RSession, "/Density_", glib, "_valid.pdf"), width=11, height=8.5)
  plot(density(unlist(log2(validTable+1))), 
       main = paste0("Density ", glib, ": reads that mapped to focus library"), 
       col=cols[1], lwd = 3)
  #dev.off()
  
  #pdf(paste0(RSession, "/Density_", glib, "_NonValid.pdf"), width=11, height=8.5)
  plot(density(unlist(log2(nonValidTable+1))), 
       main = paste0("Density ", glib, ": unmapped reads that mapped to non-valid libraries"),
       col=cols[2], lwd = 3)
  #dev.off()
  
  # together
  #pdf(paste0(RSession, "/CountDistrib_", glib, "_valid_norm.pdf"), width=11, height=8.5)
  par(mfrow = c(1, 1))
  # You can change the xlim and ylim to fit your distribution
  plot(NA, ylim = c(0, 1), xlim = c(2, 18), ylab="Distribution density",
       main=paste0("Merged mapping reads: ", glib))
  lines (density(unlist(log2(validTable+1))), col=cols[1], lwd=3)
  lines (density(unlist(log2(nonValidTable+1))), col=cols[2], lwd=3)
  legend("topright", legend=c("MapLib", "unMapLib"), lwd=2, col=cols)

  
  #pdf(paste0(RSession, "/Density_", glib, "_NonValid_after1000.pdf"), width=11, height=8.5)
  #plot(density(unlist(log2(nonValidTable[apply(nonValidTable,1,max)>1000,]+1))),
  #     main = paste0("Density", glib, "_nonValid_min1000"))
  #dev.off()
  
  #Normalize data ##
  # Calculate "factor.sizes"
  factor.sizes<-1/(colSums(na.omit(validTable))/ mean(colSums(na.omit(validTable)))) 
  # Apply normalisation
  validTable.norm<-t(t(validTable) * factor.sizes) 
  validTable.norm<-as.data.frame(validTable.norm)
  
  write.table(data.frame("ID"=rownames(validTable.norm),validTable.norm), 
              paste0(RSession, "/", runName, "_", glib, 
                     "_crisprTable_normOld_MapLib.tsv"), 
              row.names=FALSE, quote=FALSE, sep='\t')
  
  ################################################################################
  ## Step . Basic stats
  ################################################################################
  ## For valid reads ##
  ## 1) Read counts/sample
  #Check the sgRNA count distribution:
  # RAW
  #pdf(paste0(RSession, "/readCount_", glib, "_valid_raw.pdf"), width=11, height=8.5)
  par(mar=c(15, 4.1, 4.1, 2.1), cex=1)
  barplot(colSums(validTable), las=2, main=paste0("Total counts per sample: ", glib))
  #dev.off()
  # NORM
  #pdf(paste0(RSession, "/CountBoxplot_", glib, "_valid_norm.pdf"), width=11, height=8.5)
  #par(mar=c(15, 4.1, 4.1, 2.1), cex=0.5)
  #barplot(colSums(validTable.norm), las=2, main=paste0("Total normalised counts per sample: ", glib))
  #dev.off()
  
  
  ## 2) boxplots (in log scale)
  # RAW
  #pdf(paste0(RSession, "/CountBoxplot_", glib, "_valid_raw.pdf"), width=11, height=8.5)
  par(mar=c(15, 4.1, 4.1, 2.1), cex=1)
  boxplot(log2(validTable+1),las=2, main=paste0("Count boxplot per sample: ", glib))
  #dev.off()
  # NORM
  #pdf(paste0(RSession, "/CountBoxplot_", glib, "_valid_norm.pdf"), width=11, height=8.5)
  par(mar=c(15, 4.1, 4.1, 2.1), cex=1)
  boxplot(log2(validTable.norm+1),las=2, main=paste0("Normalised count boxplot per sample: ", glib))
  #dev.off()
  
  
  ## 3) density plots for each of them
  # RAW
  #pdf(paste0(RSession, "/CountDistrib_", glib, "_valid_raw.pdf"), width=11, height=8.5)
  cols <- rainbow(dim(validTable)[2]) 
  par(mfrow = c(1, 1))
  # You can change the xlim and ylim to fit your distribution
  plot(NA, ylim = c(0, 1), xlim = c(2, 18), 
       main=paste0("sgRNA Distribution before normalization: ", glib),
       ylab="Distribution density", xlab='log2(library counts)') 
  for (i in 1:ncol(validTable)){
    lines (density(log2(na.omit(validTable[,i]))+1), col=cols[i], lwd=2)
  }
  legend("topright", legend=colnames(validTable), lwd=2, col=cols)
  #dev.off()
  # NORM
  #pdf(paste0(RSession, "/CountDistrib_", glib, "_valid_norm.pdf"), width=11, height=8.5)
  cols <- rainbow(dim(validTable.norm)[2]) 
  par(mfrow = c(1, 1))
  # You can change the xlim and ylim to fit your distribution
  plot(NA, ylim = c(0, 1), xlim = c(2, 18), 
       main=paste0("sgRNA Distribution after normalization: ", glib),
       ylab="Distribution density", xlab='log2(normalized library counts)') 
  for (i in 1:ncol(validTable.norm)){
    lines (density(log2(na.omit(validTable.norm[,i]))+1), col=cols[i], lwd=2)
  }
  legend("topright", legend=colnames(validTable.norm), lwd=2, col=cols)
  dev.off()
  
  # read percentyles valid vs non valid
  print('Max read count per guide percentyle in valid')
  print(quantile(apply(validTable,1,max))) # Explore distribution
  print('Max read count per guide percentyle in non-valid')
  print(quantile(apply(nonValidTable,1,max))) # Explore distribution
  
  # show the non valids with more than 10000 counts
  #p95 <- quantile(apply(nonValidTable,1,max), probs=c(0.95))
  print('The following non-valid guides have more thn 10k counts')
  print(nonValidTable[rowSums(nonValidTable > 10000) > 0, ])
  
}
