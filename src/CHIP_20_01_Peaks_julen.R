source("src/00_init.R")

# Settings ----------------------------------------------------------------
out <- dirout("CHIP_20_01_Peaks_julen/")

ff <- list.files(paste(PATHS$LOCATIONS$DATA, "ChIP_Peaks_Julen","v2", sep="/"), recursive = TRUE, full.names = TRUE)

fx <- "Brd9"
res <- list()
res.fc <- data.table()
for(fx in gsub("_.+", "", basename(grep("extended.motif.fc.txt$", ff, value = TRUE)))){
  message(fx)
  peaks <- fread(ff[grepl(fx, ff) & grepl("extended.motif.fc.txt$", ff)], check.names = FALSE)
  colnames(peaks) <- make.names(gsub("\\+", "plus", colnames(peaks)))
  # motifs <- fread(ff[grepl(fx, ff) & grepl("motifsInPeak", ff)])
  
  xc <- grep(".+.bool", colnames(peaks),value=TRUE)
  xx <- xc[1]
  for(xx in xc){
    
    #xnam <- gsub("_.+$", "", gsub("\\..+$", "", sub("_", " ", xx)))
    xnam <- gsub(paste0(fx, ".+$"), fx, xx)
    print(xx)
    message(xnam)
    
    if(xnam %in% names(res)){
      message("Already DONE")
      next
    }
    
    xx.fc <- gsub("\\.bool$", ".fc", xx)
    
    x <- peaks[,c("chr", "start", "end", "interval_id", "Distance.to.TSS", "Annotation", "Gene.Name", xx, xx.fc),with=F]
    x <- x[get(xx) == TRUE]
    x <- x[grepl("promoter-TSS", Annotation)]
    
    x[[xx.fc]] <- sapply(strsplit(as.character(x[[xx.fc]]), ";"), function(j){
      max(as.numeric(j))
    })
    
    # Export FCs
    ret1 <- setNames(x[, c("Gene.Name", xx.fc), with=F], c("Gene", "FC"))
    ret1$ChIP <- xnam
    res.fc <- rbind(res.fc, ret1)
    
    # Export target genes
    x <- x[get(xx.fc) > 5,]
    res[[xnam]] <- unique(x$Gene.Name)
  }
}
sapply(res, length)
chip.targets <- res
chip.targets.fc <- res.fc
chip.targets.fc$FC <- unlist(chip.targets.fc$FC)
save(chip.targets, chip.targets.fc, file=out("ChIP.Targets.RData"))
