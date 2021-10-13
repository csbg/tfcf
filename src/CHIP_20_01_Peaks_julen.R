source("src/00_init.R")

# Settings ----------------------------------------------------------------
out <- dirout("CHIP_20_01_Peaks_julen/")

ff <- list.files(paste(Sys.getenv("DATA"), "ChIP_Peaks_Julen", sep="/"), recursive = TRUE, full.names = TRUE)

fx <- "Brd9"
res <- list()
for(fx in gsub("_.+", "", basename(grep("extended.txt$", ff, value = TRUE)))){
  peaks <- fread(ff[grepl(fx, ff) & grepl("extended.txt$", ff)], check.names = FALSE)
  colnames(peaks) <- make.names(gsub("\\+", "plus", colnames(peaks)))
  # motifs <- fread(ff[grepl(fx, ff) & grepl("motifsInPeak", ff)])
  
  xc <- grep(".+.bool", colnames(peaks),value=TRUE)
  xx <- xc[1]
  for(xx in xc){
    #xnam <- gsub("_.+$", "", gsub("\\..+$", "", sub("_", " ", xx)))
    xnam <- gsub(paste0(fx, ".+$"), fx, xx)
    print(xx)
    message(xnam)
    
    x <- peaks[,c("chr", "start", "end", "interval_id", "Distance.to.TSS", "Annotation", "Gene.Name", xx),with=F]
    x <- x[get(xx) == TRUE]
    #unique(x[abs(Distance.to.TSS) < 500]$Gene.Name)
    res[[xnam]] <- unique(x[grepl("promoter-TSS", Annotation)]$Gene.Name)
  }
}
sapply(res, length)
