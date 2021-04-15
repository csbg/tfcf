source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
baseDir <- "POOLED_01_CollectData/"
out <- dirout(baseDir)

require(readxl)
require(edgeR)
require(limma)
HM.COLORS.FUNC <- colorRampPalette(c("#6a3d9a", "#a6cee3", "white", "#fdbf6f", "#e31a1c"))



# read data files ---------------------------------------------------------------
ff <- list.files(paste0(Sys.getenv("DATA"), "raw_pooled/"), pattern=".txt$", full.names = TRUE)
names(ff) <- basename(ff)
ff <- lapply(ff, fread)
fnam <- names(ff)[5]
dDT <- setNames(lapply(names(ff), function(fnam){
  fx <- ff[[fnam]]
  fx <- fx[,!grepl("_score", colnames(fx)) & !grepl(".Norm", colnames(fx)) & !grepl("_norm", colnames(fx)), with=F]
  names(fx) <- paste(names(fx), fnam)
  names(fx)[grepl("^V1", names(fx))] <- "V1"
  fx
}), names(ff))


# Clean up FcgR -----------------------------------------------------------
fnam <- "LibAs_Mye-Nov2019.txt"
for(fnam in names(dDT)){
  dt <- dDT[[fnam]]
  if(any(grepl("FcgR", colnames(dt)))){
    print(paste("Fixing FcgR in ", fnam))
    fcgrColumns <- grep("FcgR", colnames(dt), value=TRUE)
    myeColumns <- gsub("FcgR\\d*_", "Mye_", fcgrColumns)
    stopifnot(all(fcgrColumns %in% colnames(dt)))
    stopifnot(all(myeColumns %in% colnames(dt)))
    res <- dt[,-fcgrColumns, with=F][,-myeColumns, with=F]
    for(i in 1:length(fcgrColumns)){
      res[[myeColumns[i]]] <- dt[[myeColumns[i]]] * 0.2 + dt[[fcgrColumns[i]]] * 0.8
    }
    dDT[[fnam]] <- res
  }
}

lapply(data.table(do.call(rbind,strsplit(gsub(" .+", "", grep("^V1$", do.call(c, lapply(dDT, colnames)), invert = TRUE, value = TRUE)), "_"))), unique)

# Extract guides into one matrix ----------------------------------------------------------
lapply(dDT, function(dt) head(dt$V1))
dDT2 <-lapply(dDT, function(dt){
  dt[,V1 := gsub("(\\d+)_(\\d{6})$", "\\1", V1)]
  dt[,V1 := gsub("MGLib", "", V1)]
  dt <- dt[!V1 %in% c("unmapped", "unmaped")]
  dt
})
#& !grepl("NonTargetingControl", V1),
dMT <- dDT2[[1]]
for(i in 2:length(dDT2)){
  dMT <- merge(dMT, dDT2[[i]], by="V1", all=TRUE)
}
m <- as.matrix(dMT[,-"V1",with=F])
row.names(m) <- dMT$V1
unique(gsub("_\\d+$", "", row.names(m)))
unique(gsub("^.+_(.+)_\\d+$", "\\1", row.names(m)))


# Read Annotation --------------------------------------------------------------
ann <- data.table(sample=colnames(m), do.call(rbind, strsplit(gsub("Wt", "WT", gsub(" ", "_", gsub("\\-", "_", colnames(m)))), "_")))
ann[V4 == "Jul2020", V2 := V3]
ann[V4 == "Jul2020", V3 := "B"]
lapply(ann, unique)
ann.col = data.frame(row.names=ann$sample, ann[,c("V1", "V2", "V3", "V6"), with=F])

# Check that names are unique (?)
x <- unique(data.table(apply(!is.na(m), 2, sum), gsub("^.+? ", "", colnames(m))))
x$V3 <- sapply(dDT2, nrow)[x$V2]
stopifnot(nrow(x[V1 != V3]) == 0)


# Correlation and guide overlaps ------------------------------------------
cleanDev(); pdf(out("Correlation.pdf"), w=12, h=11)
#cMT <- corS(m[,ann[V3 !="B" & !grepl("\\d{4}19", V4)]$sample], use="pairwise.complete.obs")
cMT <- corS(m[!grepl("NonTargetingControl", row.names(m)),ann$sample], use="pairwise.complete.obs")
diag(cMT) <- NA
pheatmap(cMT, cluster_rows = F, cluster_cols = F, annotation_col = ann.col, 
         breaks=seq(-1,1, 0.01), color=HM.COLORS.FUNC(200))
dev.off()


cleanDev(); pdf(out("Jaccard.pdf"), w=12, h=11)
pheatmap(jaccard(lapply(data.table(m), function(x) row.names(m)[!is.na(x)])), 
         cluster_rows = F, cluster_cols = F, annotation_col = ann.col, 
         breaks=seq(-1,1, 0.01), color=HM.COLORS.FUNC(200),
         main="Guide overlaps, jaccard coefficient")
dev.off()


x <- melt(data.table(m, keep.rownames = TRUE), id.vars = "rn")
x <- merge(x, ann, by.x="variable", by.y="sample")

cleanDev(); pdf(out("AllValues.pdf"), w=15, h=29)
pheatmap(m, cluster_rows = F, cluster_cols = F, fontsize_row = 2)
dev.off()

ann$sample <- make.names(ann$sample)
colnames(m) <- make.names(colnames(m))
stopifnot(all(ann$sample %in% colnames(m)))
ann <- setNames(ann, c("sample", "Genotype", "Population", "Library", "Date", "Library2", "System", "Date2"))
ann[,Date2 := gsub("\\.txt", "", Date2)]
write.table(m[,ann$sample], quote = F, sep = ",", row.names = TRUE, col.names = TRUE, file = out("Matrix.csv"))
write.tsv(ann, out("Annotation.tsv"))
