source("src/00_init.R")
baseDir <- "POOLED_01_CollectData/"
out <- dirout(baseDir)

require(readxl)
require(edgeR)
require(limma)
HM.COLORS.FUNC <- colorRampPalette(c("#6a3d9a", "#a6cee3", "white", "#fdbf6f", "#e31a1c"))

cleanNames <- function(x){
  x <- gsub("GMP.CD11bGr1", "GMP.CD11bGr1", x)
  x <- gsub("GMP.DN", "GMP.DN", x)
  x <- gsub("GMP.Mye", "GMP.Mye", x)
  x <- gsub("\\-", "_", x)
}

# read data files ---------------------------------------------------------------
ff <- list.files(paste0(PATHS$LOCATIONS$RAWDATA, "POOLED/v4_final/"), recursive = TRUE, pattern="", full.names = TRUE)
ff <- ff[!grepl("\\/DM\\/", ff) | grepl("Time0", ff)] # From the first DM folder keep only Time 0
ff <- ff[grepl(".txt$", ff) | grepl(".tsv$", ff)]
stopifnot(sum(duplicated(basename(ff))) == 0)
names(ff) <- cleanNames(basename(ff))
ff <- lapply(ff, fread)
fnam <- names(ff)[2]
dDT <- setNames(lapply(names(ff), function(fnam){
  print(fnam)
  fx <- ff[[fnam]]
  print(colnames(fx))
  print("------------")
  fx <- fx[,!grepl("_score", colnames(fx)) & !grepl(".Norm", colnames(fx)) & !grepl("_norm", colnames(fx)), with=F]
  fx <- fx[,!colnames(fx) %in% c("n", "l"), with=F]
  names(fx) <- cleanNames(paste(names(fx), fnam))
  names(fx)[grepl("^V1", names(fx))] <- "V1"
  names(fx)[grepl("^ID", names(fx))] <- "V1"
  names(fx)[grepl("^sgRNA", names(fx))] <- "V1"
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

#lapply(data.table(do.call(rbind,strsplit(gsub(" .+", "", grep("^V1$", do.call(c, lapply(dDT, colnames)), invert = TRUE, value = TRUE)), "_"))), unique)



# Correct Cas9 and WT assignments -----------------------------------------------------------
x <- dDT[["LibA_Mye_Nov2019.txt"]]
colnames(x) <- gsub("^Cas9", "XXX", colnames(x))
colnames(x) <- gsub("^WT", "Cas9", colnames(x))
colnames(x) <- gsub("^XXX", "WT", colnames(x))
dDT[["LibA_Mye_Nov2019.txt"]] <- x

# x <- dDT[["LibTF1_GMP.Mye_July2020.txt"]]
# colnames(x) <- gsub("^CAS9", "XXX", colnames(x))
# colnames(x) <- gsub("^WT", "Cas9", colnames(x))
# colnames(x) <- gsub("^XXX", "WT", colnames(x))
# dDT[["LibTF1_GMP.Mye_July2020.txt"]] <- x


# Extract guides into one matrix ----------------------------------------------------------
lapply(dDT, function(dt) head(dt$V1))
lapply(dDT, function(dt) unique(dt[grepl("NonTargetingControlGuideForMouse", V1)]$V1))
dt <- copy(dDT[[1]])
dDT2 <-lapply(dDT, function(dt){
  # Cleaning up guide names
  dt[grep("^NonTargetingControlGuideForMouse", V1), V1 := gsub("^(NonTargetingControlGuideForMouse_\\d+)_.+$", "\\1", V1)]
  dt[!grep("^NonTargetingControlGuideForMouse", V1),V1 := gsub("(\\d+)_(\\d{6})$", "\\1", V1)]
  dt[!grep("^NonTargetingControlGuideForMouse", V1),V1 := gsub("^(.+?)_(.+?)_(\\d+)$", "\\1_\\3", V1)]
  dt[,V1 := gsub("MGLib", "", V1)]
  dt <- dt[!V1 %in% c("unmapped", "unmaped")]
  dt
})
#& !grepl("NonTargetingControl", V1),
# Merge DTs
dMT <- dDT2[[1]]
for(i in 2:length(dDT2)){
  dMT <- merge(dMT, dDT2[[i]], by="V1", all=TRUE)
}
# Convert to matrix
m <- as.matrix(dMT[,-"V1",with=F])
row.names(m) <- dMT$V1
unique(gsub("_\\d+$", "", row.names(m)))
unique(gsub("^.+_(.+)_\\d+$", "\\1", row.names(m)))

# Clean up DM names
names.to.clean <- colnames(m)
names.to.clean <- gsub("^Cas9DM_", "Cas9_DM.", names.to.clean, ignore.case = TRUE)
names.to.clean <- gsub("^Cas9_DM_", "Cas9_DM.", names.to.clean, ignore.case = TRUE)
names.to.clean <- gsub("^Cas9_DM\\.LSC_CD11b_", "Cas9_DM.LSC.CD11b_", names.to.clean, ignore.case = TRUE)
# Fix colnames for specific samples
idx <- grep("LibR1_DM_Oct2021", names.to.clean)
names.to.clean[idx] <- gsub("^DM_R1_(d\\d+) ", "Cas9_DM.\\1_R1_Oct2021 ", names.to.clean[idx])
colnames(m) <- names.to.clean



# Read Annotation --------------------------------------------------------------
# First those that do not start with "Lib" --> actual data
names1 <- grep("^Lib", colnames(m),invert = TRUE, value=TRUE)
# names1 <- gsub("^Cas9(.*?)_", "Cas9_\\1_", names1, ignore.case = TRUE)
# names1 <- gsub("^WT(.*?)_", "WT_\\1_", names1, ignore.case = TRUE)
ann <- data.table(sample=names1, do.call(rbind, strsplit(gsub("Wt", "WT", gsub(" ", "_", gsub("\\-", "_", names1))), "_")))
#ann[1:100]
ann[V4 == "Jul2020", VX := V2]
ann[V4 == "Jul2020", V2 := V3]
ann[V4 == "Jul2020", V3 := VX]
ann$VX <- NULL
ann[,V1 := gsub("CAS9", "Cas9", V1, ignore.case = TRUE)]
ann[,V1 := gsub("Cas9.+", "Cas9", V1, ignore.case = TRUE)]
table(ann$V1)
table(ann$V2)
table(ann$V3)
table(ann$V4)
table(ann$V5)

unique(ann[,c("V3", "V5"),with=F])
ann[,V5 := gsub("Lib", "", V5)]

names2 <- grep("^Lib", colnames(m),invert = FALSE, value=TRUE)
ann2 <- data.table(sample=names2, do.call(rbind, strsplit(gsub("_", "", gsub("_1 ", " ", gsub(".txt", "", names2))), " ")))
names(ann2)[3] <- "V5"
ann2$V1 <- "Library"
ann <- rbind(ann, ann2, fill=TRUE)


# Cleanup -----------------------------------------------------------------

# Add those or not?
# names2 <- grep("^Lib", colnames(m),invert = FALSE, value=TRUE)
# ann <- data.table(sample=names1, do.call(rbind, strsplit(gsub("_.+? ", " ", gsub("Lib", "", gsub(".txt", "", names2))), " ")))
# ann[V4 == "Jul2020", VX := V2]
# ann[V4 == "Jul2020", V2 := V3]
# ann[V4 == "Jul2020", V3 := VX]
# ann$VX <- NULL
# lapply(ann, unique)

# REMOVE v2 samples (aggregated by David thus duplicated)
ann[grep("v2.t..$", V7)]
ann <- ann[!grep("v2.t..$", V7)]

ann.col = data.frame(row.names=ann$sample, ann[,c("V1", "V2", "V5", "V6"), with=F])

# Check that names are unique (?)
x <- unique(data.table(apply(!is.na(m), 2, sum), gsub("^.+? ", "", colnames(m))))
x$V3 <- sapply(dDT2, nrow)[x$V2]
stopifnot(nrow(x[V1 != V3]) == 0)



# Limit matrix to proper samples ------------------------------------------
m <- m[,ann$sample]

# Correlation and guide overlaps ------------------------------------------
libs <- split(ann$sample, ann$V5)
lnam <- "BBr"
for(lnam in names(libs)){
  m2 <- m[,libs[[lnam]]]
  
  cleanDev(); pdf(out("AllData_",lnam,".pdf"), w=12, h=50)
  pheatmap(m2[apply(!is.na(m2), 1, sum) > 0,], cluster_cols = F, cluster_rows=F, fontsize_row = 5)
  dev.off()
  
  cleanDev(); pdf(out("Correlation_",lnam,".pdf"), w=12, h=12)
  #cMT <- corS(m[,ann[V5 !="B" & !grepl("\\d{4}19", V4)]$sample], use="pairwise.complete.obs")
  cMT <- corS(m2[!grepl("NonTargetingControl", row.names(m2)),], use="pairwise.complete.obs")
  diag(cMT) <- NA
  pheatmap(cMT, cluster_rows = F, cluster_cols = F, annotation_col = ann.col,
           breaks=seq(-1,1, 0.01), color=HM.COLORS.FUNC(200)
           )
  dev.off()
}



# 
# 
cleanDev(); pdf(out("Jaccard.pdf"), w=21, h=20)
pheatmap(jaccard(lapply(data.table(m), function(x) row.names(m)[!is.na(x)])),
         cluster_rows = F, cluster_cols = F, annotation_col = ann.col,
         breaks=seq(-1,1, 0.01), color=HM.COLORS.FUNC(200),
         main="Guide overlaps, jaccard coefficient")
dev.off()

# 
# 
# x <- melt(data.table(m, keep.rownames = TRUE), id.vars = "rn")
# x <- merge(x, ann, by.x="variable", by.y="sample")
# 
# cleanDev(); pdf(out("AllValues.pdf"), w=25, h=50)
# pheatmap(m, cluster_rows = F, cluster_cols = F, fontsize_row = 2)
# dev.off()


# Check different guides --------------------------------------------------
rn <- row.names(m)
rn[grepl("^NonTargetingControlGuideForMouse_\\d+$", rn)] <- paste0(rn[grepl("^NonTargetingControlGuideForMouse_\\d+$", rn)], "_X_0")
rn <- gsub("NonTargetingControlGuideForMouse_(\\d+)_(.+?)_\\d+$", "NTC_\\2", rn)
unique(gsub("^.+?_(.+?)_.+$", "\\1", rn))

# grp.cnts <- sapply(split(row.names(m), gsub("^.+?_.+$", "\\1", rn)), function(rows){
#   apply(!is.na(m[rows,]), 2, sum)
# })
# cleanDev(); pdf(out("GroupCounts.pdf"), w=10, h=25)
# pheatmap(grp.cnts)
# dev.off()


# Finalize and export data ------------------------------------------------
ann$sample <- make.names(ann$sample)
colnames(m) <- make.names(colnames(m))
stopifnot(all(ann$sample == colnames(m)))
ann <- setNames(ann, c("sample", "Genotype", "Population", "Library2", "Date", "Library", "System", "Date2"))
ann[,Date2 := gsub("\\.txt", "", Date2)]
ann[Date2 == "Feb2021" & Library == "R2", Date := "10022021"]

# Clean up system for DM
ann[System == "DM" & grepl("^CFSE", Population), System := "DM.CFSE"]
ann[System == "DM" & grepl("^CD34", Population), System := "DM.CD34"]

# Rename Non-targeting controls
row.names(m) <- gsub("NonTargetingControlGuideForMouse_", "NTC",  row.names(m))

# Export data
write.table(m[,ann$sample], quote = F, sep = ",", row.names = TRUE, col.names = TRUE, file = out("Matrix.csv"))
write.tsv(ann, out("Annotation.tsv"))
