source("src/00_init.R")
out <- dirout("FULLINT_08_01_Markers/")


# LOAD DATA ---------------------------------------------------------------
# Human/mouse map
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)

# Our dataset
(load(PATHS$FULLINT$Monocle))


marker.lists <- list()

# Download reference data LARRY -------------------------------------------------
geo.file <- out("GEO.txt.gz")
if(!file.exists(geo.file)){
  system(paste("wget","--no-check-certificate","-O",geo.file,"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60101/suppl/GSE60101_1256271tableS2.txt.gz"))
}
refx <- fread(geo.file, skip = 1)
refx <- refx[,-"HSC"]
refx <- refx[,-c("MF", "CD4", "CD8", "NK")]
m <- as.matrix(refx[,3:ncol(refx)])
row.names(m) <- refx$NAME
# Data should be normalized to 10 M reads but I think it's not:
barplot(colSums(m), las=2)
m <- t(t(m)/colSums(m))* 1e6
barplot(colSums(m), las=2)
m <- log1p(m)

# Get markers LARRY
cx <- colnames(m)[1]
rx <- row.names(m)[1]
markers <- lapply(colnames(m), function(cx){
  m2 <- m[,-which(colnames(m) == cx)]
  delta <- m[,cx] - apply(m2, 1, max)
  unique(names(tail(sort(delta[delta > 0.5]), 100)))
})
names(markers) <- colnames(m)
write.tsv(melt(markers), out("Markers_Larry.tsv"))
marker.lists[["Larry"]] <- markers



# Larry dataset - processed by David --------------------------------------
marker.lists[["LarrySelf"]] <- with(fread("metadata/markers.LarryEtAl.david.tsv"), split(Gene, CellType))


# This dataset (bulk) -----------------------------------------------------
for(typex in c("Old", "New")){
  infile <- if(typex == "Old")"metadata/INVITRO_BulkExpression.csv" else "metadata/INVITRO_BulkExpression_NEW.csv"
  str(bulk.data <- as.matrix(read.table(infile, row.names = 1, sep=",", header = TRUE)))
  bulk.data <- sapply(split(colnames(bulk.data), gsub(".\\d+$", "", colnames(bulk.data))), function(sx){rowMeans(bulk.data[,sx,drop=F])})
  bulk.data <- bulk.data[apply(bulk.data, 1, max) > median(apply(bulk.data, 1, max)),]
  bulk.data.fc <- sapply(colnames(bulk.data), function(cx){bulk.data[,cx] / apply(bulk.data[,colnames(bulk.data) != cx], 1, max)})
  x <- lapply(apply(bulk.data.fc, 2, function(cx) which(cx > 3)), names)
  x <- lapply(x, function(ii) gsub("\\.", "-", ii))
  marker.lists[[paste0("Bulk", typex)]] <- x
}



# Enrichr ----------------------------------------------------------------
enr.file <- out("EnrichR.RData")
if(file.exists(enr.file)){
  load(enr.file)
} else {
  enr.terms <- enrichrGetGenesets(c("CellMarker_Augmented_2021", "PanglaoDB_Augmented_2021"))
  save(enr.terms, file=enr.file)
}

# Convert to mouse --------------------------------------------------------
hm <- unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
names(hm) <- c("Mouse", "Human")
enr.terms <- lapply(enr.terms, function(dbl){
  dbl <- lapply(dbl, function(gs){
    unique(hm[Human %in% gs]$Mouse)
  })
  dbl[sapply(dbl, length) > 0]
})
# THESE MARKERS ARE CURRENTLY NOT USED!
# marker.lists <- c(marker.lists, enr.terms)

# x <- colnames(mx)
# sort(table(gsub(".+\\.", "", x)))
# x <- x[grepl("Marrow$", x) | grepl("Blood$", x) | grepl("Undefined$", x)]
# x <- x[!grepl("T\\.+cell", x)]
# x <- x[!grepl("T.Helper", x)]
# x <- x[!grepl("T\\.reg", x)]

# Panglao DB --------------------------------------------------------------
panglaoDB.file <- out("PanglaoDB.txt.gz")
if(!file.exists(panglaoDB.file)){
  refx <- fread("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")
  write.tsv(refx, panglaoDB.file)
}
refx <- fread(panglaoDB.file, check.names = TRUE)
refx <- refx[ubiquitousness.index < 0.1]
refx <- refx[grepl("Mm", species)]
refx <- refx[organ %in% c("Blood", "Bone", "Immune system")]
refx <- refx[!grepl("^Osteo", official.gene.symbol)]
refx <- refx[canonical.marker == 1]
m <- with(refx, split(official.gene.symbol, cell.type))
m <- lapply(m, function(mx) unique(hm.map[Human.gene.name %in% mx]$Gene.name))
m <- m[sapply(m, length) >= 10]
marker.lists[["PanglaoDB"]] <- m


# Izzo et al., https://doi.org/10.1038/s41588-020-0595-4 -----------------------------------------------------------------
izzo.file <- out("IzzoEtAl.xls")
if(!file.exists(izzo.file)){
  url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-020-0595-4/MediaObjects/41588_2020_595_MOESM3_ESM.xlsx"
  download.file(url, izzo.file)
}
x <- readxl::read_xlsx(izzo.file,sheet = 2, skip = 2)
x <- data.table(x, check.names = TRUE)[,1:8]
x <- setNames(x, c("p", "logFC", "perc1", "perc2", "padj", "cl", "gene", "cell"))
x <- x[!is.na(p),]
x <- x[logFC > 0.5 & padj < 0.05]
marker.lists[["IzzoEtAl"]] <- with(x, split(gene, cell))



# Store markers -----------------------------------------------------------
exportDT <- do.call(rbind, lapply(names(marker.lists), function(lnam){
  data.table(melt(marker.lists[[lnam]]), db=lnam)
}))
write.tsv(exportDT, out("Markers.used.tsv"))

# Calculate signatures ----------------------------------------------------
dat <- SCRNA.TPXToLog(SCRNA.RawToTPX(counts(monocle.obj), scale.factor = 1e6))
min.reads <- ncol(dat) * 0.001
dat <- dat[Matrix::rowSums(dat) > min.reads,]

mnam <- "PanglaoDB"
mnam <- names(marker.lists)[5]
for(mnam in names(marker.lists)){
  mfile <- out("Signatures_",mnam,".csv")
  if(file.exists(mfile)) next
  
  markers <- marker.lists[[mnam]]
  markers <- lapply(markers, function(x) x[x %in% row.names(dat)])
  markers <- markers[sapply(markers, length) >= 5]
  
  dat2 <- dat[row.names(dat) %in% unique(do.call(c, markers)),]
  dat2 <- dat2 - rowMins(dat2)
  dat2 <- dat2 / rowMaxs(dat2)
  # markers <- lapply(markers, function(gg) gg[gg %in% row.names(dat2)])
  # markers <- markers[sapply(markers, length) > 0]
  # str(markers)
  sigs <- sapply(markers, function(gg) Matrix::colMeans(dat2[gg,,drop=F]))
  sigs <- round(sigs, 3)
  #ggplot(data.table(melt(sigs)), aes(x=Var2, y=value)) + geom_violin()
  
  write.table(sigs, mfile, quote=F, row.names = TRUE, col.names = TRUE, sep=",")
}

