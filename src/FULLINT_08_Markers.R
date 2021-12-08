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
marker.lists[["Larry"]] <- markers


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

# Calculate signatures ----------------------------------------------------
dat <- SCRNA.TPXToLog(SCRNA.RawToTPX(counts(monocle.obj), scale.factor = 1e6))
mnam <- "PanglaoDB"
for(mnam in names(marker.lists)){
  markers <- marker.lists[[mnam]]
  dat2 <- dat[row.names(dat) %in% unique(do.call(c, markers)),]
  dat2 <- dat2[apply(dat2, 1, max) > 0,]
  dat2 <- dat2 - apply(dat2, 1, min)
  dat2 <- dat2 / apply(dat2, 1, max)
  markers <- lapply(markers, function(x) x[x %in% row.names(dat2)])
  markers <- markers[sapply(markers, length) >= 5]
  # markers <- lapply(markers, function(gg) gg[gg %in% row.names(dat2)])
  # markers <- markers[sapply(markers, length) > 0]
  # str(markers)
  sigs <- sapply(markers, function(gg) Matrix::colMeans(dat2[gg,,drop=F]))
  #ggplot(data.table(melt(sigs)), aes(x=Var2, y=value)) + geom_violin()
  
  write.table(sigs, out("Signatures_",mnam,".csv"), quote=F, row.names = TRUE, col.names = TRUE, sep=",")
}

