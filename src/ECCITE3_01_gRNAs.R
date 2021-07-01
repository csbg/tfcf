source(paste0(Sys.getenv("CODE"), "src/00_init.R"))

library(Matrix)
out <- dirout("ECCITE3_01_Guides/")
require(data.table)


source("~/code/resources/RFunctions/scRNA_Basics.R")

ff <- list.files("~/GFS/PROJECTS/TfCf/Data/", pattern="ECCITE3")
ff <- ff[!grepl("\\.log$", ff)]
ff <- ff[ff != "ECCITE3_OLD"]
ff

all.data <- lapply(ff, function(fx){
  print(fx)
  SCRNA.read_10Xh5.610(paste0("~/GFS/PROJECTS/TfCf/Data/",fx,"/outs/filtered_feature_bc_matrix.h5"))
})
names(all.data) <- ff
str(all.data)

datax <- all.data[[1]]
gMT <- lapply(all.data, function(datax){
  gDT <- datax$features[feature_type != "Gene Expression"]
  datax$matrix[gDT$id,]
})

gfpMT <- lapply(all.data, function(datax){
  gDT <- datax$features[feature_type == "Gene Expression"][id %in% c("GFP", "BFP")]
  datax$matrix[gDT$id,]
})



# Export counts -----------------------------------------------------------
mat <- gMT[[1]]
mnam <- names(gMT)[1]
for(mnam in names(gMT)){
  mat <- gMT[[mnam]]
  res <- melt(data.table(as.matrix(mat), keep.rownames = TRUE), id.vars = "rn", variable.name = "barcode")
  res[,barcode := gsub("\\-.+", "", barcode)]
  res <- res[,.(number_of_barcodes_with_UMIs = .N), by=c("rn", "value")]
  colnames(res) <- c("TAG", "UMIs", "number_of_barcodes_with_UMIs")
  res <- res[UMIs != 0][order(TAG, UMIs)]
  write.table(res, sep="\t", quote=F, out("guides_",mnam,".tsv"), row.names = FALSE)
  print(res[UMIs > 5][,sum(number_of_barcodes_with_UMIs), by="TAG"])
  print(sum(res[UMIs > 5][,sum(number_of_barcodes_with_UMIs), by="TAG"]$V1))
  
  mat <- gfpMT[[mnam]]
  res <- melt(data.table(as.matrix(mat), keep.rownames = TRUE), id.vars = "rn", variable.name = "barcode")
  res[,barcode := gsub("\\-.+", "", barcode)]
  res <- res[,.(number_of_barcodes_with_UMIs = .N), by=c("rn", "value")]
  colnames(res) <- c("TAG", "UMIs", "number_of_barcodes_with_UMIs")
  res <- res[UMIs != 0][order(TAG, UMIs)]
  write.table(res, sep="\t", quote=F, out("gfp_",mnam,".tsv"), row.names = FALSE)
}



# Export mapping of cells to guides ----------------------------------------------------------
cutoff <- 1

mnam <- names(gMT)[1]
for(mnam in names(gMT)){
  print(mnam)
  mat2 <- gMT[[mnam]]
  print(dim(mat2))
  
  mat2 <- as.matrix(mat2[,apply(mat2, 2, max) >= cutoff]) # only keep cells with any guide counted above cutoff
  ii <- apply(mat2, 2, function(col) which(col >= cutoff)) # for each cell (column) get the rows that are above the cutoff
  ii[1:10]
  table(sapply(ii, length)) # a few cells have more than one guide
  
  # Data table to be exported
  guideDT <- data.table(
    barcode=names(ii),
    guide=sapply(ii, function(x) paste(row.names(mat2)[sort(x)], collapse=" "))
  )
  
  # Cells per (combination of) guides
  guideDT[,.N, by="guide"]
  (triplets <- guideDT[sapply(ii, length) > 2])
  mat2[,triplets$barcode] # seems reliable
  
  print(guideDT[,.N, by="guide"][order(guide)])
  write.table(guideDT, out("guides2cells_",mnam,".csv"), sep=",", quote=FALSE, row.names = FALSE)
  write.table(guideDT[!grepl(" ", guide)], out("guides2cells_withoutCombinations_",mnam,".csv"), sep=",", quote=FALSE, row.names = FALSE)
}



# Export mapping of cells to gfp ----------------------------------------------------------
cutoff <- 1

mnam <- names(gMT)[1]
for(mnam in names(gMT)){
  print(mnam)
  mat2 <- gfpMT[[mnam]]
  print(dim(mat2))
  mat2 <- as.matrix(mat2[,apply(mat2, 2, max) >= cutoff]) # only keep cells with any guide counted above cutoff
  ii <- apply(mat2, 2, function(col) which(col >= cutoff)) # for each cell (column) get the rows that are above the cutoff
  ii[1:10]
  table(sapply(ii, length)) # a few cells have more than one guide
  
  # Data table to be exported
  guideDT <- data.table(
    barcode=names(ii),
    guide=sapply(ii, function(x) paste(row.names(mat2)[sort(x)], collapse=" "))
  )
  
  # Cells per (combination of) guides
  guideDT[,.N, by="guide"]
  (triplets <- guideDT[sapply(ii, length) > 2])
  mat2[,triplets$barcode] # seems reliable
  
  print(guideDT[,.N, by="guide"][order(guide)])
  write.table(guideDT, out("gfp2cells_",mnam,".csv"), sep=",", quote=FALSE, row.names = FALSE)
  write.table(guideDT[!grepl(" ", guide)], out("gfp2cells_withoutCombinations_",mnam,".csv"), sep=",", quote=FALSE, row.names = FALSE)
}