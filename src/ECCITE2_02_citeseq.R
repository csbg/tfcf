source(paste0(Sys.getenv("CODE"), "src/00_init.R"))

library(Matrix)
out <- dirout("ECCITE4_02_Guides/")
require(data.table)


source("~/code/resources/RFunctions/scRNA_Basics.R")

#all.data <- SCRNA.read_10Xh5.610("~/GFS/PROJECTS/TfCf/Data/ECCITE3_low_7d/outs/filtered_feature_bc_matrix.h5")
all.data <- SCRNA.read_10Xh5.610("~/GFS/PROJECTS/TfCf/Data/ECCITE2/outs/filtered_feature_bc_matrix.h5")

str(all.data)
gDT <- all.data$features[feature_type != "Gene Expression"]
gMT <- all.data$matrix[gDT$id,]

# Export counts -----------------------------------------------------------
mat <- gMT
res <- melt(data.table(as.matrix(mat), keep.rownames = TRUE), id.vars = "rn", variable.name = "barcode")
res[,barcode := gsub("\\-.+", "", barcode)]
res <- res[,.(number_of_barcodes_with_UMIs = .N), by=c("rn", "value")]
colnames(res) <- c("TAG", "UMIs", "number_of_barcodes_with_UMIs")
res <- res[UMIs != 0][order(TAG, UMIs)]
write.table(res, sep="\t", quote=F, out("citeseq.tsv"), row.names = FALSE)
res[UMIs > 5][,sum(number_of_barcodes_with_UMIs), by="TAG"]
sum(res[UMIs > 5][,sum(number_of_barcodes_with_UMIs), by="TAG"]$V1)

# Export mapping of cells to guides ----------------------------------------------------------
cutoff <- 5
mat2 <- mat
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

sort(table(guideDT$guide))
write.table(guideDT, out("citeseq_guides2cells.csv"), sep=",", quote=FALSE, row.names = FALSE)
write.table(guideDT[!grepl(" ", guide)], out("citeseq_guides2cells_withoutCombinations.csv"), sep=",", quote=FALSE, row.names = FALSE)
