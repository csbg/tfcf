# Analyze citeseq output

# install.packages('Seurat')
library(Seurat)
library(Matrix)
require(data.table)

setwd("/Volumes/GFS_MBIO_AGFORTELNY/PROJECTS/TfCf/")

# barcodes <- fread("barcodes.tsv.gz", header = F)$V1
# features <- fread("features.tsv.gz", header = F)$V1
# mtx <- fread("matrix.mtx.gz")


matrix_dir = "citeseq_all_barcodes/umi_count/"
matrix_dir = "citeseq_combined//umi_count/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

apply(mat[1:6,], 1, table)



# Export counts -----------------------------------------------------------
res <- data.table(melt(as.matrix(mat[1:6,])))
res[,Var1 := gsub("\\-.+", "", Var1)]
res <- res[,.(number_of_barcodes_with_UMIs = .N), by=c("Var1", "value")]
colnames(res) <- c("TAG", "UMIs", "number_of_barcodes_with_UMIs")
res <- res[UMIs != 0][order(TAG, UMIs)]
write.table(res, sep="\t", quote=F, "citeseq.tsv", row.names = FALSE)
res[UMIs > 5][,sum(number_of_barcodes_with_UMIs), by="TAG"]
sum(res[UMIs > 5][,sum(number_of_barcodes_with_UMIs), by="TAG"]$V1)

# Export mapping of cells to guides ----------------------------------------------------------
cutoff <- 5
mat2 <- mat
row.names(mat2) <- gsub("\\-.+", "", row.names(mat2))
colnames(mat2) <- paste0(colnames(mat), "-1") # To fit to the barcodes from cellranger
mat2 <- mat2[row.names(mat2) != 'unmapped',]
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

# compare to the barcodes from the RNA
barcodes <- fread("RNA_cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", header = F)$V1
#barcodes <- fread("~/Downloads/citeseq_guides2cells.csv")$barcode

# for 5 cells we have guides assigned but no transcriptome
guideDT[!barcode %in% barcodes]
guideDT <- guideDT[barcode %in% barcodes] 

guideDT <- rbind(
  guideDT,
  data.table(barcode=setdiff(barcodes, guideDT$barcode), guide="None")
)

write.table(guideDT, "citeseq_guides2cells2.csv", sep=",", quote=FALSE, row.names = FALSE)
