# Analyze citeseq output

# install.packages('Seurat')
library(Seurat)
library(Matrix)
require(data.table)

setwd("/Volumes/GFS_MBIO_AGFORTELNY/DATA_David/Raw_data_ECCITE/citeseq/umi_count/")

# barcodes <- fread("barcodes.tsv.gz", header = F)$V1
# features <- fread("features.tsv.gz", header = F)$V1
# mtx <- fread("matrix.mtx.gz")


matrix_dir = ""
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

apply(mat[1:6,], 1, table)


res <- data.table(melt(as.matrix(mat[1:6,])))
res <- res[,.(number_of_barcodes_with_UMIs = .N), by=c("Var1", "value")]
colnames(res) <- c("TAG", "UMIs", "number_of_barcodes_with_UMIs")
res <- res[UMIs != 0][order(TAG, UMIs)]
write.table(res, sep="\t", quote=F, "/Volumes/GFS_MBIO_AGFORTELNY/DATA_David/Raw_data_ECCITE/citeseq.tsv", row.names = FALSE)

res[UMIs > 10][,sum(number_of_barcodes_with_UMIs), by="TAG"]