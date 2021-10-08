source("src/00_init.R")

inF <- "metadata/ECCITE7_Features.csv"
feat <- fread(inF)
gmap <- fread(dirout_load("02_EnsemblForGenes_10x")("GMAP.tsv"))

# Clean file
feat[target_gene_name == "Gltscr1", target_gene_name := "Bicra"]
feat[,id := gsub("_(.+?)_", "_\\U\\1_", id, perl = TRUE)]
feat[,id := gsub("Men_", "Men1_", id, perl = TRUE)]
feat[,name := id]

feat$target_gene_id <- gmap[match(feat$target_gene_name, Gene),]$ENSG
feat[target_gene_name == "Non-targeting control", target_gene_id := "Non-Targeting"]
feat[target_gene_name == "Non-targeting control", target_gene_name := "Non-Targeting"]

feat <- unique(feat)

stopifnot(!any(is.na(feat$target_gene_id)))
# feat[is.na(target_gene_id)]

stopifnot(!any(duplicated(feat$id)))

stopifnot(!any(duplicated(feat$name)))

stopifnot(!any(duplicated(feat$sequence)))
# feat[sequence %in% feat[,.N, by=c("sequence")][N>1]$sequence][order(sequence)]

write.table(feat, gsub("\\.csv", "_ensgs.csv", inF), sep=",", quote = F, col.names = TRUE, row.names = FALSE)
