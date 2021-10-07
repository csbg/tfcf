inF <- "metadata/ECCITE7_Features.csv"
feat <- fread(inF)
gmap <- fread(dirout_load("02_EnsemblForGenes_10x")("GMAP.tsv"))

feat[target_gene_name == "Gltscr1", target_gene_name := "Bicra"]

feat$target_gene_id <- gmap[match(feat$target_gene_name, Gene),]$ENSG
feat[target_gene_name == "Non-targeting control", target_gene_id := "Non-targeting control"]

feat[is.na(target_gene_id)]
stopifnot(!any(is.na(feat$target_gene_id)))

write.table(feat, gsub("\\.csv", "_ensgs.csv", inF), sep=",", quote = F, col.names = TRUE, row.names = FALSE)
