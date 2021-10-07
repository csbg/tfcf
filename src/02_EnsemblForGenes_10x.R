source("src/00_init.R")
out <- dirout("02_EnsemblForGenes_10x/")


require(data.table)
x <- fread(paste(Sys.getenv("GFS"), "RESOURCES", "Genomes", "refdata-gex-mm10-2020-A", "genes", "genes.gtf", sep="/"))
x <- x[V3 == "gene"]

gsub("^.*?gene_id .(.+?).\\;*$", "\\1", x[1]$V9)

x[, gene_id := gsub('^.*?gene_id \\"(.+?)\\".*$', "\\1", V9)]
x[, gene_name := gsub('^.*?gene_name \\"(.+?)\\".*$', "\\1", V9)]
gmap <- unique(setNames(x[,c("gene_id", "gene_name")], c("ENSG", "Gene")))

stopifnot(all(gmap[,.N, by="ENSG"]$N == 1))
gmap[,.N, by="Gene"][N>1][order(N)]
write.tsv(gmap, out("GMAP.tsv"))
