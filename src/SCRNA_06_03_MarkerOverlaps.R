source("src/00_init.R")
out <- dirout("SCRNA_06_03_MarkerOverlaps/")

mm <- fread(dirout_load("SCRNA_06_01_Markers")("Markers.used.tsv"))

# define datasets
query <- mm[grepl("Bulk", db)]
ref <- mm[!grepl("Bulk", db)]
bg <- unique(mm$value)

# Setup enrichment test
ref.dbs <- unique(ref$db)
geneSets <- setNames(lapply(ref.dbs, function(dbx){
  with(ref[db == dbx], split(value, L1))
}), ref.dbs)
gene.list <- with(query, split(value, L1))

# Run enrichment
fish <- fisher.test.enrichment(geneSets = geneSets, gene.list = gene.list, bg = bg, cores = 10)

# modify output
fish[, log2OR := log2(oddsRatio)]
fish[, log2OR_cap := pmin(5, abs(log2OR)) * sign(log2OR)]
fish[, padj := p.adjust(pval, method="BH")]

# plot
ggplot(fish, aes(x=list, y=geneset, color=log2OR_cap, size=pmin(5, -log10(padj)))) + 
  theme_bw(12) +
  scale_color_gradient2(low="blue", high="red") +
  geom_point() +
  facet_grid(database ~ ., space = "free", scales = "free")
ggsave(out("Enrichments.pdf"),w=8, h=15)
write.tsv(fish, out("Enrichments.tsv"))
