source("src/00_init.R")

base.dir <- "SCRNA_60_01_Leukemia_Kmt2d/"
out <- dirout(base.dir)

ann.kos <- fread(dirout_load("SCRNA_20_Summary/leukemia_monocle")("Annotation.tsv"))
ann.kos2 <- fread(dirout_load("SCRNA_20_Summary/leukemia_monocle.singleR")("Annotation.tsv"))
stopifnot(all(ann.kos2$rn == ann.kos$rn))

ann.kos$celltype <- ann.kos2$Clusters
ann.kos[, gene := gsub("_.+$", "", CRISPR_Cellranger)]


# Summarize by cluster and celltype ---------------------------------------
xDT <- ann.kos[gene %in% c("Kmt2d", "NTC")]
xDT <- xDT[,.N, by=c("Clusters", "gene", "celltype")]
xDT <- dcast.data.table(xDT, Clusters + celltype ~ gene, value.var="N")
xDT[is.na(Kmt2d), Kmt2d := 0]
xDT[is.na(NTC), NTC := 0]
xDT[, fraction := Kmt2d / NTC]
#xDT[is.na(fraction), fraction := 0]

# Numbers in clusters
ggplot(xDT, aes(x=factor(as.numeric(Clusters)), y=fraction, size=NTC)) + 
  theme_bw(12) +
  geom_point() + 
  scale_size_continuous(range = c(1,5)) +
  facet_grid(.~ celltype, space="free", scale="free") +
  xRot() +
  xlab("Cluster")
ggsave(out("Numeric_Clusters.pdf"), w=8,h=4)



# Summary by cluster for UMAP ---------------------------------------------
xDT <- ann.kos[gene %in% c("Kmt2d", "NTC")]
xDT <- xDT[,.N, by=c("Clusters", "gene")]
xDT <- dcast.data.table(xDT, Clusters ~ gene, value.var="N")
xDT[is.na(Kmt2d), Kmt2d := 0]
xDT[is.na(NTC), NTC := 0]
xDT[, fraction := Kmt2d / NTC]
xDT <- merge(xDT, ann.kos[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by="Clusters"], by="Clusters")

# Numbers in clusters
ggplot(xDT, aes(x=UMAP1, y=UMAP2)) + 
  geom_hex(data=ann.kos, fill="lightgrey", bins=50) +
  theme_bw(12) +
  geom_point(aes(size=NTC, color=fraction)) + 
  scale_size_continuous(name = "Nr of NTCs", range = c(1,5)) +
  geom_text_repel(aes(label=Clusters), color="black") + 
  scale_color_gradientn(name = "Kmt2d/NTC", colors=c("black", "purple", "red", "orange")) + 
  xRot()
ggsave(out("Clusters_UMAP.pdf"), w=5,h=4)
