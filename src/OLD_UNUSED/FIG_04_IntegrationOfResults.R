source("src/00_init.R")
inDirPerturb <- "FIG_02_scRNA_UMAPs/"
inDirFACS <- "FIG_01_POOLED_vsWT/"
out <- dirout("FIG_04_Integration")



# Load data ---------------------------------------------------------
perturb=list()

# in vivo
perturb[["in.vivo"]] <- fread(dirout_load(inDirPerturb)("cluster.enrichments/Cluster_enrichments_broad_in.vivo_noMixscape_14d",".tsv"))

# ex vivo
fish.enrich <- list(
  day7=fread(dirout_load(inDirPerturb)("cluster.enrichments/Cluster_enrichments_broad_ex.vivo_noMixscape_7d.tsv")),
  day9=fread(dirout_load(inDirPerturb)("cluster.enrichments/Cluster_enrichments_broad_ex.vivo_noMixscape_9d.tsv"))
)
fish.enrich <- rbindlist(fish.enrich, idcol="day")
perturb[["ex.vivo"]] = fish.enrich[day == "day7"]

# leukemia
#ann <- annList[tissue == "leukemia"]
fish.enrich <- fread(dirout_load(inDirPerturb)("cluster.enrichments/Cluster_enrichments_broad_leukemia_noMixscape_6d.tsv"))
#fish.enrich$Clusters <- unique(ann[,c("Cluster.number", "Clusters"),with=F])[match(fish.enrich$Clusters, paste("cl", Cluster.number))]$Clusters
perturb[["leukemia"]] = fish.enrich


perturb <- rbindlist(perturb, idcol = "system", fill = TRUE)
for(typex in c("HSC", "Mye", "Ery", "Lymphoid", "Baso")){
  perturb[Clusters %in% CLEAN.CELLTYPES[Type == typex]$Name, Clusters := typex]
}
perturb[Clusters %in% c("Mega", "Ery", "MkP"), Clusters := "Mega/Ery"]
unique(perturb$Clusters)
unique(perturb$sample)
perturb[, quantile(sig.perc), by=c("system", "sample")]
perturb.summarized <- perturb[sig.perc > 0.5][,.(log2OR = max(log2OR_cap)), by=c("system", "sample", "Clusters", "gene")]

ggplot(perturb.summarized[log2OR > 0], aes(x=system, y=gene, size=pmax(log2OR, 0))) + 
  theme_bw(12) +
  geom_point() + 
  facet_grid(. ~ Clusters, space = "free", scales = "free") + 
  xRot()
ggsave(out("Attempt1.pdf"),w=10,h=10)

pDT <- perturb[!Clusters %in% c("unclear", "B-cell")][system != "leukemia"]
pDT <- pDT[,.(log2OR = max(log2OR_cap), sig.perc=max(sig.perc)), by=c("system", "sample", "Clusters", "gene")]
pDT[,log2OR_cap := sign(log2OR) * pmin(abs(log2OR), 2)]
pDT.select <- dcast.data.table(pDT, gene + Clusters ~ system, value.var = "log2OR")
gg <- unique(pDT.select[order(-abs(ex.vivo - in.vivo))]$gene)[1:10]

ggplot(pDT[gene %in% gg], aes(x=gsub("\\.", " ", system), y=Clusters, color=log2OR_cap, size=sig.perc)) + 
  themeNF(rotate = TRUE, grid = FALSE) +
  geom_point() + 
  facet_grid(. ~ gene, space = "free", scales = "free") + 
  scale_color_gradient2(name="log2(odds ratio)", low="#1f78b4", high="#e31a1c") +
  scale_size_continuous(name="% significant") +
  ylab("Cell type") + xlab("Experimental system")
ggsaveNF(out("Summarized.pdf"),w=3,h=0.7, guides = TRUE)
ggsaveNF(out("Summarized_guide.pdf"),w=3,h=2, guides = TRUE)

