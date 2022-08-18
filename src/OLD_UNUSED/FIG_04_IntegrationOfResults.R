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
  day7=fread(dirout_load(inDirPerturb)("cluster.enrichments/Cluster_enrichments_basic_ex.vivo_noMixscape_7d.tsv")),
  day9=fread(dirout_load(inDirPerturb)("cluster.enrichments/Cluster_enrichments_basic_ex.vivo_noMixscape_9d.tsv"))
)
fish.enrich <- rbindlist(fish.enrich, idcol="day")
perturb[["ex.vivo"]] = fish.enrich[day == "day7"]

# leukemia
ann <- annList[tissue == "leukemia"]
fish.enrich <- fread(dirout_load(inDirPerturb)("cluster.enrichments/Cluster_enrichments_numeric_leukemia_noMixscape_6d.tsv"))
fish.enrich$Clusters <- unique(ann[,c("Cluster.number", "Clusters"),with=F])[match(fish.enrich$Clusters, paste("cl", Cluster.number))]$Clusters
perturb[["leukemia"]] = fish.enrich


perturb <- rbindlist(perturb, idcol = "system", fill = TRUE)
for(typex in c("HSC", "Mye", "Ery", "Lymphoid", "Baso")){
  perturb[Clusters %in% CLEAN.CELLTYPES[Type == typex]$Name, Clusters := typex]
}
perturb[Clusters %in% c("Mega", "Ery", "MkP"), Clusters := "Mega/Ery"]
unique(perturb$Clusters)
unique(perturb$sample)
perturb[, quantile(sig.perc), by=c("system", "sample")]
perturb <- perturb[sig.perc > 0.5][,.(log2OR = max(log2OR_cap)), by=c("system", "sample", "Clusters", "gene")]

ggplot(perturb[log2OR > 0], aes(x=system, y=gene, size=pmax(log2OR, 0))) + 
  theme_bw(12) +
  geom_point() + 
  facet_grid(. ~ Clusters, space = "free", scales = "free") + 
  xRot()
ggsave(out("Attempt1.pdf"),w=10,h=10)
