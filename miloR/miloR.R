
# SETUP -------------------------------------------------------------------
path <- "/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SCRNA_09_miloR/"
dir.create(path)
setwd(path)


# FUNCTIONS ---------------------------------------------------------------
source("miloR_test_function_NF.R")

require(miloR)
require(scater)
require(pheatmap)
require(dplyr)
require(data.table)
require(igraph)



# Create miloR Object -----------------------------------------------------
(load("/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SCRNA_02_01_Integration/in.vivo/MonocleObject.RData"))

table(monocle.obj$guide)
colnames(colData(monocle.obj))

mobj <- monocle.obj[,!is.na(monocle.obj$guide)]
traj_milo <- Milo(mobj)
reducedDim(traj_milo, "UMAP") <- reducedDim(mobj, "UMAP")
traj_milo <- buildGraph(traj_milo, k = 10, d = 30)
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d=30, refined = TRUE)
plotNhoodSizeHist(traj_milo)

traj_milo$miloR_gene <- gsub("_.+$", "", traj_milo$CRISPR_Cellranger)
traj_milo$miloR_guide <- traj_milo$CRISPR_Cellranger

traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="miloR_guide")
as.matrix(nhoodCounts(traj_milo))[1:10, 1:10]

traj_milo <- calcNhoodDistance(traj_milo, d=30)

traj_design <- data.frame(colData(traj_milo))[,c("miloR_gene", "miloR_guide"),drop=F]
traj_design <- distinct(traj_design)
traj_design$miloR_gene <- relevel(factor(traj_design$miloR_gene), ref = "NTC")
row.names(traj_design) <- traj_design$miloR_guide

da_results <- testNhoods_NF(traj_milo, design=~miloR_gene, design.df = traj_design)
saveRDS(da_results, "DF_res.RDS")
saveRDS(traj_milo, "miloR_obj.RDS")



# Plot --------------------------------------------------------------------
da_results <- readRDS("DF_res.RDS")
traj_milo <- readRDS("miloR_obj.RDS")
da_results <- data.table(da_results)

str(da_results)
table(da_results$coef)
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50) + facet_wrap(~ coef)
ggsave("P.Value_Distributions.pDF", h=30, w=30)

da_results[, coef := gsub("miloR_gene", "", coef)]

traj_milo <- buildNhoodGraph(traj_milo)


nh_graph <- nhoodGraph(traj_milo)
nh_graph <- permute(nh_graph, order(V(nh_graph)$size, decreasing = TRUE))
node.idx <- unlist(nhoodIndex(traj_milo)[unique(da_results$Nhood)])


pDT <- lapply(setdiff(unique(traj_milo$miloR_gene), "NTC"), function(gx){
  cbind(
    reducedDim(traj_milo, "UMAP")[node.idx,],
    da_results[coef==gx])
})

pDT <- bind_rows(pDT)

ggplot(pDT[coef != "(Intercept)"], aes(x=V1, V2, color=logFC, alpha=-log10(PValue))) + 
  geom_point() + 
  scale_color_gradient2() +
  theme_bw() +
  facet_wrap(~coef)
ggsave("UMAP.pdf", h=30, w=30)


sort(unique(pDT$coef))
pDT2 <- pDT[coef != "(Intercept)"][coef %in% c("Cebpa","Smarcd2","Rcor1","Wdr82","Brd9","Rbbp4","Hdac3","Chd4")]
ggplot(pDT2, aes(x=V1, V2, color=logFC, alpha=-log10(FDR), size=-log10(FDR))) + 
  geom_point() + 
  scale_color_gradient2(name=expression(log(FC)), low="#1f78b4", high="#e31a1c") +
  scale_alpha_continuous(name=expression(p[adj])) +
  scale_size_continuous(name=expression(p[adj]), range=c(0,3)) +
  theme_bw() +
  facet_wrap(~coef, ncol = 4)
ggsave("UMAP_selection.pdf", h=6, w=12)


