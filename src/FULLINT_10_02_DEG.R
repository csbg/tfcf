source("src/00_init.R")
out <- dirout("FULLINT_10_02_DEG")

require(monocle3)

(load(PATHS$FULLINT$Monocle))
(load(PATHS$FULLINT$DEG))
res.clean <- fread(PATHS$FULLINT$DEG.clean)
ann <- fread(PATHS$FULLINT$DEG.ann)

cds <- monocle.obj[,ann$rn]

res.ex <- res.clean[interaction == FALSE & tissue == "in vitro" & q_value < 0.05 & estimate > 0]
res.ex[,.N, by="guide"]

genex <- "Kmt2d"
gg <- res.ex[guide == genex]$gene_id

annx <- ann[GuideDE %in% c("NTC", genex) & tissueDE == "in vitro"]
eMT <- NF_TPM_Matrix(cds, gg)
eMT <- eMT[,annx$rn]

eMT <- SCRNA.TPXToLog(eMT)
eMT <- t(scale(t(as.matrix(eMT))))
eMT[eMT > 5] <- 5
eMT[eMT < -5] <- -5

annC <- data.frame(
  row.names = annx$rn,
  guide=annx$GuideDE,
  cluster=annx$ClusterDE
)

cleanDev(); pdf(out("HM.pdf"),w=10,h=10)
pheatmap(
  cluster_cols = FALSE,
  annotation_col = annC,
  eMT[,annx[order(GuideDE, ClusterDE)]$rn], 
  show_rownames = FALSE, 
  show_colnames = FALSE
  )
dev.off()

str(colData(cds))
pObj <- cds[gg,annx$rn]
pObj$cluster <- annx$ClusterDE
pObj$group <- paste(pObj$guide, pObj$tissue, paste0("c", pObj$cluster))
plot_genes_violin(
  pObj,
  log_scale = FALSE,
  label_by_short_name=FALSE,
  group_cells_by="group", 
  ncol = 10) + 
  scale_y_continuous(trans = "log1p") +
  xRot()
ggsave(out("Violins.pdf"), w=20,h=10)
