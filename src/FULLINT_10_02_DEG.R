source("src/00_init.R")
out <- dirout("FULLINT_10_02_DEG")

require(monocle3)
require(fgsea)
source("src/FUNC_Monocle_PLUS.R")


# Load full Monocle object and differential expression analysis -----------
(load(PATHS$FULLINT$Monocle))
(load(PATHS$FULLINT$DEG))
res.clean <- fread(PATHS$FULLINT$DEG.clean)
ann <- fread(PATHS$FULLINT$DEG.ann)
logFCMT <- as.matrix(read.csv(PATHS$FULLINT$DEG.logFCMT))
cds <- monocle.obj[,ann$rn]


# Plots some genes (only in vitro) -------------------------------------------
res.ex <- res.clean[interaction == FALSE & tissue == "in vitro" & q_value < 0.05 & estimate > 0]
res.ex[,.N, by="guide"]

genex <- "Kmt2d"

# Plot analysis for one guide
gg <- res.ex[guide == genex]$gene_id
annx <- ann[GuideDE %in% c("NTC", genex) & tissueDE == "in vitro"]

# Prepare exprsesion matrix
eMT <- NF_TPM_Matrix(cds, gg)
eMT <- eMT[,annx$rn]
eMT <- SCRNA.TPXToLog(eMT)
eMT <- t(scale(t(as.matrix(eMT))))
eMT[eMT > 5] <- 5
eMT[eMT < -5] <- -5

# column annotation
annC <- data.frame(
  row.names = annx$rn,
  guide=annx$GuideDE,
  cluster=annx$ClusterDE
)

# Heatmap
cleanDev(); pdf(out("HM.pdf"),w=10,h=10)
pheatmap(
  cluster_cols = FALSE,
  annotation_col = annC,
  eMT[,annx[order(GuideDE, ClusterDE)]$rn], 
  show_rownames = FALSE, 
  show_colnames = FALSE
  )
dev.off()

# Violin plots
str(colData(cds))
if(length(gg) > 40) gg <- sample(gg, 40)
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
ggsave(out("Violins.pdf"), w=30,h=30)



# Recalculate interaction coefficients -------------------------------------------
res.invitro <- res.clean[interaction == FALSE & tissue == "in vitro"]
res.interaction.leukemia <- res.clean[interaction == TRUE & tissue == "leukemia"]
res.leukemia <- fread(dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")("DEG_Results.tsv"))
gg <- intersect(intersect(res.invitro$guide, res.interaction.leukemia$guide), res.leukemia$guide)

pDT <- merge(res.invitro, res.leukemia, by=c("guide", "gene_id"), suffixes=c("_vitro", "_leuk"))
pDT <- merge(pDT, res.interaction.leukemia, by=c("guide", "gene_id"))
ggplot(pDT, aes(x=estimate_leuk, y=estimate_vitro + estimate)) + 
  theme_bw(12) +
  facet_wrap(~guide) +
  geom_hex() +
  xlab("Measured estimate (leukemia)") +
  ylab("Calculated estimate (vitro + interaction)")
ggsave(out("InteractionCoefficients.pdf"), w=10,h=10)

ggplot(pDT[,cor(estimate_leuk, estimate_vitro + estimate), by="guide"],
       aes(x=guide, y=V1)) +
  theme_bw(12) +
  geom_bar(stat='identity') +
  xRot()
ggsave(out("InteractionCoefficients_Corr.pdf"), w=6,h=4)




# Enrichments -------------------------------------------------------------
fgsea
