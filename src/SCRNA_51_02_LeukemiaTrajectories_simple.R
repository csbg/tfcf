source("src/00_init.R")

inp.dir <- "SCRNA_51_01_LeukemiaTrajectories/"
base.dir <- "SCRNA_51_02_LeukemiaTrajectories_simple/"
out <- dirout(base.dir)
inDir <- dirout_load(inp.dir)


# Samples annotation --------------------------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)

# Read leukemia data -------------------
mobjs.ann <- readRDS(dirout_load("SCRNA_06_02_MergeMarkers")("CellTypes_", "leukemia", ".RDS"))
mobjs.umap <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle.RDS"))



# Read prediction data ----------------------------------------------------
in.acc <- fread(inDir("Accuracy.tsv"))
in.pred <- fread(inDir("Predictions_prob.tsv"))


# ct centroids ------------------------------------------------------------
ct.centroids <- merge(in.acc, mobjs.umap, by="rn")[, .(UMAP_1 = median(UMAP_1), UMAP_2=median(UMAP_2)), by=c("celltype")]
centMT <- t(ct.centroids[,-"celltype",with=F])
colnames(centMT) <- ct.centroids$celltype

ctx <- "Eo/Ba"
trajDT <- data.table()
for(ctx in setdiff(colnames(centMT), "LSC")){
  M <- centMT[,ctx] - centMT[,"LSC"]
  rns <- in.pred[cell.type == ctx]$rn
  P = as.matrix(mobjs.umap[match(rns, rn)][,c("UMAP_1", "UMAP_2"),with=F])
  t0 = (t(t(P) - centMT[,ctx]) %*% M) / (M %*% M)[1]
  trajDT <- rbind(trajDT, data.table(mobjs.umap[match(rns, rn)], traj=setNames(t0[,1], rns), cell.type = ctx))
}

ggplot(trajDT, aes(x=UMAP_1, y=UMAP_2)) + 
  theme_bw(12) +
  geom_hex(data=mobjs.umap[tissue == "leukemia"], fill="lightgrey", color=NA) +
  geom_point(aes(color=traj), size=0.3) +
  scale_color_gradient2(mid="#fdbf6f") +
  facet_grid(. ~ cell.type) +
  geom_point(data=ct.centroids)

table(in.pred$cell.type)


# test KOs ---------------------------------------------------------------
ann.kos <- fread(dirout_load("SCRNA_20_Summary/leukemia_monocle.singleR")("Annotation.tsv"))
koDT <- merge(ann.kos, fread(out("Predictions_prob.tsv")), by="rn")

pDTx <- koDT[!is.na(CRISPR_Cellranger)]
pDTx[, celltype := cell.type]
pDTx[, gene := gsub("_.+$", "", CRISPR_Cellranger)]
pDTx[, traj := cell.type.prob]


# . test ------------------------------------------------------------------
typex <- "Ery"
gx <- "Rcor1"
res <- data.table()
for(typex in unique(pDTx$celltype)){
  pDT1 <- pDTx[celltype == typex]
  for(gx in unique(pDT1[mixscape_class.global != "NTC"]$gene)){
    x1 <- pDT1[gene == gx]$traj
    x2 <- pDT1[gene == "NTC"]$traj
    if(length(x1) > 10 & length(x2) > 10){
      res <- rbind(res, data.table(
        p.wx=wilcox.test(x1, x2)$p.value,
        p.ks=ks.test(x1, x2)$p.value,
        d=median(x1) - median(x2),
        type=typex,
        gene=gx
      ))
    }
  }
}
res[, padj.wx := p.adjust(p.wx, method="BH")]
res[, padj.ks := p.adjust(p.ks, method="BH")]
write.tsv(res, out("Statistics.tsv"))


# . load ------------------------------------------------------------------
#res <- fread(outS("Statistics.tsv"))

# . plot stats ------------------------------------------------------------
ggplot(res, aes(y=-log10(p.wx+1e-10), x=-log10(p.ks+1e-10), color=d)) + 
  geom_point() +
  geom_point(shape=1, color="black") +
  theme_bw(12)+ 
  geom_text_repel(aes(label=paste(gene)), color="black")+
  scale_color_gradient2(low="blue", high="red") +
  facet_grid(. ~ type)
ggsave(out("Statistics_Comparison.pdf"),w=20,h=6)

# xDT <- melt(res, id.vars = c("type", "gene"))
# xDT[, measurement := gsub("\\..+$", "", variable)]
# xDT[, type := gsub("^.+?\\.", "", variable)]
ggplot(res, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
  geom_point() +
  theme_bw(12)+ 
  scale_color_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("Statistics.pdf"), w=4,h=10)

# . plot distributions ----------------------------------------------------
pDT.distr <- copy(pDTx)
pDT.sum <- pDT.distr[, .(traj = median(traj), q1 = quantile(traj, 0.25), q2 = quantile(traj, 0.75)),by=c("gene", "celltype")]
pDT.stats <- copy(res)
pDT.stats[, celltype := type]
pDT.stats[, type := "not.sig"]
pDT.stats[padj.ks < 0.1, type := "sig.low"]
pDT.stats[padj.ks < 0.01, type := "sig.high"]
pDT.distr <- merge(pDT.distr, pDT.stats[,c("gene", "celltype", "type"),with=F], by=c("gene", "celltype"), all.x=TRUE)
pDT.distr[gene == "NTC", type := "NTC"]
ggplot(pDT.distr, aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum, color="black") + 
  geom_errorbarh(data=pDT.sum, color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ celltype, scale="free_x") +
  xRot()
ggsave(out("Statistics_Distribution.pdf"), w=20,h=10)
