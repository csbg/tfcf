source("src/00_init.R")

base.dir <- "SCRNA_51_01_LeukemiaTrajectories/"
out <- dirout(base.dir)

require(glmnet)
require(doMC)
require(ggrepel)
registerDoMC(cores = 10)

# Samples annotation --------------------------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)

# Read leukemia data -------------------
mobjs.ann <- readRDS(dirout_load("SCRNA_06_02_MergeMarkers")("CellTypes_", "leukemia", ".RDS"))
mobjs.umap <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle.RDS"))
(load(PATHS$SCRNA$MONOCLE.DIR("leukemia")))
mobjs <- monocle.obj



# Modify matrix -----------------------------------------------------------
mat <- counts(mobjs)
mat.log <- SCRNA.TPXToLog(SCRNA.RawToTPX(mat, 1e6))
non.zero.rows <- tabulate(mat.log@i + 1)
stopifnot(sum(mat.log[1, ] != 0) == non.zero.rows[1])
stopifnot(sum(mat.log[3, ] != 0) == non.zero.rows[3])
quantile(non.zero.rows)
mat.log <- mat.log[non.zero.rows > ncol(mat.log)/100,]
str(mat.log)
mat.log <- t(as.matrix(mat.log))

# Data for training the model ---------------------------------------------
en.ann <- mobjs.ann[labels %in% c("Eo/Ba", "LSC", "GMP (late)", "Mono")]#[sample(1:nrow(en.ann), 1000)]
en.dat <- mat.log[en.ann$cellname,]


# Model fit ---------------------------------------------------------------
fit <- cv.glmnet(x = en.dat, y = factor(en.ann$labels), family="multinomial", parallel=TRUE, keep=TRUE)


# Examing results ---------------------------------------------------------
str(fit$foldid)
str(fit$fit.preval)
fit$lambda.min
fit$lambda
plot(fit$cvm)


# Plot accuracy -----------------------------------------------------------
fit.predicted <- fit$fit.preval[, , which(fit$lambda==fit$lambda.min)]
fit.predicted <- data.table(data.frame(fit.predicted), keep.rownames = TRUE)
fit.predicted$celltype <- en.ann[match(fit.predicted$rn, cellname)]$labels
write.tsv(fit.predicted, out("Accuracy.tsv"))
pDT <- melt(fit.predicted, id.vars = c("celltype", "rn"))

# Plot valuse
ggplot(pDT, aes(x=variable, y=value)) + 
  theme_bw() +
  geom_boxplot(coef=Inf) + 
  facet_grid(. ~ celltype) +
  xRot()
ggsave(out("Accuracy.pdf"), w=8,h=5)

# Plot UMAP
pDT <- merge(pDT, mobjs.umap, by="rn")
ggplot(pDT, aes(x=UMAP_1, y=UMAP_2)) + 
  stat_summary_hex(bins = 50, aes(z=value),fun=mean) +
  facet_grid(. ~ variable) +
  scale_fill_gradient2() +
  theme_bw(12)
ggsave(out("Accuracy_UMAP.pdf"), w=20,h=5)


# Create new predictions --------------------------------------------------
pred.dat <- mat.log#[!row.names(mat.log) %in% en.ann$cellname,]
pred.result <- predict(object = fit, newx = pred.dat, type = "response")

# write
pDT <- data.table(pred.result[,,1], keep.rownames = TRUE)
write.tsv(pDT, out("Predictions.tsv"))

# Plot
pDT <- melt(pDT, id.vars = "rn")
pDT <- merge(pDT, mobjs.umap, by="rn")
ggplot(pDT, aes(x=UMAP_1, y=UMAP_2)) + 
  stat_summary_hex(bins = 50, aes(z=value),fun=mean) +
  facet_grid(. ~ variable) +
  theme_bw(12)
ggsave(out("Predictions.pdf"), w=20,h=5)




# Predict cell types ------------------------------------------------------
predDT <- fread(out("Predictions.tsv"))
predMT <- as.matrix(predDT[,-c("rn", "LSC")])
predLabels <- colnames(predMT)[apply(predMT, 1, which.max)]
predClassProb <- apply(predMT, 1, function(row) row[which.max(row)])
predDT$cell.type <- predLabels
predDT$cell.type.prob <- predClassProb
write.tsv(predDT, out("Predictions_prob.tsv"))

# Plot predicted celltype on UMAP
predDT <- merge(predDT, mobjs.umap, by="rn")
ggplot(predDT, aes(x=UMAP_1, y=UMAP_2)) + 
  geom_hex(bins = 50) +
  facet_grid(. ~ cell.type) +
  theme_bw(12)
ggsave(out("Predictions_PredictedCelltype.pdf"), w=16,h=5)

# Plot
ggplot(predDT, aes(x=UMAP_1, y=UMAP_2)) + 
  stat_summary_hex(bins = 50, aes(z=cell.type.prob),fun=mean) +
  facet_grid(. ~ cell.type) +
  theme_bw(12)
ggsave(out("Predictions_ClassProb.pdf"), w=15,h=5)




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
