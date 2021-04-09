source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("POOLED_01_CollectData/")

require(readxl)
require(edgeR)
require(limma)
HM.COLORS.FUNC <- colorRampPalette(c("#6a3d9a", "#a6cee3", "white", "#fdbf6f", "#e31a1c"))

# fm <- list.files(paste0(Sys.getenv("DATA"), "raw_pooled/"), pattern=".xlsx$", full.names = TRUE)
# read_xlsx(fm,skip = 2)

ff <- list.files(paste0(Sys.getenv("DATA"), "raw_pooled/"), pattern=".txt$", full.names = TRUE)
names(ff) <- basename(ff)
ff <- lapply(ff, fread)
fnam <- names(ff)[5]
dDT <- setNames(lapply(names(ff), function(fnam){
  fx <- ff[[fnam]]
  fx <- fx[,!grepl("_score", colnames(fx)) & !grepl(".Norm", colnames(fx)) & !grepl("_norm", colnames(fx)), with=F]
  names(fx) <- paste(names(fx), fnam)
  names(fx)[grepl("^V1", names(fx))] <- "V1"
  fx
}), names(ff))

lapply(dDT, function(dt) head(dt$V1))
dDT <-lapply(dDT, function(dt){
  dt[,V1 := gsub("(\\d{5})_(\\d{6})$", "\\1", V1)]
  dt <- dt[!V1 %in% "unmapped" & !grepl("NonTargetingControl", V1),]
  dt
})

dMT <- dDT[[1]]
for(i in 2:length(dDT)){
  dMT <- merge(dMT, dDT[[i]], by="V1", all=TRUE)
}


m <- as.matrix(dMT[,-"V1",with=F])
row.names(m) <- dMT$V1

ann <- data.table(sample=colnames(m), do.call(rbind, strsplit(gsub("Wt", "WT", gsub(" ", "_", gsub("\\-", "_", colnames(m)))), "_")))
ann[V4 == "Jul2020", V2 := V3]
ann[V4 == "Jul2020", V3 := "B"]
ann.col = data.frame(row.names=ann$sample, ann[,c("V1", "V2", "V3", "V6"), with=F])

x <- unique(data.table(apply(!is.na(m), 2, sum), gsub("^.+? ", "", colnames(m))))
x$V3 <- sapply(dDT, nrow)[x$V2]
stopifnot(nrow(x[V1 != V3]) == 0)

cleanDev(); pdf(out("Correlation.pdf"), w=12, h=11)
#cMT <- corS(m[,ann[V3 !="B" & !grepl("\\d{4}19", V4)]$sample], use="pairwise.complete.obs")
cMT <- corS(m[,ann$sample], use="pairwise.complete.obs")
diag(cMT) <- NA
pheatmap(cMT, cluster_rows = F, cluster_cols = F, annotation_col = ann.col, 
         breaks=seq(-1,1, 0.01), color=HM.COLORS.FUNC(200))
dev.off()


cleanDev(); pdf(out("Jaccard.pdf"), w=12, h=11)
pheatmap(jaccard(lapply(data.table(m), function(x) row.names(m)[!is.na(x)])), 
         cluster_rows = F, cluster_cols = F, annotation_col = ann.col, 
         breaks=seq(-1,1, 0.01), color=HM.COLORS.FUNC(200),
         main="Guide overlaps, jaccard coefficient")
dev.off()



cleanCoeff <- function(x){
  x[x == "Mye_Cas9"] <- "Effect Mye"
  x[x == "V1Cas9"] <- "Effect Und"
  x[x == "V1Cas9:V2Mye"] <- "Interaction"
  x
}

lAnn <- ann[V3 == "As" & V6 == "Mye" & V7 != "Nov2019.txt"]
lAnn$V1 <- factor(lAnn$V1, levels=c("WT", unique(lAnn[V1 != "WT"]$V1)))
lAnn$V2 <- factor(lAnn$V2, levels=c("Und", unique(lAnn[V2 != "Und"]$V2)))

lMT <- m[,lAnn$sample]
table(apply(!is.na(lMT), 1, sum))
table(apply(!is.na(lMT), 2, sum))
lMT <- lMT[apply(!is.na(lMT), 1, sum) == ncol(lMT),]

desMat <- model.matrix(~ V1 * V2 + V4, data=data.frame(lAnn, row.names = lAnn$sample))

cleanDev(); pdf(out("Design.pdf"), w=8, h=6)
pheatmap(desMat, cluster_cols = F)
dev.off()

lMT <- calcNormFactors(DGEList(lMT))

cleanDev(); pdf(out("Voom_", "_Before.pdf"), w=6,h=6)
voomRes <- voom(lMT, design = desMat, plot=TRUE)
dev.off()


fit <- lmFit(voomRes, design=desMat)
fit <- eBayes(fit)

res <- data.table()
coefs <- colnames(coef(fit))
coefx <- coefs[1]
for(coefx in coefs){
  res <- rbind(res, data.table(topTable(fit, coef = coefx, number = nrow(lMT)),
                               coef=coefx, keep.rownames = TRUE))
}



# Contrast fit ------------------------------------------------------------
coef.all <- colnames(coef(fit))
coef.used <- c("V1Cas9:V2Mye", "V1Cas9")
contr.use <- data.frame(
  "Mye_Cas9"=c(c(1,1), rep(0, length(coef.all) -2)), 
  row.names = c(coef.used, setdiff(coef.all, coef.used)))
fit2 <- contrasts.fit(fit, contrasts = contr.use[coef.all,,drop=F])
fit2 <- eBayes(fit2)
coefx <- colnames(contr.use)
res <- rbind(res, data.table(topTable(fit2, coef = coefx, number = nrow(lMT)),
                             coef=coefx, keep.rownames = TRUE))



res[,padj := adj.P.Val]
res[,binExp := factor(floor(AveExpr))]
res[,significant := ifelse(padj < 0.05 & abs(logFC) > 1, "significant", "non-significant")]
res[,significant_TF := significant == "significant"]
table(res[padj < 0.05 & abs(logFC) > 1]$coef)



# Diasgnostic plots -------------------------------------------------------
# p-value histogram
ggplot(res, aes(x=P.Value, fill=binExp)) +
  geom_histogram() +
  facet_wrap(~coef, scale="free_y") +
  theme_bw(12)
ggsave(out("PValue_histogram.pdf"), w=10,h=10)


# Numbers
ggplot(res[significant_TF == TRUE], aes(x=coef)) +
  geom_bar() +
  theme_bw(12) +
  xRot() + 
  scale_y_log10() 
ggsave(out("Number.pdf"), w=5,h=5)

# Vulcano plot
ggplot(res, aes(x=logFC, y=-log10(P.Value))) +
  geom_hex() +
  facet_wrap(~coef, scale="free") +
  theme_bw(12)
ggsave(out("VulcanoPlots.pdf"), w=10,h=10)

write.tsv(res, out("Results.tsv"))
res <- fread(out("Results.tsv"))

xx <- res[coef %in% c("V1Cas9", "V1Cas9:V2Mye")][,c("rn", "logFC", "padj", "coef")]
xx <- merge(xx[coef == "V1Cas9"], xx[coef=="V1Cas9:V2Mye"], suffixes = c("_g", "_i"), by="rn")
xx$label <- with(xx, paste(ifelse(padj_g < 0.05, "guide", ""), ifelse(padj_i < 0.05, "interaction", "")))

ggplot(xx, aes(x=logFC_g, y=logFC_i, color=label)) + 
  geom_point()
ggsave(out("Interaction_Analysis.pdf"), w=5, h=5)


sig.res <- res[gsub("_.+", "", rn) %in% gsub("_.+", "", xx[label != " "]$rn)][coef %in% c("V1Cas9", "V1Cas9:V2Mye", "Mye_Cas9")]
sig.res[,gene := gsub("_.+", "", rn)]
sig.res[, keep := sum(padj < 0.05) >= 2 & (all(sign(logFC[padj < 0.05]) > 0) | all(sign(logFC[padj < 0.05]) < 0)), by=c("gene", "coef")]
(p_coef <- ggplot(sig.res[gene %in% sig.res[keep==TRUE]$gene], 
                  aes(x=coef, y=rn, fill=logFC, size=-log10(padj), color=padj < 0.05)) + 
  geom_point(shape=21) +
  scale_color_manual(values=c("TRUE"="black", "FALSE"="white")) +
  scale_fill_gradient2(low="blue", high="red") +
  facet_grid(gene~"x"+"y",scales = "free", space = "free") +
  theme_bw(12))
ggsave(out("Interaction_Analysis_Coefficients.pdf"), w=4,h=15, plot=p_coef)


pDT <- do.call(rbind, lapply(sig.res$rn, function(guide) data.table(lAnn, guide=guide, log2cpm=voomRes$E[guide, lAnn$sample])))
pDT[,gene := gsub("_.+", "", guide)]
pDT[,zscores := scale(log2cpm), by=c("guide")]
(p_vals <- ggplot(pDT[V2 %in% c("Und", "Mye")][gene %in% sig.res[keep == TRUE]$gene], aes(x=gsub("^.+-(.+).txt$", "\\1", sample), y=guide, fill=zscores)) + 
  theme_bw(12) + 
  geom_tile() + 
  facet_grid(gene~V1 + V2, scales ="free", space = "free") +
  scale_fill_gradient2(low="blue", high="red"))
ggsave(out("Interaction_Analysis_HM.pdf"), w=5,h=15)


p <- gridExtra::grid.arrange(
  p_vals,
  p_coef, 
  nrow=1, ncol=2, widths=c(5,4))
ggsave(out("Interaction_Analysis_combined.pdf"), h=15, w=9, plot=p)




(p_coef2 <- ggplot(sig.res, 
                  aes(x=cleanCoeff(coef), y=rn, fill=logFC, size=-log10(padj), color=padj < 0.05)) + 
    geom_point(shape=21) +
    scale_color_manual(values=c("TRUE"="black", "FALSE"="white")) +
    scale_fill_gradient2(low="blue", high="red") +
    facet_wrap(~gene,scales = "free_y", ncol=10) +
    theme_bw(12) +
    xRot() +
    theme(axis.text.y = element_blank())
    )
ggsave(out("AllGenes_Coefficients.pdf"), w=12,h=8, plot=p_coef2)

write.tsv(sig.res, out("Interaction_model_results.tsv"))
