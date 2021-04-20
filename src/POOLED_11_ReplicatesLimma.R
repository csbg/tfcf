source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
baseDir <- "POOLED_11_Replicates/"
out <- dirout(baseDir)

require(limma)
require(edgeR)


# SETUP -------------------------------------------------------------------



# Load data ---------------------------------------------------------------
m <- as.matrix(read.csv(PATHS$DATA$matrix))
str(m)
ann <- fread(PATHS$DATA$annotation)
ann[,Date := paste0(Date, "_",Date2)]
stopifnot(all(ann$sample == colnames(m)))



# Function 2 Clean names ----------------------------------------------------
cleanCoeff <- function(x){
  x[x == "Mye_Cas9"] <- "Effect Mye"
  x[x == "GenotypeCas9"] <- "Effect Und"
  x[x == "GenotypeCas9:PopulationMye"] <- "Interaction"
  x
}


# Prep Annotation file ----------------------------------------------------

BASELINE <- "Und"
lAnn <- ann[Library == "As" & System == "Mye"]
lAnn$Genotype <- factor(lAnn$Genotype, levels=c("WT", unique(lAnn[Genotype != "WT"]$Genotype)))
lAnn$Population <- factor(lAnn$Population, levels=c("Und", unique(lAnn[Population != "Und"]$Population)))


analysis.type <- "guides"
for(analysis.type in c("guides", "controls")){
  outStats <- dirout(paste0(baseDir, "/analysis_",analysis.type, "/"))
  
  # Prepare expression matrix -----------------------------------------------
  lMT <- m[,lAnn$sample]
  table(apply(!is.na(lMT), 1, sum))
  table(apply(!is.na(lMT), 2, sum))
  table(apply(!is.na(m[grepl("NonTargetingControl", row.names(m)),lAnn$sample]), 1, sum))
  lMT <- lMT[apply(!is.na(lMT), 1, sum) > ncol(lMT) * 0.8,]
  #lMT[grep("Brd9", row.names(lMT)),]
  
  # Control analysis -----------------------------------------------------------
  control.guides <- grepl("NonTargetingControl", row.names(lMT))
  lMT <- if(analysis.type == "controls") lMT else lMT[!control.guides,]
  wiAdj <- if(analysis.type == "controls") 2 else 1
  
  
  # Design matrix -----------------------------------------------------------
  desMat <- model.matrix(~ Genotype * Population + Date, data=data.frame(lAnn, row.names = lAnn$sample))
  
  cleanDev(); pdf(outStats("Design.pdf"), w=8, h=6)
  pheatmap(desMat, cluster_cols = F)
  dev.off()
  
  
  # Normalize data ----------------------------------------------------------
  lMT <- calcNormFactors(DGEList(lMT))
  cleanDev(); pdf(outStats("Voom_", "_Before.pdf"), w=6,h=6)
  voomRes <- voom(lMT, design = desMat, plot=TRUE)
  dev.off()
  
  voomRes.batchRemoval <- removeBatchEffect(
    x=voomRes$E,
    batch=lAnn$Date
  )
  
  # Model fit ---------------------------------------------------------------
  fit <- lmFit(voomRes, design=desMat)
  fit <- eBayes(fit)
  
  # Get results -------------------------------------------------------------
  res <- data.table()
  coefs <- colnames(coef(fit))
  coefx <- coefs[1]
  for(coefx in coefs){
    res <- rbind(res, data.table(topTable(fit, coef = coefx, number = nrow(lMT)),
                                 coef=coefx, keep.rownames = TRUE))
  }

  
  
  # Contrast fit ------------------------------------------------------------
  coef.all <- colnames(coef(fit))
  coef.used <- c(grep("\\:", coef.all, value=TRUE), "GenotypeCas9")
  contr.use <- data.frame(row.names = c(coef.used, setdiff(coef.all, coef.used)))
  contr.use[[paste0("Cas9_", lAnn[Population != BASELINE]$Population[1])]] <- c(c(1,1), rep(0, length(coef.all) -2))
  fit2 <- contrasts.fit(fit, contrasts = contr.use[coef.all,,drop=F])
  fit2 <- eBayes(fit2)
  coefx <- colnames(contr.use)
  res <- rbind(res, data.table(topTable(fit2, coef = coefx, number = nrow(lMT)),
                               coef=coefx, keep.rownames = TRUE))
  

  # Clean coefficients ------------------------------------------------------
  for(xx in names(attr(desMat, "contrasts"))){
    res[,coef := gsub(xx, "", coef)]
  }
  
  
  res[,padj := adj.P.Val]
  res[,gene := gsub("_.+", "", rn)]
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
  ggsave(outStats("PValue_histogram.pdf"), w=10,h=10)
  
  
  # Numbers
  ggplot(res[significant_TF == TRUE], aes(x=coef)) +
    geom_bar() +
    theme_bw(12) +
    xRot() + 
    scale_y_log10() 
  ggsave(outStats("Number.pdf"), w=5,h=5)
  
  # Vulcano plot
  ggplot(res, aes(x=logFC, y=-log10(P.Value))) +
    geom_hex() +
    facet_wrap(~coef, scale="free") +
    theme_bw(12)
  ggsave(outStats("VulcanoPlots.pdf"), w=10,h=10)
  
  write.tsv(res, outStats("AllResults.tsv"))
  res <- fread(outStats("AllResults.tsv"))
  
  xx <- res[coef %in% c("Cas9", "Cas9:Mye")][,c("rn", "logFC", "padj", "coef")]
  xx <- merge(xx[coef == "Cas9"], xx[coef=="Cas9:Mye"], suffixes = c("_g", "_i"), by="rn")
  xx$label <- with(xx, paste(ifelse(padj_g < 0.05, "guide", ""), ifelse(padj_i < 0.05, "interaction", "")))
  
  ggplot(xx, aes(x=logFC_g, y=logFC_i, color=label)) + 
    geom_point()
  ggsave(outStats("Interaction_Analysis.pdf"), w=5, h=5)
  
  
  sig.res <- res[gsub("_.+", "", rn) %in% gsub("_.+", "", xx[label != " "]$rn)][coef %in% c("Cas9", "Cas9:Mye", "Cas9_Mye")]
  sig.res[, keep := sum(padj < 0.05) >= 2 & (all(sign(logFC[padj < 0.05]) > 0) | all(sign(logFC[padj < 0.05]) < 0)), by=c("gene", "coef")]
  (p_coef <- ggplot(sig.res[gene %in% sig.res[keep==TRUE]$gene], 
                    aes(x=coef, y=rn, fill=logFC, size=-log10(padj), color=padj < 0.05)) + 
      geom_point(shape=21) +
      scale_color_manual(values=c("TRUE"="black", "FALSE"="white")) +
      scale_fill_gradient2(low="blue", high="red") +
      facet_grid(gene~"x"+"y",scales = "free", space = "free") +
      theme_bw(12))
  ggsave(outStats("Interaction_Analysis_Coefficients.pdf"), w=4*wiAdj,h=15, plot=p_coef)
  
  
  pDT <- do.call(rbind, lapply(sig.res$rn, function(guide) data.table(lAnn, guide=guide, log2cpm=voomRes.batchRemoval[guide, lAnn$sample])))
  pDT[,gene := gsub("_.+", "", guide)]
  pDT[,zscores := scale(log2cpm), by=c("guide")]
  (p_vals <- ggplot(pDT[Population %in% c("Und", "Mye")][gene %in% sig.res[keep == TRUE]$gene], aes(x=gsub("^.+-(.+).txt$", "\\1", sample), y=guide, fill=zscores)) + 
      theme_bw(12) + 
      geom_tile() + 
      facet_grid(gene~Genotype + Population, scales ="free", space = "free") +
      scale_fill_gradient2(low="blue", high="red"))
  ggsave(outStats("Interaction_Analysis_HM.pdf"), w=5*wiAdj,h=15)
  
  
  p <- gridExtra::grid.arrange(
    p_vals,
    p_coef, 
    nrow=1, ncol=2, widths=c(5,4))
  ggsave(outStats("Interaction_Analysis_combined.pdf"), h=15, w=9*wiAdj, plot=p)
  
  
  write.tsv(sig.res, outStats("Interaction_model_results.tsv"))
  
  
  
  (p_coef2 <- ggplot(res[coef %in% c("Cas9", "Cas9:Mye", "Cas9_Mye")], 
                     aes(x=cleanCoeff(coef), y=rn, fill=logFC, size=-log10(padj), color=padj < 0.05)) + 
      geom_point(shape=21) +
      scale_color_manual(values=c("TRUE"="black", "FALSE"="white")) +
      scale_fill_gradient2(low="blue", high="red") +
      facet_wrap(~gene,scales = "free_y", ncol=10) +
      theme_bw(12) +
      xRot() +
      theme(axis.text.y = element_blank())
  )
  ggsave(outStats("AllGenes_Coefficients.pdf"), w=12,h=10*wiAdj, plot=p_coef2)
}