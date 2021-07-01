source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
baseDir <- "POOLED_11_Replicates/"
out <- dirout(baseDir)

require(limma)
require(edgeR)


# SETUP -------------------------------------------------------------------



# Load data ---------------------------------------------------------------
m <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix))
str(m)
ann <- fread(PATHS$POOLED$DATA$annotation)
ann[,Date := paste0(Date, "_",Date2)]
stopifnot(all(ann$sample == colnames(m)))



# Function 2 Clean names ----------------------------------------------------
# cleanCoeff <- function(x){
#   x[x == "Mye_Cas9"] <- "Effect Mye"
#   x[x == "GenotypeCas9"] <- "Effect Und"
#   x[x == "GenotypeCas9:PopulationMye"] <- "Interaction"
#   x
# }

# summary of replicates, by library and system
ann[System == "DM" & grepl("^CFSE", Population), System := "DM.CFSE"]
ann[System == "DM" & grepl("^CD34", Population), System := "DM.CD34"]
(ann_replicates <- ann[,.N, by=c("Library", "System", "Genotype", "Population")][N>=2][order(N)])
(ann_replicates <- ann_replicates[,.N, by=c("Library", "System")][N>=2])


# Prep Annotation file ----------------------------------------------------
ai <- 2
for(ai in 1:nrow(ann_replicates)){
  message(ai)
  sysx <- ann_replicates[ai]$System
  libx <- ann_replicates[ai]$Library
  
  lAnn <- ann[Library == libx & System == sysx]
  
  baselines.vec <- c("Und", "LSK", "CD34neg")
  BASELINE <- intersect(lAnn$Population, baselines.vec)
  if(length(BASELINE) == 0) BASELINE <- NA
  print(BASELINE)
  
  lAnn$Genotype <- factor(lAnn$Genotype, levels=c("WT", unique(lAnn[Genotype != "WT"]$Genotype)))
  if(!is.na(BASELINE)){
    lAnn$Population <- factor(lAnn$Population, levels=c(BASELINE, unique(lAnn[Population != BASELINE]$Population)))
  } else {
    lAnn$Population <- factor(lAnn$Population)
  }
  BASELINE <- levels(lAnn$Population)[1]
  
  wtx <- "WT" %in% lAnn$Genotype

  analysis.type <- "guides"
  for(analysis.type in c("guides", "controls")){
    print(analysis.type)
    analysisx <- paste0(sysx,"_",libx, "_", analysis.type)
    outStats <- dirout(paste0(baseDir, "/analysis_",analysisx, "/"))
    
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
    if(wtx){
      desMat <- model.matrix(~ Genotype * Population + Date, data=data.frame(lAnn, row.names = lAnn$sample))
    } else {
      desMat <- model.matrix(~ Population + Date, data=data.frame(lAnn, row.names = lAnn$sample))
    }
    
    for(xx in c("Population", "Genotype", "Date")){
      colnames(desMat) <- gsub(xx, "", colnames(desMat))
    }
    
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
    write.tsv(data.table(voomRes.batchRemoval, keep.rownames = TRUE), outStats("Data_DateRemoved.csv"))
    
    # Model fit ---------------------------------------------------------------
    fit <- lmFit(voomRes, design=desMat)
    fit <- eBayes(fit)
    
    # Get results -------------------------------------------------------------
    res <- data.table()
    coefs <- colnames(coef(fit))
    coefx <- coefs[1]
    for(coefx in coefs){
      res <- rbind(res, data.table(topTable(fit, coef = coefx, number = nrow(lMT)),
                                   coef=coefx, analysis=analysisx, keep.rownames = TRUE))
    }

    
    # Add additional comparisons and clean names ------------------------------------------------------------
    not.baseline <- setdiff(levels(lAnn$Population), BASELINE)
    
    # setup contrasts
    coef.all <- colnames(coef(fit))
    coef.used <- c(grep("\\:", coef.all, value=TRUE), "Cas9")
    contr.use <- data.frame(row.names = c(coef.used, setdiff(coef.all, coef.used)))
    if(wtx){
      
      # Clean coefficients
      res[coef == "Cas9", coef := paste0("Cas9_",BASELINE)]
      res[grepl("\\:", coef), coef := paste0(gsub("Cas9:", "", coef), "vs", BASELINE)]
      
      
      # Define contrasts (Cas9 vs WT in each population)
      for(popx in not.baseline){
        des.vec <- rep(0, length(coef.all))
        des.vec[row.names(contr.use) == "Cas9"] <- 1
        des.vec[row.names(contr.use) == paste0("Cas9:", popx)] <- 1
        contr.use[[paste0("Cas9_", popx)]] <- des.vec
      }
      
      # Define contrasts (compare Cas9 effects between populations)
      if(length(not.baseline) > 1){
        for(i1 in 1:(length(not.baseline)-1)){
          for(i2 in (i1+1):length(not.baseline)){
            pop1 <- not.baseline[i1]
            pop2 <- not.baseline[i2]
            des.vec <- rep(0, nrow(contr.use))
            des.vec[which(row.names(contr.use) == paste0("Cas9:", pop1))] <- 1
            des.vec[which(row.names(contr.use) == paste0("Cas9:", pop2))] <- -1
            contr.use[[paste0(pop1, "vs", pop2)]] <- des.vec
          }
        }
      }
    } else {
      for(bsx in not.baseline){
        res[coef == bsx, coef := paste0(bsx, "vs",BASELINE)]
      }
      
      if(length(not.baseline) > 1){
        for(i1 in 1:(length(not.baseline)-1)){
          for(i2 in (i1+1):length(not.baseline)){
            pop1 <- not.baseline[i1]
            pop2 <- not.baseline[i2]
            des.vec <- rep(0, nrow(contr.use))
            des.vec[which(row.names(contr.use) == pop1)] <- 1
            des.vec[which(row.names(contr.use) == pop2)] <- -1
            contr.use[[paste0(pop1, "vs", pop2)]] <- des.vec
          }
        }
      }
    }
    
    if(ncol(contr.use) > 0){
      # Plot contrasts
      cleanDev(); pdf(outStats("Design_Contrasts.pdf"), w=8, h=6)
      pheatmap(contr.use, cluster_cols = F)
      dev.off()
      
      # Fit contrasts
      fit2 <- contrasts.fit(fit, contrasts = contr.use[coef.all,,drop=F])
      fit2 <- eBayes(fit2)
      for(coefx in colnames(contr.use)){
        res <- rbind(res, data.table(topTable(fit2, coef = coefx, number = nrow(lMT)),
                                     coef=coefx,analysis=analysisx, keep.rownames = TRUE))
      }
    }
  
    # Clean coefficients ------------------------------------------------------
    res[,padj := adj.P.Val]
    res[,gene := gsub("_.+", "", rn)]
    res[,binExp := factor(floor(AveExpr))]
    res[,significant := ifelse(padj < 0.05 & abs(logFC) > 1, "significant", "non-significant")]
    res[,significant_TF := significant == "significant"]
    table(res[padj < 0.05 & abs(logFC) > 1]$coef)
    
    print(unique(res$coef))
    
    write.tsv(res, outStats("AllResults.tsv"))
  }
}

ff <- list.files(out(""), recursive = TRUE, pattern="AllResults.tsv", full.names = TRUE)
res <- do.call(rbind, lapply(ff, fread))
table(res$coef, res$analysis)

# Diagnostic plots -------------------------------------------------------
# p-value histogram
ggplot(res, aes(x=P.Value, fill=binExp)) +
  geom_histogram() +
  facet_wrap(~analysis + coef, scale="free_y") +
  theme_bw(12)
ggsave(out("PValue_histogram.pdf"), w=30,h=30)
  
  
# Numbers
ggplot(res[significant_TF == TRUE][coef != "(Intercept)"], aes(x=coef)) +
  geom_bar() +
  theme_bw(12) +
  facet_wrap(~analysis, scale="free") +
  xRot() + 
  scale_y_log10() 
ggsave(out("Number.pdf"), w=12,h=12)

# Vulcano plot
ggplot(res[coef != "(Intercept)"], aes(x=logFC, y=-log10(P.Value))) +
  geom_hex() +
  facet_wrap(~analysis + coef, scale="free") +
  theme_bw(12)
ggsave(out("VulcanoPlots.pdf"), w=30,h=30)
  


# res <- fread(outStats("AllResults.tsv"))
#   
#   xx <- res[coef %in% c("Cas9", "Cas9:Mye")][,c("rn", "logFC", "padj", "coef")]
#   xx <- merge(xx[coef == "Cas9"], xx[coef=="Cas9:Mye"], suffixes = c("_g", "_i"), by="rn")
#   xx$label <- with(xx, paste(ifelse(padj_g < 0.05, "guide", ""), ifelse(padj_i < 0.05, "interaction", "")))
#   
#   ggplot(xx, aes(x=logFC_g, y=logFC_i, color=label)) + 
#     geom_point()
#   ggsave(outStats("Interaction_Analysis.pdf"), w=5, h=5)

ax <- unique(res$analysis)[2]
for(ax in unique(res[!grepl("_controls", analysis)]$analysis)){
  message(ax)
  resx <- res[analysis == ax][grepl("vs", coef) | grepl("^Cas9_", coef)]
  
  data.br <- fread(out("analysis_", ax, "/", "Data_DateRemoved.csv"))
  rns <- data.br$rn
  data.br <- as.matrix(data.br[, -"rn"])
  row.names(data.br) <- rns
  
  sig.res <- resx[!coef %in% unique(ann$Date)][!coef == "(Intercept)"]
  sig.res[, keep := sum(padj < 0.05) >= 2 & (all(sign(logFC[padj < 0.05]) > 0) | all(sign(logFC[padj < 0.05]) < 0)), by=c("gene", "coef")]
  if(nrow(sig.res[keep == TRUE]) == 0) next
  
  hx <- length(unique(sig.res[keep==TRUE]$gene)) * 1 + 1
  
  (p_coef <- ggplot(sig.res[gene %in% sig.res[keep==TRUE]$gene], 
                    aes(x=coef, y=rn, fill=logFC, size=-log10(padj), color=padj < 0.05)) + 
      geom_point(shape=21) +
      scale_color_manual(values=c("TRUE"="black", "FALSE"="white")) +
      scale_fill_gradient2(low="red", high="blue") +
      facet_grid(gene~.,scales = "free", space = "free") +
      theme_bw(12))
  ggsave(out("Results_", ax, "_Coefficients.pdf"), w=4*wiAdj,h=hx, plot=p_coef)
  
  
  annx <- ann[paste0(System, "_", Library) == gsub("_[a-z]+$", "", ax)]
  pDT <- do.call(rbind, lapply(sig.res$rn, function(guide){data.table(annx, guide=guide, log2cpm=data.br[guide, annx$sample])}))
  pDT[,gene := gsub("_.+", "", guide)]
  pDT[,zscores := scale(log2cpm), by=c("guide")]
  (p_vals <- ggplot(pDT[gene %in% sig.res[keep == TRUE]$gene], aes(x=gsub("^.+-(.+).txt$", "\\1", sample), y=guide, fill=zscores)) + 
      theme_bw(12) + 
      geom_tile() + 
      facet_grid(gene~Genotype + Population, scales ="free", space = "free") +
      scale_fill_gradient2(low="red", high="blue"))
  ggsave(out("Results_", ax, "_HM.pdf"), w=5*wiAdj,h=hx)
  
  
  p <- gridExtra::grid.arrange(
    p_vals,
    p_coef, 
    nrow=1, ncol=2, widths=c(5,4))
  ggsave(out("Results_", ax, "_combined.pdf"), h=hx, w=9*wiAdj, plot=p)
  
  
  write.tsv(sig.res, out("Results_", ax, "_Model_results.tsv"))
  
  (p_coef2 <- ggplot(resx[!coef %in% unique(ann$Date)][!coef == "(Intercept)"], 
                     aes(x=coef, y=rn, fill=logFC, size=-log10(padj), color=padj < 0.05)) + 
      geom_point(shape=21) +
      scale_color_manual(values=c("TRUE"="black", "FALSE"="white")) +
      scale_fill_gradient2(low="red", high="blue") +
      facet_wrap(~gene,scales = "free_y", ncol=10) +
      theme_bw(12) +
      xRot() +
      theme(axis.text.y = element_blank()))
  ggsave(out("Results_", ax, "_AllGenes_Coefficients.pdf"), w=12,h=10*wiAdj, plot=p_coef2)
}