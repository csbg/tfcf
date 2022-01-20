source("src/00_init.R")
out <- dirout("POOLED_10_03_IndividualAnalysis_NormFactors_Controls/")

require(limma)
require(igraph)
require(ggrepel)
require(latex2exp)
require(edgeR)

# Load data ---------------------------------------------------------------
m2 <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix.aggregated))
ann <- fread(PATHS$POOLED$DATA$annotation.aggregated)


# AGGREGATED VALUES for genes compared to WT -------------------------------------------------------------------------
rMT <- m2[,unique(ann[Library != "A"]$id)]
rMT <- rMT[apply(!is.na(rMT), 1, sum) > 0,]
tpmMT <- t(t(rMT) / colSums(rMT, na.rm=TRUE)*1e6)
stopifnot(all(is.na(rMT) == is.na(tpmMT)))
stopifnot(all(tpmMT[,1] == rMT[,1]/sum(rMT[,1],na.rm=T)*1e6, na.rm = TRUE))

# Compare to WT WT
cDT <- copy(ann[Library != "A"])
cDT[,id2 := paste(Population, Library)]
cDT <- dcast.data.table(unique(cDT[,c("id", "id2", "Genotype")]), id2 ~ Genotype, value.var = "id")
cDT <- cDT[!is.na(WT) & !is.na(Cas9)]
normMT <- sapply(cDT$id2, function(i){
  tpmMT[,cDT[id2 == i]$Cas9] / tpmMT[,cDT[id2 == i]$WT]
})

# Aggregate across libraries
agMT <- sapply(split(colnames(normMT), gsub(" .+$", "", colnames(normMT))), function(ids){
  rowMeans(normMT[,ids], na.rm=TRUE)
})

# Aggregate for genes
agMT <- t(sapply(split(row.names(agMT), gsub("_.+$", "", row.names(agMT))), function(guides){
  colMeans(agMT[guides,,drop=F], na.rm=TRUE)
}))
stopifnot(all(apply(!is.na(agMT), 1, sum) > 0))

agDT <- setNames(melt(data.table(agMT, keep.rownames = TRUE), id.vars = "rn")[!is.na(value)], c("Gene", "Population", "FC"))
agDT[,log2FC := log2(FC)]
write.tsv(agDT, out("ComparisonToWT.tsv"))
# ggplot(agDT[rn %in% hit.genes], aes(x=variable, y=rn, color=log2FC)) +
#   geom_point() +
#   theme_bw(12) + 
#   xRot() +
#   scale_color_gradient2()
# ggsave(out("Aggregated_log2TPM_vsWT.pdf"), w=4,h=15)



# CALCULATE P-VALUES ----------------------------------------------------

# . Get raw scores (log2FCs) ------------------------------------------
res <- data.table()
datex <- "10022021_Feb2021"
libx <- "Br"
for(libx in unique(ann$Library)){
  
  xAnn <- ann[Library == libx]
  xMT <- m2[,unique(xAnn$id),drop=F]
  xMT <- xMT[apply(!is.na(xMT), 1, sum) > ncol(xMT) * 0.8,, drop=F]
  
  # Normalize data
  normfacs.controls <- calcNormFactors(DGEList(xMT[grepl("^NTC", row.names(xMT)),]))
  xMT <- calcNormFactors(DGEList(xMT))
  stopifnot(all(colnames(xMT) == colnames(normfacs.controls)))
  stopifnot(all(row.names(xMT$samples) == row.names(normfacs.controls$samples)))
  xMT$samples$norm.factors <- normfacs.controls$samples$norm.factors
  xMT <- voom(xMT, plot=FALSE)$E
  
  # For each genotype (WT and Cas9), compare log2FC for each pair of populations
  gtx <- "WT"
  for(gtx in unique(xAnn$Genotype)){
    gAnn <- xAnn[Genotype == gtx]
    #stopifnot(nrow(gAnn) == 2)
    #stopifnot(c("Cas9", "WT") %in% pAnn$Genotype)
    
    popx <- sort(unique(gAnn$Population))
    if(length(popx) < 2) next
    
    for(cx in names(COMPARISONS)){
      pop1 <- unique(gAnn[Population == COMPARISONS[[cx]][1]]$id)
      pop2 <- unique(gAnn[Population == COMPARISONS[[cx]][2]]$id)
      if(length(pop1) > 1 | length(pop2) > 1) stop("ERROR")
      if(length(pop1) == 1 & length(pop2) == 1){
        res <- rbind(res, data.table(
          Genotype = gtx, 
          Score=xMT[,pop1] - xMT[,pop2], # Score is basically a log2FC (subtraction in log2 space)
          Guide=row.names(xMT),
          Population1=COMPARISONS[[cx]][1],
          Population2=COMPARISONS[[cx]][2],
          Comparison=cx,
          Library=libx
        ))
      }
    }
  }
}

# ggplot(res, aes(x=Score, color=Genotype)) + 
#   geom_density() +
#   scale_color_manual(values=COLOR.Genotypes) +
#   #scale_color_manual(values=RColorBrewer::brewer.pal(name = "Dark2", n=length(unique(res$Population)))) +
#   facet_wrap(~Library + Date + Population1 + Population2, scales = "free", ncol=6)
# ggsave(out("PopulationScores_Density.pdf"), w=20,h=20)
# 
# ggplot(res, aes(x=Score, color=Genotype)) + 
#   stat_ecdf() +
#   scale_color_manual(values=COLOR.Genotypes) +
#   #scale_color_manual(values=RColorBrewer::brewer.pal(name = "Dark2", n=length(unique(res$Population)))) +
#   facet_wrap(~Library + Date + Population1 + Population2, scales = "free", ncol=6)
# ggsave(out("PopulationScores_ECDF.pdf"), w=20,h=20)

# Add second parallel comparisons
# res2 <- copy(res)
# res2[,PopulationX := Population1]
# res2[,Population1 := Population2]
# res2[,Population2 := PopulationX]
# res2[,Score := -Score]
# res2$PopulationX <- NULL
# res <- rbind(data.table(res, Type=1), data.table(res2, Type =2))
# rm(list = c("res2"))
# 

write.tsv(res, out("Results_Scores.tsv"))
#res <- fread(out("Results_Scores.tsv"))


# . Get mean and sd summary statistics ------------------------------------------
res[,Analysis := paste(Library, Comparison)]
stats <- data.table()
ax <- res$Analysis[1]
for(ax in unique(res$Analysis)){
  wtDT <- res[Genotype == "WT" & Analysis == ax]
  if(nrow(wtDT) ==0) next
  x <- wtDT$Score
  # x <- x[x > quantile(x, probs = c(0.01))]
  # x <- x[x < quantile(x, probs = c(0.99))]
  stats <- rbind(stats, data.table(unique(wtDT[,-c("Genotype", "Score", "Guide"), with=F]), mean=mean(x), sd=sd(x)))
  # ggplot(wtDT, aes(x=Score)) + geom_density()
  # ggplot(res[Analysis == ax], aes(x=Score, color=Genotype)) + 
  #   geom_density() + 
  #   stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = sd(x)),color="red")
}

write.tsv(stats, out("Results_SummaryStats.tsv"))
#stats <- fread(out("Results_SummaryStats.tsv"))


# . Plot calculation of p-values ------------------------------------------
stat.rnorm <- data.table()
for(ax in unique(stats$Analysis)){
  x <- stats[Analysis == ax]
  stat.rnorm <- rbind(stat.rnorm, data.table(x, Score=rnorm(10000, mean=x$mean, sd=x$sd)))
}


ggplot(res, aes(x=Score)) + 
  geom_density(data=stat.rnorm, fill="#b2df8a", color=NA) +
  scale_color_manual(values=COLOR.Genotypes) +
  geom_density(aes(color=Genotype)) +
  facet_grid(Library ~ Comparison, scales = "free") + 
  theme_bw(12)
ggsave(out("PopulationScores_Density_BgNormal.pdf"), w=20,h=15)


# . Calculate z-scores and p-values ------------------------------------------
res.stats <- data.table()
ax <- res$Analysis[1]
for(ax in unique(res$Analysis)){
  x <- stats[Analysis == ax]
  if(nrow(x) == 0) next
  stopifnot(nrow(x) == 1)
  ret <- res[Analysis == ax]
  ret[, z := (Score-x$mean)/x$sd]
  ret[, p := 2*pnorm(-abs(z))]
  ret[, padj := p.adjust(p, method = "BH")]
  res.stats <- rbind(res.stats, ret)
}

res.stats[, GuideType := ifelse(grepl("^NTC", Guide), "CTRL", "Targeted")]
res.stats[,Gene := gsub("_.+", "", Guide)]
res.stats[, hit := sum(padj < 0.05) >= 2 & (all(sign(Score[padj < 0.05]) > 0) | all(sign(Score[padj < 0.05]) < 0)), by=c("Gene", "Analysis")]
write.tsv(res.stats, out("Results_Pvalues.tsv"))
#res.stats <- fread(out("Results_Pvalues.tsv"))


# . Plot Summaries of z-scores and number of hits ------------------------------------------

# Z scores
ggplot(res.stats, aes(y=z, x=paste(Genotype, GuideType), fill=paste(Genotype, GuideType))) + 
  geom_violin(color=NA) + geom_boxplot(coef=Inf, fill=NA) +
  facet_grid(Library~Comparison) +
  theme_bw(12) +
  ylab("z scores") + 
  scale_fill_manual(values=c("Cas9 Targeted"="#ff7f00", "WT Targeted"="#a6cee3", "Cas9 CTRL"="#a6cee3", "WT CTRL"="#a6cee3")) +
  xRot() +
  guides(fill=FALSE)
ggsave(out("PopulationScores_Z_Boxplots_all.pdf"), w=20,h=15)

ggplot(res.stats[Genotype == "Cas9" & GuideType == "Targeted"], 
       aes(y=z, x=Comparison)) + 
  geom_violin(color=NA, fill="#ff7f00") + geom_boxplot(coef=Inf, fill=NA) +
  facet_grid(Library~.) +
  theme_bw(12) +
  ylab("z scores") + 
  xRot()
ggsave(out("PopulationScores_Z_Boxplots_Cas9_Targeted.pdf"), w=5,h=10)

# Number of hits
pDT <- res.stats[Genotype == "Cas9" & GuideType == "Targeted"]
pDT[, direction := ifelse(z < 0, "down", "up")]
# make sure we are counting unique guides in the next step
stopifnot(all(pDT[,.N, by=c("Comparison", "Guide", "Library")][order(N)]$N == 1)) 
pDT <- pDT[,.(sig=sum(padj < 0.05), total=.N), by=c("Comparison", "direction", "Library")]
pDT[direction == "down", sig := -sig]
pDT[,percent := sig/total * 100]

ggplot(pDT, aes(y=percent, x=Comparison, fill=direction)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c(up="#e31a1c", down="#1f78b4")) + 
  facet_grid(Library ~ .) +
  theme_bw(12) +
  ylab("Significant guides (% of tested guides)") + 
  xRot()
ggsave(out("PopulationScores_N_Significant_Cas9_Targeted.pdf"), w=6,h=15)



# . Plot specific genes and guides ------------------------------------------
(libx <- res$Library[1])
for(libx in unique(res.stats$Library)){
  lDT <- res.stats[Library == libx][!grepl("^NTC", Gene)]
  lDT <- lDT[Gene %in% lDT[hit == TRUE]$Gene]
  if(nrow(lDT) == 0) next
  ggplot(lDT, aes(x=cleanComparisons(Comparison), y=Guide, color=Score, size=pmin(5, -log10(padj)))) + 
    geom_point() +
    geom_point(data=lDT[padj < 0.05], color="black", shape=1) +
    scale_color_gradient2(name="log2FC", low="red", high="blue") +
    scale_size_continuous(name="padj (cap = 5)") + 
    facet_grid(Gene ~ Library + Genotype,  space = "free", scales = "free") + 
    theme_bw(12) +
    xlab("") +
    xRot()
  ggsave(out("PopulationScores_Results_",libx,".pdf"), w=length(unique(lDT$Analysis)) * 0.3 + 4, h=length(unique(lDT$Guide))*0.2 + 4)
}




# . Plot examples -----------------------------------------------------------
m.norm <- t(t(m)/colSums(m, na.rm = TRUE)) * 1e-6
m.norm <- log2(m.norm + 1)
table(is.na(m.norm))

genex <- "Kmt2d"
for(genex in c("Spi1", "Kmt2d")){
  
  if(!genex %in% res.stats$Gene){
    message(genex, " not found")
    next
  }
  
  guides <- unique(res.stats[Gene == genex]$Guide)
  gx <- guides[1]
  pDT <- data.table()
  for(gx in guides){
    ret <- copy(ann)
    ret$NormCounts <- m.norm[gx, ret$sample]
    pDT <- rbind(pDT, data.table(ret, Guide=gx))
  }
  pDT <- pDT[!is.na(NormCounts)]
  #pDT <- pDT[Genotype == "Cas9"]
  pDT[, NormCounts2 := scale(NormCounts), by="Guide"]
  ggplot(pDT, aes(x=sample,y=Guide, fill=NormCounts2)) + 
    geom_tile() +
    facet_grid(~ Genotype + Library,scales="free", space = "free") + 
    scale_fill_gradient2(name="Normalized\ncounts", low="red", high="blue") +
    theme_bw(12) +
    xRot()
  ggsave(out("Example_", genex, ".pdf"), w=length(unique(pDT$sample)) * 0.2 + 2,h=length(unique(pDT$Guide)) * 0.4 +4)
}



# . Blacklist of guides that didn't work ------------------------------------
# res.stats <- fread(out("Results_Pvalues.tsv"))
# x <- res.stats[GuideType == "Targeted" & Genotype == "Cas9"]
# x[Gene == "Kmt2d"]
# x2 <- x[Comparison == "UND.MYE"][Gene == "Kmt2d"][Library == "Br"]
# gxx <- unique(x2$Guide)
# res.stats[Guide %in% gxx][GuideType == "Targeted" & Genotype == "Cas9"]
# corS(t(m2[gxx,]), use="pairwise.complete.obs")
# length(unique(res.stats$Comparison))
# res.stats[GuideType == "Targeted" & Genotype == "Cas9"][,.(sum(padj < 0.05), .N), c("Guide", "Gene")][order(Gene)][1:100]


# . Summarize MDS -----------------------------------------------------------
mdsDT <- fread(out("Results_Pvalues.tsv"))[Library != "A"][GuideType == "Targeted" & Genotype == "Cas9"]
mdsDT.hits <- copy(mdsDT)

# Prepare matrix
mdsDT[,Comparison := gsub("^.+? ", "", Analysis)]
mdsDT <- mdsDT[,.(mean(z, na.rm=TRUE)), by=c("Gene", "Comparison")]
mdsMT <- toMT(dt=mdsDT, row="Gene", col="Comparison", val="V1")
mdsMT <- mdsMT[apply(!is.na(mdsMT), 1, sum)  == ncol(mdsMT),]

# Define hit genes
mdsDT.hits[, keep := sum(padj < 0.05) >= 2 & (all(sign(Score[padj < 0.05]) > 0) | all(sign(Score[padj < 0.05]) < 0)), by=c("Gene", "Analysis")]
hit.genes <- unique(mdsDT.hits[keep==TRUE & Gene %in% row.names(mdsMT)]$Gene)

for(typex in c("all", "hits")){
  mt <- mdsMT
  if(typex == "hits") mt <- mt[hit.genes,]
  
  # Correlation of guides
  cMT <- corS(t(mt))
  diag(cMT) <- NA
  dd <- as.dist(2-cMT)
  xx <- nrow(cMT)/7
  cleanDev(); pdf(out("Correlation_",typex,"_Guides_HM.pdf"), w=xx + 2, h=xx+1)
  pheatmap(cMT, clustering_distance_rows = dd, clustering_distance_cols = dd, breaks=seq(-1,1, 0.01), color=COLORS.HM.FUNC(200))
  dev.off()
  # MDS
  mds <- data.table(cmdscale(dd, k=2), keep.rownames = TRUE)
  write.tsv(mds, out("Correlation_",typex,"_MDS.tsv"))
  
  ggplot(mds, aes(x=V1, y=V2, label=rn)) + 
    geom_point(color="lightblue") + 
    geom_point(data=mds[rn %in% hit.genes], color="red") + 
    geom_text_repel() +
    theme_bw(12)
  ggsave(out("Correlation_",typex,"_MDS.pdf"), w=8,h=8)
  
  # Correlation of comparisons
  cMT <- corS(mdsMT)
  diag(cMT) <- NA
  dd <- as.dist(2-cMT)
  cleanDev(); pdf(out("Correlation_",typex,"_Comparisons_HM.pdf"), w=6, h=6)
  pheatmap(cMT, clustering_distance_rows = dd, clustering_distance_cols = dd, breaks=seq(-1,1, 0.01), color=COLORS.HM.FUNC(200))
  dev.off()
}