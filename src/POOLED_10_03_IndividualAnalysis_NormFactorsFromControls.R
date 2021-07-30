source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("POOLED_10_03_IndividualAnalysis_NormFactors_Controls/")

require(limma)
require(igraph)
require(ggrepel)
require(latex2exp)

# Load data ---------------------------------------------------------------
m <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix))
str(m)
ann <- fread(PATHS$POOLED$DATA$annotation)
ann[grepl("^LSK", System) & Population == "LSK", Population := System]
#ann[,Genotype := gsub("CAS9", "Cas9", Genotype, ignore.case = TRUE)]
table(ann$Population, ann$System)
ann[,Date := paste0(Date, "_",Date2)]
stopifnot(all(ann$sample == colnames(m)))
ann[, id := paste(Genotype, Population, Library, sep="_")]

#m2 <- 
m2 <- sapply(with(ann, split(sample, id)), function(sx){
  apply(m[,sx, drop=F], 1, function(row){
    if(all(is.na(row))){
      NA
    } else {
      sum(row, na.rm=TRUE)
    }
  })
})


# Get raw scores ------------------------------------------
res <- data.table()
datex <- "10022021_Feb2021"
libx <- "Br"
for(libx in unique(ann$Library)){
  
  xAnn <- ann[Library == libx]
  xMT <- m2[,unique(xAnn$id),drop=F]
  #table(apply(!is.na(m[grepl("NonTargetingControl", row.names(m)),xAnn$sample]), 1, sum))
  xMT <- xMT[apply(!is.na(xMT), 1, sum) > ncol(xMT) * 0.8,, drop=F]
  
  # Normalize data
  normfacs.controls <- calcNormFactors(DGEList(xMT[grepl("NonTargetingControl", row.names(xMT)),]))
  xMT <- calcNormFactors(DGEList(xMT))
  stopifnot(all(colnames(xMT) == colnames(normfacs.controls)))
  stopifnot(all(row.names(xMT$samples) == row.names(normfacs.controls$samples)))
  xMT$samples$norm.factors <- normfacs.controls$samples$norm.factors
  xMT <- voom(xMT, plot=FALSE)$E
  
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
          Score=xMT[,pop1] - xMT[,pop2], 
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


# Get mean and sd summary statistics ------------------------------------------
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


# Plot calculation of p-values ------------------------------------------
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
ggsave(out("PopulationScores_Density_BgNormal.pdf"), w=15,h=15)


# Calculate z-scores and p-values ------------------------------------------
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

res.stats[, GuideType := ifelse(grepl("NonTargeting", Guide), "CTRL", "Targeted")]
res.stats[,Gene := gsub("_.+", "", Guide)]
res.stats[, hit := sum(padj < 0.05) >= 2 & (all(sign(Score[padj < 0.05]) > 0) | all(sign(Score[padj < 0.05]) < 0)), by=c("Gene", "Analysis")]
write.tsv(res.stats, out("Results_Pvalues.tsv"))
#res.stats <- fread(out("Results_Pvalues.tsv"))


# Plot Summaries of z-scores and number of hits ------------------------------------------

# Z scores
ggplot(res.stats, aes(y=z, x=paste(Genotype, GuideType), fill=paste(Genotype, GuideType))) + 
  geom_violin(color=NA) + geom_boxplot(coef=Inf, fill=NA) +
  facet_grid(Library~Comparison) +
  theme_bw(12) +
  ylab("z scores") + 
  scale_fill_manual(values=c("Cas9 Targeted"="#ff7f00", "WT Targeted"="#a6cee3", "Cas9 CTRL"="#a6cee3", "WT CTRL"="#a6cee3")) +
  xRot() +
  guides(fill=FALSE)
ggsave(out("PopulationScores_Z_Boxplots_all.pdf"), w=15,h=15)

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



# Plot specific genes and guides ------------------------------------------
(libx <- res$Library[1])
for(libx in unique(res.stats$Library)){
  lDT <- res.stats[Library == libx][Gene != "NonTargetingControlGuideForMouse"]
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




# Plot examples -----------------------------------------------------------
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



# Blacklist of guides that didn't work ------------------------------------
# res.stats <- fread(out("Results_Pvalues.tsv"))
# x <- res.stats[GuideType == "Targeted" & Genotype == "Cas9"]
# x[Gene == "Kmt2d"]
# x2 <- x[Comparison == "UND.MYE"][Gene == "Kmt2d"][Library == "Br"]
# gxx <- unique(x2$Guide)
# res.stats[Guide %in% gxx][GuideType == "Targeted" & Genotype == "Cas9"]
# corS(t(m2[gxx,]), use="pairwise.complete.obs")
# length(unique(res.stats$Comparison))
# res.stats[GuideType == "Targeted" & Genotype == "Cas9"][,.(sum(padj < 0.05), .N), c("Guide", "Gene")][order(Gene)][1:100]


# Summarize MDS -----------------------------------------------------------
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




# Aggregate values -------------------------------------------------------------------------
stopifnot(all(colnames(m) == ann$sample))
rMT <- m[,ann[Library != "A"]$sample]
tpmMT <- t(t(rMT) / colSums(rMT, na.rm=TRUE)*1e6)
stopifnot(all(is.na(rMT) == is.na(tpmMT)))
stopifnot(all(tpmMT[,1] == rMT[,1]/sum(rMT[,1],na.rm=T)*1e6, na.rm = TRUE))
cDT <- copy(ann[Library != "A"])
cDT[,id := paste(Population, Library, Date, Library2, System, Date2)]
cDT <- dcast.data.table(cDT, id ~ Genotype, value.var = "sample")
cDT <- cDT[!is.na(WT)]
normMT <- sapply(cDT$id, function(i){
  tpmMT[,cDT[id == i]$Cas9] / tpmMT[,cDT[id == i]$WT]
})
normMT[normMT == Inf] <- 8

cDT2 <- data.table(id=cDT$id, do.call(rbind, strsplit(cDT$id, " ")))
cDT2[V1 == "LSK", V1 := V5]
agMT <- sapply(with(cDT2, split(id, V1)), function(ids){
  rowMeans(normMT[,ids], na.rm=TRUE)
})
agMT[is.nan(agMT)] <- NA
#agMT[grepl("Kmt2d", row.names(agMT)),]
agMT <- t(sapply(split(row.names(agMT), gsub("_.+$", "", row.names(agMT))), function(guides){
  colMeans(agMT[guides,,drop=F], na.rm=TRUE)
}))
agMT[is.nan(agMT)] <- NA

agDT <- melt(data.table(agMT, keep.rownames = TRUE), id.vars = "rn")[!is.na(value)]
agDT[,log2FC := log2(value)]
write.tsv(agDT, out("Aggregated.tsv"))
# ggplot(agDT[rn %in% hit.genes], aes(x=variable, y=rn, color=log2FC)) +
#   geom_point() +
#   theme_bw(12) + 
#   xRot() +
#   scale_color_gradient2()
# ggsave(out("Aggregated_log2TPM_vsWT.pdf"), w=4,h=15)



# Hopefully Cool plot ---------------------------------------------------------------
agDT <- fread(out("Aggregated.tsv"))
compDT <- unique(res.stats[Library != "A"][,c("Population1", "Population2", "Comparison"), with=F])
compDT[Comparison %in% c("GMP.LSK", "MEP.LSK","UND.MEP", "MYE.GMP"), Comparison.Group := "Main branch"]
compDT[is.na(Comparison.Group), Comparison.Group := Comparison]
compDT2 <- unique(melt(compDT[,-"Comparison",with=F], id.vars = "Comparison.Group")[,-"variable"])
mainBranchOrdering <- c("Und", "MEP", "LSKd7", "GMP", "Mye")
agDT$Population <- factor(agDT$variable, levels=c(mainBranchOrdering, setdiff(agDT$variable, mainBranchOrdering)))
agDT <- merge(agDT, compDT2, by.y="value", by.x="variable", allow.cartesian=TRUE)

unique(agDT$Comparison.Group)
cleanComparisons <- function(x){
  orderX <- c("Main branch", "CKIT.LSK", "GMP.MEP", "UND.MYE", "GMPcd11.DN")
  x <- factor(x, levels=orderX)
}

hit.genes2 <- hit.genes[1:10]

pDT.stats <- res.stats[Gene %in% hit.genes][Library != "A"]
pDT.stats <- merge(pDT.stats, unique(compDT[,c("Comparison", "Comparison.Group")]), by="Comparison")
# Filter 1: two significant guides and all significant guides going into the same direction
pDT.stats[, keep := sum(padj < 0.05) >= 2 & (all(sign(Score[padj < 0.05]) > 0) | all(sign(Score[padj < 0.05]) < 0)), by=c("Gene", "Analysis")]
# Filter 2:  50 % of guides being significant
pDT.stats <- pDT.stats[keep == TRUE][, .(mean(z), length(unique(Guide[padj < 0.05])),n=length(unique(Guide))), by=c("Gene", "Comparison",  "Comparison.Group", "Population1", "Population2")]
pDT.stats[,percSig := V2/n*100]
pDT.stats <- pDT.stats[percSig > 50]

mds <- fread(out("Correlation_hits_MDS.tsv"))
#agDT <- hierarch.ordering(agDT, toOrder = "rn", orderBy = "variable", value.var = "log2FC", aggregate = TRUE)
agDT$rn <- factor(agDT$rn, levels=mds$rn[hclust(dist(as.matrix(data.frame(mds[,c("V1", "V2"), with=F]))))$order])
ggplot(agDT[rn %in% pDT.stats$Gene], aes(x=Population, y=rn)) +
  theme_bw(12) + 
  geom_point(aes(fill=log2FC), shape=21, color="white", size=5) +
  facet_grid(. ~ cleanComparisons(Comparison.Group), scales = "free", space = "free") +
  geom_segment(data=pDT.stats, aes(xend=Population1, x=Population2, y=Gene, yend=Gene, color=V1), arrow=arrow(type="closed", length = unit(0.3, "cm"))) + 
  scale_fill_gradient2(name=TeX(r'($\\overset{\Delta_{Cas9-WT}}{(dots)}$)')) +
  #geom_point(aes(fill=log2FC), shape=21, color="white", size=2) +
  scale_color_gradient2(name=TeX(r'($\\overset{\Delta_{Populations}}{(arrows)}$)')) +
  xRot()
ggsave(out("Aggregated_Edges.pdf"), w=8,h=13)




# Selected comparisons (David) --------------------------------------------
pDT.stats <- res.stats[Gene %in% hit.genes][Library != "A"]
pDT.stats <- merge(pDT.stats, unique(compDT[,c("Comparison", "Comparison.Group")]), by="Comparison")
# Filter 1: two significant guides and all significant guides going into the same direction
pDT.stats[, keep := sum(padj < 0.05) >= 2 & (all(sign(Score[padj < 0.05]) > 0) | all(sign(Score[padj < 0.05]) < 0)), by=c("Gene", "Analysis")]
# Filter 2:  50 % of guides being significant
pDT.stats <- pDT.stats[Gene %in% pDT.stats[keep == TRUE]$Gene][, .(mean(z), length(unique(Guide[padj < 0.05])),n=length(unique(Guide))), by=c("Gene", "Comparison",  "Comparison.Group", "Population1", "Population2")]
pDT.stats[,percSig := V2/n*100]
pDT.stats <- pDT.stats[Gene %in% pDT.stats[percSig > 50]$Gene]
pDT.stats <- pDT.stats[Comparison %in% COMPARISONS.USE]

mx <- toMT(pDT.stats, row = "Gene", col = "Comparison", val = "V1")
mx[is.na(mx)] <- 0
pDT.stats$Cluster <- cutree(hclust(dist(mx)), k = 7)[pDT.stats$Gene]

ggplot(pDT.stats, aes(x=cleanComparisons(Comparison), y=Gene)) +
  theme_bw(12) + 
  geom_point(aes(color=V1), size=4) +
  scale_color_gradient2(name=TeX(r'($\\overset{\Delta_{Cas9-WT}}$)')) +
  facet_grid(Cluster ~ ., scales = "free", space = "free") +
  xRot() + 
  xlab("")
ggsave(out("SimpleHM.pdf"), w=3.5,h=13)



# Plot network ------------------------------------------------------------
genex <- "Kmt2d"
for(genex in hit.genes){
  COLORS.graph <- c("#fb9a99", "lightgrey", "#a6cee3")
  
  agx <- agDT[rn == genex]
  agx <- unique(agx[,c("variable", "log2FC"), with=F])
  
  statx <- res.stats[Gene == genex][Genotype == "Cas9"]
  statx[, keep := sum(padj < 0.05) >= 2 & (all(sign(Score[padj < 0.05]) > 0) | all(sign(Score[padj < 0.05]) < 0)), by=c("Gene", "Analysis")]
  statx <- statx[, .(
    z=mean(z), 
    up=length(unique(Guide[padj < 0.05 & z > 0])),
    dn=length(unique(Guide[padj < 0.05 & z < 0])),
    n=length(unique(Guide))), 
    by=c("Gene", "Comparison")]
  
  el <- do.call(rbind, COMPARISONS)
  statx <- statx[match(row.names(el), Comparison)][!is.na(Gene)]
  el <- el[statx$Comparison,]
  
  g <- graph.edgelist(el[,2:1])
  V(g)$log2FC <- agx[match(V(g)$name, variable),]$log2FC
  E(g)$z <- statx$z
  E(g)$up <- statx$up
  E(g)$dn <- statx$dn
  E(g)$n <- statx$n
  E(g)$sig.label <- with(statx, paste(up, dn, n, sep="/"))
  E(g)$cnt <- ifelse(E(g)$z > 0, E(g)$up, E(g)$dn)
  E(g)$perc <- round(E(g)$cnt/E(g)$n) * 100
  
  g <- delete.vertices(g, is.na(V(g)$log2FC))
  V(g)$color <- mapNumericToColors(V(g)$log2FC, cols = COLORS.graph)
  V(g)$frame.color <- NA
  V(g)$label.color <- "black"
  E(g)$color <- mapNumericToColors(E(g)$z, cols = COLORS.graph)
  E(g)$width <- 2+E(g)$perc/100*3
  E(g)$arrow.width <- 2
  E(g)$arrow.size <- 1
  E(g)$label <- paste0("", round(E(g)$z, 1),"\n", "(",E(g)$sig.label, ")")
  E(g)$label.color <- "black"
  
  
  layout <- list(
    "LSKd7" = c(0,3),
    "GMP" = c(-1,1.5),
    "MEP" = c(1,1.5),
    "Mye" = c(-1,0),
    "Und" = c(1,0),
    "cKit" = c(3,3),
    "LSKd9" = c(3,2),
    "GMP.CD11bGr1" = c(3,1),
    "GMP.DN" = c(3,0))
  
  cleanDev(); pdf(out("Graph_", genex, ".pdf"), w=5,h=5)
  plot.igraph(g, layout=do.call(rbind, layout)[V(g)$name,], main=genex)
  dev.off()
}