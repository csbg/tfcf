source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("POOLED_10_01_IndividualAnalysis/")

require(limma)
require(ggrepel)

# Load data ---------------------------------------------------------------
m <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix))
str(m)
ann <- fread(PATHS$POOLED$DATA$annotation)
#ann[,Genotype := gsub("CAS9", "Cas9", Genotype, ignore.case = TRUE)]
ann[,Date := paste0(Date, "_",Date2)]
stopifnot(all(ann$sample == colnames(m)))


# # Comparing Genotypes within Populations -----------------------------------------------------
# res <- data.table()
# 
# datex <- "Nov2019_Nov2019"
# libx <- "A"
# for(datex in unique(ann$Date)){
#   for(libx in unique(ann[Date == datex]$Library)){
#     
#     xAnn <- ann[Date == datex & Library == libx]
#     xMT <- m[,xAnn$sample]
#     #table(apply(!is.na(m[grepl("NonTargetingControl", row.names(m)),xAnn$sample]), 1, sum))
#     xMT <- xMT[apply(!is.na(xMT), 1, sum) > ncol(xMT) * 0.8,]
#     xMT <- voom(xMT, plot=FALSE)$E
#     
#     popx <- "Mye"
#     for(popx in unique(xAnn$Population)){
#       pAnn <- xAnn[Population == popx]
#       stopifnot(nrow(pAnn) == 2)
#       stopifnot(c("Cas9", "WT") %in% pAnn$Genotype)
#       res <- rbind(res, data.table(
#         Population = popx, 
#         Score=xMT[,pAnn[Genotype == "Cas9"]$sample] - xMT[,pAnn[Genotype == "WT"]$sample], 
#         Guide=row.names(xMT),
#         Date = datex,
#         Library=libx
#       ))
#     }
#   }
# }
# 
# 
# ggplot(res, aes(x=Score, color=Population)) + 
#   geom_density() +
#   scale_color_manual(values=RColorBrewer::brewer.pal(name = "Dark2", n=length(unique(res$Population)))) +
#   facet_wrap(~Library + Date, scales = "free")
# ggsave(out("GenotypeScores_Density.pdf"), w=20,h=20)
# 
# ggplot(res, aes(x=Score, color=Population)) + 
#   stat_ecdf() +
#   scale_color_manual(values=RColorBrewer::brewer.pal(name = "Dark2", n=length(unique(res$Population)))) +
#   facet_wrap(~Library + Date, scales = "free")
# ggsave(out("GenotypeScores_ECDF.pdf"), w=20,h=20)




# Get raw scores ------------------------------------------
res <- data.table()
datex <- "10022021_Feb2021"
libx <- "R2.Br"
for(datex in unique(ann$Date)){
  for(libx in unique(ann[Date == datex]$Library)){
    
    xAnn <- ann[Date == datex & Library == libx]
    xMT <- m[,xAnn$sample,drop=F]
    #table(apply(!is.na(m[grepl("NonTargetingControl", row.names(m)),xAnn$sample]), 1, sum))
    xMT <- xMT[apply(!is.na(xMT), 1, sum) > ncol(xMT) * 0.8,, drop=F]
    xMT <- voom(xMT, plot=FALSE)$E
    
    gtx <- "WT"
    for(gtx in unique(xAnn$Genotype)){
      gAnn <- xAnn[Genotype == gtx]
      #stopifnot(nrow(gAnn) == 2)
      #stopifnot(c("Cas9", "WT") %in% pAnn$Genotype)
      
      popx <- sort(unique(gAnn$Population))
      if(length(popx) < 2) next
      
      for(p1 in 1:(length(popx)-1)){
        for(p2 in (p1 + 1):(length(popx))){
      # for(p1 in (1):(length(popx))){
      #   for(p2 in (1):(length(popx))){
          if(p1 != p2){
            res <- rbind(res, data.table(
                Genotype = gtx, 
                Score=xMT[,gAnn[Population == popx[p1]]$sample] - xMT[,gAnn[Population == popx[p2]]$sample], 
                Guide=row.names(xMT),
                Population1=popx[p1],
                Population2=popx[p2],
                Date = datex,
                Library=libx
              ))
          }
        }
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
res2 <- copy(res)
res2[,PopulationX := Population1]
res2[,Population1 := Population2]
res2[,Population2 := PopulationX]
res2[,Score := -Score]
res2$PopulationX <- NULL
res <- rbind(data.table(res, Type=1), data.table(res2, Type =2))
rm(list = c("res2"))
res[,Analysis := paste(Population1, Population2, Date, Library)]

write.tsv(res, out("Results_Scores.tsv"))
#res <- fread(out("Results_Scores.tsv"))


# Get mean and sd summary statistics ------------------------------------------
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
  facet_wrap(~Library + Date + Population1 + Population2, scales = "free", ncol = 10) + 
  theme_bw(12)
ggsave(out("PopulationScores_Density_BgNormal.pdf"), w=30,h=30)


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
write.tsv(res.stats, out("Results_Pvalues.tsv"))
#res.stats <- fread(out("Results_Pvalues.tsv"))


# Plot Summaries of z-scores and number of hits ------------------------------------------

# Z scores
ggplot(res.stats[Type == 1], aes(y=z, x=paste(Genotype, GuideType), fill=paste(Genotype, GuideType))) + 
  geom_violin(color=NA) + geom_boxplot(coef=Inf, fill=NA) +
  facet_wrap(~Library + Date + Population1 + Population2, scales = "free", ncol = 10) +
  theme_bw(12) +
  ylab("z scores") + 
  scale_fill_manual(values=c("Cas9 Targeted"="#ff7f00", "WT Targeted"="#a6cee3", "Cas9 CTRL"="#a6cee3", "WT CTRL"="#a6cee3")) +
  xRot() +
  guides(fill=FALSE)
ggsave(out("PopulationScores_Z_Boxplots_all.pdf"), w=30,h=30)

ggplot(res.stats[Genotype == "Cas9" & GuideType == "Targeted"], 
       aes(y=z, x=paste(Population2, Date))) + 
  geom_violin(color=NA, fill="#ff7f00") + geom_boxplot(coef=Inf, fill=NA) +
  facet_grid(. ~ Library + Population1, scales = "free", space = "free") +
  theme_bw(12) +
  ylab("z scores") + 
  xRot()
ggsave(out("PopulationScores_Z_Boxplots_Cas9_Targeted.pdf"), w=15,h=8)

# Number of hits
pDT <- res.stats[Genotype == "Cas9" & GuideType == "Targeted"]
pDT[, direction := ifelse(z < 0, "down", "up")]
# make sure we are counting unique guides in the next step
stopifnot(all(pDT[,.N, by=c("Population1", "Population2", "Date", "Guide", "Library")][order(N)]$N == 1)) 
pDT <- pDT[,.(sig=sum(padj < 0.05), total=.N), by=c("Population1", "Population2", "Date", "direction", "Library")]
pDT[direction == "down", sig := -sig]
pDT[,percent := sig/total * 100]

ggplot(pDT, aes(y=percent, x=paste(Population2, Date), fill=direction)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c(up="#e31a1c", down="#1f78b4")) + 
  facet_grid(. ~ Library+Population1, scales = "free_x", space = "free_x") +
  theme_bw(12) +
  ylab("Significant guides (% of tested guides)") + 
  xRot()
ggsave(out("PopulationScores_N_Significant_Cas9_Targeted.pdf"), w=15,h=8)



# Plot specific genes and guides ------------------------------------------
libx <- res$Library[1]
for(libx in unique(res.stats$Library)){
  lDT <- res.stats[Library == libx][Gene != "NonTargetingControlGuideForMouse"]
  lDT[, keep := sum(padj < 0.05) >= 2 & (all(sign(Score[padj < 0.05]) > 0) | all(sign(Score[padj < 0.05]) < 0)), by=c("Gene", "Analysis")]
  lDT <- lDT[Gene %in% lDT[keep == TRUE]$Gene]
  if(nrow(lDT) == 0) next
  ggplot(lDT, aes(x=paste(Population2, Date), y=Guide, color=Score, size=pmin(5, -log10(padj)))) + 
    geom_point() +
    geom_point(data=lDT[padj < 0.05], color="black", shape=1) +
    scale_color_gradient2(name="log2FC", low="blue", high="red") +
    scale_size_continuous(name="padj (cap = 5)") + 
    facet_grid(Gene ~ Library + Genotype + Population1,  space = "free", scales = "free") + 
    theme_bw(12) +
    xlab("") +
    xRot()
  ggsave(out("PopulationScores_Results_",libx,".pdf"), w=length(unique(lDT$Analysis)) * 0.5 + 6, h=length(unique(lDT$Guide))*0.2 + 4)
}




# Plot examples -----------------------------------------------------------
m.norm <- t(t(m)/colSums(m, na.rm = TRUE)) * 1e-6
m.norm <- log2(m.norm + 1)
table(is.na(m.norm))

genex <- "Spi1"

guides <- unique(res.stats[Gene == genex]$Guide)
gx <- guides[1]
pDT <- data.table()
for(gx in guides){
  ret <- copy(ann)
  ret$NormCounts <- m.norm[gx, ret$sample]
  pDT <- rbind(pDT, data.table(ret, Guide=gx))
}
pDT <- pDT[!is.na(NormCounts)]
pDT <- pDT[Genotype == "Cas9"]
pDT[, NormCounts2 := scale(NormCounts), by="Guide"]
ggplot(pDT, aes(x=paste(Population, Date),y=Guide, fill=NormCounts2)) + 
  geom_tile() +
  facet_grid(. ~ Library,scales="free", space = "free") + 
  scale_fill_gradient2(name="Normalized\ncounts", low="blue", high="red") +
  theme_bw(12) +
  xRot()
ggsave(out("Example_", genex, ".pdf"), w=length(unique(pDT$sample)) * 0.2 + 2,h=length(unique(pDT$Guide)) * 0.2 + 3)




# Summarize MDS -----------------------------------------------------------
res.stats <- fread(out("Results_Pvalues.tsv"))
mdsDT <- res.stats[GuideType == "Targeted" & Genotype == "Cas9"]
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








# Plot individual values -------------------------------------------------------------------------
rMT <- m
tpmMT <- t(t(rMT) / colSums(rMT, na.rm=TRUE)*1e6)
stopifnot(all(is.na(rMT) == is.na(tpmMT)))
stopifnot(all(tpmMT[,1] == rMT[,1]/sum(rMT[,1],na.rm=T)*1e6, na.rm = TRUE))
cDT <- copy(ann)
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
agMT <- t(sapply(split(row.names(agMT), gsub("_.+$", "", row.names(agMT))), function(guides){
  colMeans(agMT[guides,,drop=F], na.rm=TRUE)
}))
agMT[is.nan(agMT)] <- NA

agDT <- melt(data.table(agMT, keep.rownames = TRUE), id.vars = "rn")[!is.na(value)]
agDT[,log2FC := log2(value)]
write.tsv(agDT, out("Aggregated.tsv"))
ggplot(agDT[rn %in% hit.genes], aes(x=variable, y=rn, color=log2FC)) +
  geom_point() +
  theme_bw(12) + 
  xRot() +
  scale_color_gradient2()
ggsave(out("Aggregated_log2TPM_vsWT.pdf"), w=4,h=15)



# Hopefully Cool plot ---------------------------------------------------------------
agDT <- fread(out("Aggregated.tsv"))

compDT <- unique(res.stats[,c("Population1", "Population2", "Comparison"), with=F])
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
ggplot(agDT[rn %in% hit.genes], aes(x=Population, y=rn, color=log2FC)) +
  geom_point() +
  theme_bw(12) + 
  facet_grid(. ~ cleanComparisons(Comparison.Group), scales = "free", space = "free") +
  xRot() +
  geom_segment(x=1.2,xend=1.8, y=5,yend=5, color="black") +
  scale_color_gradient2()
ggsave(out("Aggregated_log2TPM_vsWT_arrows.pdf"), w=4,h=15)



# Plot network ------------------------------------------------------------
agx <- agDT[rn == "Spi1"]
statx <- res.stats[Gene == "Spi1"][Genotype == "Cas9"]
statx[, keep := sum(padj < 0.05) >= 2 & (all(sign(Score[padj < 0.05]) > 0) | all(sign(Score[padj < 0.05]) < 0)), by=c("Gene", "Analysis")]

unique(agx[,c("variable", "log2FC"), with=F])
statx[, .(mean(z), sum(keep),n=.N), by=c("Gene", "Comparison")]

c("cKit", "LSK9", "MEP", "LSK7", "GMP", "LSK7", "MEP", "", "", "", "", "")
