source("src/00_init.R")
out <- dirout("POOLED_12_01_Controls/")

require(limma)
require(igraph)
require(ggrepel)
require(latex2exp)

# Load data ---------------------------------------------------------------
# m <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix))
# str(m)
# ann <- fread(PATHS$POOLED$DATA$annotation)
# ann[grepl("^LSK", System) & Population == "LSK", Population := System]
# table(ann$Population, ann$System)
# stopifnot(all(ann$sample == colnames(m)))
# ann[, id := paste(Genotype, Population, Library, sep="_")]
# ann[, group := paste(Genotype, Library, sep="_")]
m2 <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix.aggregated))
ann <- fread(PATHS$POOLED$DATA$annotation.aggregated)


# normalize ---------------------------------------------------------------
m2.cpm <- t(t(m2) / colSums(m2, na.rm = TRUE)) * 1e6
stopifnot(all(colSums(m2.cpm,na.rm = TRUE) == 1e6))
pheatmap(cor(m2.cpm, m2, use="pairwise.complete.obs"), cluster_rows = F, cluster_cols = F)
stopifnot(all(diag(cor(m2.cpm, m2, use="pairwise.complete.obs"))==1))

# Figure out which samples to compare -------------------------------------
res.ratios <- data.table()
cx <- names(COMPARISONS)[1]
for(cx in names(COMPARISONS)){
  print(cx)
  p1 <- COMPARISONS[[cx]][1]
  p2 <- COMPARISONS[[cx]][2]
  message(p1, "-", p2)
  annx <- ann[Population %in% c(p1,p2)]
  annx <- annx[group %in% annx[,length(unique(paste(Population))), by=c("group")][V1 == 2]$group]
  gx <- annx$group[1]
  for(gx in unique(annx$group)){
    print(gx)
    annx2 <- annx[group == gx]
    x <- m2.cpm[,annx2[Population == p1]$id[1], drop=F] / m2.cpm[,annx2[Population == p2]$id[1], drop=F] # Division of CPM values - later log2 tranformation
    x <- setNames(data.table(group=gx, comparison=cx, data.frame(x), keep.rownames = TRUE), c("Group", "Comparison", "Guide", "Ratio"))
    res.ratios <- rbind(res.ratios, x)
  }
}
res.ratios[, control := grepl("NonTargetingControlGuideForMouse", Guide)]
res.ratios <- res.ratios[!is.na(Ratio)]
res.ratios[,Ratio.log2 := Ratio]
res.ratios[Ratio == 0,Ratio.log2 := min(res.ratios[Ratio !=0]$Ratio)]
res.ratios[Ratio == Inf,Ratio.log2 := max(res.ratios[Ratio !=Inf]$Ratio)]
res.ratios[,Ratio.log2 := log2(Ratio.log2)] # log 2 transformation
write.tsv(res.ratios, out("Results_Scores.tsv"))

# Get mean and sd summary statistics ------------------------------------------
stats <- data.table()
res.ratios[, Analysis := paste(Group, Comparison)]
sort(unique(res.ratios$Analysis))
ax <- "WT_Br UND.MEP"
(ax <- res.ratios$Analysis[1])
for(ax in unique(res.ratios$Analysis)){
  resx <- res.ratios[Analysis == ax][control == TRUE]
  if(nrow(resx) ==0) next
  x <- resx$Ratio.log2
  plot(density(x))
  # Remove outliers
  q <- quantile(x, probs=c(0.25,0.75))
  qr <- q[2]-q[1]
  x <- x[x > median(x) - 5 * qr & x < median(x) + 5 * qr]
  
  # Store in DF
  stats <- rbind(stats, data.table(unique(resx[,-c("Ratio", "Ratio.log2", "Guide"), with=F]), mean=mean(x), sd=sd(x)))
}
write.tsv(stats, out("Results_SummaryStats.tsv"))
#stats <- fread(out("Results_SummaryStats.tsv"))


# Plot calculation of p-values ------------------------------------------
stat.rnorm <- data.table()
for(ax in unique(stats$Analysis)){
  x <- stats[Analysis == ax]
  stat.rnorm <- rbind(stat.rnorm, data.table(x, Ratio.log2=rnorm(10000, mean=x$mean, sd=x$sd)))
}
#stat.rnorm[,Ratio.capped := Ratio]

p <- ggplot(res.ratios, aes(x=Ratio.log2)) + 
  geom_density(data=stat.rnorm, fill="#b2df8a", color=NA) +
  geom_density(aes(color=control)) +
  facet_grid(Group ~ sub("\\.", "\n", Comparison), scales = "free") + 
  theme_bw(12)
ggsave(out("PopulationScores_Density_BgNormal.pdf"), w=15,h=15, plot=p)


# Calculate z-scores and p-values ------------------------------------------
res.stats <- data.table()
ax <- res.ratios$Analysis[1]
for(ax in unique(res.ratios$Analysis)){
  x <- stats[Analysis == ax]
  if(nrow(x) == 0) next
  stopifnot(nrow(x) == 1)
  ret <- res.ratios[Analysis == ax]
  ret[, z := (Ratio.log2-x$mean)/x$sd]
  ret[, p := 2*pnorm(-abs(z))]
  ret[, padj := p.adjust(p, method = "BH")]
  res.stats <- rbind(res.stats, ret)
}

res.stats[,Gene := gsub("_.+", "", Guide)]
res.stats[, hit := sum(padj < 0.05) >= 2 & (all(sign(Ratio.log2[padj < 0.05]) > 0) | all(sign(Ratio.log2[padj < 0.05]) < 0)), by=c("Gene", "Analysis")]
write.tsv(res.stats, out("Results_Pvalues.tsv"))



# Plot specific genes and guides ------------------------------------------
(libx <- res.stats$Group[1])
#libx <- "Cas9_B"
for(libx in unique(res.stats$Group)){
  lDT <- res.stats[Group == libx]
  lDT <- lDT[Gene %in% lDT[hit == TRUE]$Gene]
  if(nrow(lDT) == 0) next
  lDT[,Ratio.log2 := pmin(Ratio.log2, 5)]
  lDT[,Ratio.log2 := pmax(Ratio.log2, -5)]
  ggplot(lDT, aes(x=cleanComparisons(Comparison), y=Guide, color=Ratio.log2, size=pmin(5, -log10(padj)))) + 
    geom_point() +
    geom_point(data=lDT[padj < 0.05], color="black", shape=1) +
    scale_color_gradient2(name="log2FC", low="red", high="blue") +
    scale_size_continuous(name="padj (cap = 5)") + 
    facet_grid(Gene ~ .,  space = "free", scales = "free") + 
    theme_bw(12) +
    xlab("") +
    xRot()
  ggsave(out("PopulationScores_Libraries_",libx,".pdf"), w=length(unique(lDT$Analysis)) * 0.3 + 4, h=length(unique(lDT$Guide))*0.2 + 4, limitsize = FALSE)
}


cx <- unique(res.stats$Comparison)[1]
for(cx in unique(res.stats$Comparison)){
  lDT <- res.stats[Comparison == cx]#[!grepl("^WT", Group)]
  lDT <- lDT[Gene %in% lDT[hit == TRUE]$Gene]
  if(nrow(lDT) == 0) next
  lDT[,Ratio.log2 := pmin(Ratio.log2, 5)]
  lDT[,Ratio.log2 := pmax(Ratio.log2, -5)]
  lDT[,GuideNr := rank(Guide, ties.method = "random"), by=c("Gene", "Group")]
  ggplot(lDT, aes(x=factor(GuideNr), y=Gene, color=Ratio.log2, size=pmin(5, -log10(padj)))) + 
    geom_point() +
    geom_point(data=lDT[padj < 0.05], color="black", shape=1) +
    scale_color_gradient2(name="log2FC", low="red", high="blue") +
    scale_size_continuous(name="padj (cap = 5)") + 
    facet_grid(. ~ Group,  space = "free", scales = "free") + 
    theme_bw(12) +
    xlab("Guide") +
    ggtitle(cx)
  ggsave(out("PopulationScores_Comparisons_",cx,".pdf"), w=length(unique(lDT$Group)) * 0.75 + 4, h=length(unique(lDT$Gene))*0.2 + 4, limitsize = FALSE)
}
