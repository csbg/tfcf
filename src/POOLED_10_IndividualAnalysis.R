source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("POOLED_10_IndividualAnalysis/")

require(limma)

# Load data ---------------------------------------------------------------
m <- as.matrix(read.csv(PATHS$DATA$matrix))
str(m)
ann <- fread(PATHS$DATA$annotation)
ann[,Date := paste0(Date, "_",Date2)]
stopifnot(all(ann$sample == colnames(m)))


COLOR.Genotypes = c(WT="#33a02c", Cas9="#6a3d9a")

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




# Comparing Genotypes within Populations -----------------------------------------------------

res <- data.table()

datex <- "Nov2019_Nov2019"
libx <- "A"
for(datex in unique(ann$Date)){
  for(libx in unique(ann[Date == datex]$Library)){
    
    xAnn <- ann[Date == datex & Library == libx]
    xMT <- m[,xAnn$sample]
    #table(apply(!is.na(m[grepl("NonTargetingControl", row.names(m)),xAnn$sample]), 1, sum))
    xMT <- xMT[apply(!is.na(xMT), 1, sum) > ncol(xMT) * 0.8,]
    xMT <- voom(xMT, plot=FALSE)$E
    
    gtx <- "WT"
    for(gtx in unique(xAnn$Genotype)){
      gAnn <- xAnn[Genotype == gtx]
      #stopifnot(nrow(gAnn) == 2)
      #stopifnot(c("Cas9", "WT") %in% pAnn$Genotype)
      
      popx <- sort(unique(gAnn$Population))
      
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


stats <- data.table()
ax <- res$Analysis[1]
for(ax in unique(res$Analysis)){
  wtDT <- res[Genotype == "WT" & Analysis == ax]
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

stat.rnorm <- data.table()
for(ax in unique(stats$Analysis)){
  x <- stats[Analysis == ax]
  stat.rnorm <- rbind(stat.rnorm, data.table(x, Score=rnorm(10000, mean=x$mean, sd=x$sd)))
}


ggplot(res, aes(x=Score)) + 
  geom_density(data=stat.rnorm, fill="#b2df8a", color=NA) +
  scale_color_manual(values=COLOR.Genotypes) +
  geom_density(aes(color=Genotype)) +
  facet_wrap(~Library + Date + Population1 + Population2, scales = "free", ncol = 6) + 
  theme_bw(12)
ggsave(out("PopulationScores_Density_BgNormal.pdf"), w=20,h=20)



# Calculate p-values
# Plot p-values for wt / cas9 and for non-targeting guides
res.stats <- data.table()
ax <- res$Analysis[1]
for(ax in unique(res$Analysis)){
  x <- stats[Analysis == ax]
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

ggplot(res.stats[Type == 1], aes(y=z, x=paste(Genotype, GuideType))) + 
  geom_violin(color=NA, fill="#a6cee3") + geom_boxplot(coef=Inf, fill=NA) +
  facet_wrap(~Library + Date + Population1 + Population2, scales = "free", ncol = 6) +
  theme_bw(12) +
  xRot()
ggsave(out("PopulationScores_Z_Boxplots.pdf"), w=20,h=20)

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

