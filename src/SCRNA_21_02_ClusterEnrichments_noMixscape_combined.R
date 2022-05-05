source("src/00_init.R")
out <- dirout("SCRNA_21_02_ClusterEnrichments_simple/")


require(doMC)
registerDoMC(10)


# Load annotation ---------------------------------------------------------
inDir <- dirout_load("SCRNA_20_Summary")
annList <- list()
for(fx in list.files(inDir(""), pattern="*singleR$")){
  annList[[fx]] <- fread(inDir(fx, "/Annotation.tsv"))
}
annListOriginal <- rbindlist(annList, fill = TRUE)
annList <- copy(annListOriginal)


# Cell types --------------------------------------------------------------
celltype.ann <- CLEAN.CELLTYPES

# Load cluster numbers ----------------------------------------------------
clusters <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle_Clusters.RDS"))
annList$Cluster.number <- clusters[match(annList$rn, rn)]$functional.cluster
# Remove cells without guides
annList <- annList[!is.na(mixscape_class)]
# Use only lin- in vivo?
# annList <- annList[tissue != "in.vivo" | markers == "lin-"]
annList[, gene := gsub("_.+$", "", CRISPR_Cellranger)]


# Identify and remove or label perturbed clusters (clusters with only perturbations) ---------------------
clDT <- dcast.data.table(annList[,.N, by=c("mixscape_class.global", "tissue", "Cluster.number")], tissue + Cluster.number ~ mixscape_class.global, value.var = "N")
clDT[, frac := KO/NTC]

# remove clusters from in vivo
(clDT.remove <- clDT[tissue == "in.vivo"][(NTC < 5 | is.na(NTC)) & (frac > 25 | is.na(frac))])
annList[tissue == "in.vivo" & Cluster.number %in% clDT.remove$Cluster.number, Clusters := "remove"]
write.tsv(clDT.remove, out("ClustersRemoved_in.vivo.tsv"))

# Label clusters in leukemia
(clDT.remove <- clDT[tissue == "leukemia" & frac > 30])
annList[tissue == "leukemia" & Cluster.number %in% clDT.remove$Cluster.number, Clusters := "novel"]
write.tsv(clDT.remove, out("ClustersRemoved_leukemia.tsv"))

# remove clusters to remove
annList <- annList[Clusters != "remove"]


# Define sets to analyze ------------------------------------------------------
fish.test.sets <- list()
TISSUES <- unique(annList$tissue)
CELLTYPES <- unique(annList$Clusters)

# basic analysis of everything
for(tx in TISSUES){
  fish.test.sets[[paste("basic", tx, sep="_")]] <- annList[tissue == tx]
}

for(tx in TISSUES){
  x <- annList[tissue == tx]
  x <- merge(x, celltype.ann, by.x="Clusters", by.y="Name")[,-"Clusters",with=F]
  x[,Clusters := Type]
  fish.test.sets[[paste("broad", tx, sep="_")]] <- x
}


# IN VIVO SETS
tx <- "in.vivo"
x <- annList[tissue == tx]

# without b cells
fish.test.sets[[paste("noBcells", tx, sep="_")]] <- x[!(grepl("B.cell", Clusters) | grepl("CLP", Clusters))]

# Early branching analysis
eba <- list(MEP=c("MEP (early)", "MEP"),GMP=c("GMP (early)", "GMP"))
xxx <- x[Clusters %in% do.call(c, eba)]
for(xnam in names(eba)){xxx[Clusters %in% eba[[xnam]], Clusters := xnam]}
fish.test.sets[[paste("earlyBranches", tx, sep="_")]] <- xxx

# Erythroid vs lymphoid
eba <- list(
  MEP=c(grep("Ery", CELLTYPES, value=TRUE), grep("MEP", CELLTYPES, value=TRUE)),
  GMP=c(grep("Gran", CELLTYPES, value=TRUE), grep("GMP", CELLTYPES, value=TRUE))
)
xxx <- x[Clusters %in% do.call(c, eba)]
for(xnam in names(eba)){xxx[Clusters %in% eba[[xnam]], Clusters := xnam]}
fish.test.sets[[paste("eryVsMye", tx, sep="_")]] <- xxx

# # broad groups
# xxx <- x[!(grepl("B.cell", Clusters) | grepl("CLP", Clusters) | grepl("EBMP", Clusters))]
# eba <- list(
#   MEP=c(grep("Ery", CELLTYPES, value=TRUE), grep("MEP", CELLTYPES, value=TRUE)),
#   GMP=c(grep("Gran", CELLTYPES, value=TRUE), grep("GMP", CELLTYPES, value=TRUE))
# )
# for(xnam in names(eba)){xxx[Clusters %in% eba[[xnam]], Clusters := xnam]}
# table(xxx$Clusters)
# fish.test.sets[[paste("broadBranches", tx, sep="_")]] <- xxx
# 
# # broad groups with B cells
# xxx <- x[!Clusters %in% c("B-cell", "CLP", "EBMP")]
# eba <- list(
#   MEP=c(grep("Ery", CELLTYPES, value=TRUE), grep("MEP", CELLTYPES, value=TRUE)),
#   GMP=c(grep("Gran", CELLTYPES, value=TRUE), grep("GMP", CELLTYPES, value=TRUE))
# )
# for(xnam in names(eba)){xxx[Clusters %in% eba[[xnam]], Clusters := xnam]}
# table(xxx$Clusters)
# fish.test.sets[[paste("broadBranchesWithB", tx, sep="_")]] <- xxx

# # Terminal diff Ery
# eba <- list(
#   Ery=setdiff(c(grep("Ery", CELLTYPES, value=TRUE), grep("MEP", CELLTYPES, value=TRUE)), "MEP (early)"),
#   MEP=c("MEP (early)"))
# xxx <- x[Clusters %in% do.call(c, eba)]
# for(xnam in names(eba)){xxx[Clusters %in% eba[[xnam]], Clusters := xnam]}
# fish.test.sets[[paste("diffGran", tx, sep="_")]] <- xxx
# 
# # Terminal diff Gran
# eba <- list("Gran."=c("Gran.", "Gran. P"),GMP=c("GMP (early)", "GMP", "GMP (late)"))
# xxx <- x[Clusters %in% do.call(c, eba)]
# for(xnam in names(eba)){xxx[Clusters %in% eba[[xnam]], Clusters := xnam]}
# fish.test.sets[[paste("diffGran", tx, sep="_")]] <- xxx



# W/O Mixscape ------------------------------------------------------------
fish.test.sets <- c(
  setNames(lapply(fish.test.sets, function(dt) dt[mixscape_class.global != "NP"]), paste0(names(fish.test.sets), "_withMixscape")),
  setNames(fish.test.sets, paste0(names(fish.test.sets), "_noMixscape"))
  )

# remove those with 0 rows
sapply(fish.test.sets, nrow)
fish.test.sets <- fish.test.sets[sapply(fish.test.sets, nrow) > 0]



# START ANALYSIS
(fish.test.x <- names(fish.test.sets)[1])
foreach(fish.test.x = names(fish.test.sets)) %dopar% {
  res <- data.table()
  pDT1 <- fish.test.sets[[fish.test.x]]
  write.tsv(pDT1[,"rn",with=F], out("Guides_Fisher_Mixscape_",fish.test.x,"_Cells.tsv"))
  stopifnot(sum(is.na(pDT1$CRISPR_Cellranger)) == 0)
  sx <- pDT1$timepoint[1]
  for(sx in unique(pDT1$timepoint)){
    pDT2 <- pDT1[timepoint == sx]
    gx <- pDT2$gene[1]
    for(gx in unique(pDT2[mixscape_class != "NTC"]$gene)){
      cx <- pDT2$Clusters[1]
      for(cx in unique(pDT2$Clusters)){
        ntc <- unique(pDT2[mixscape_class == "NTC"][,.N, by="CRISPR_Cellranger"][N > 20]$CRISPR_Cellranger)[1]
        for(ntc in unique(pDT2[mixscape_class == "NTC"][,.N, by="CRISPR_Cellranger"][N > 20]$CRISPR_Cellranger)){
          mx <- as.matrix(with(pDT2[gene == gx | CRISPR_Cellranger == ntc], table(Clusters == cx, gene == gx)))
          if(all(dim(mx) == c(2,2))){
            fish <- fisher.test(mx)
            res <- rbind(res, data.table(
              Clusters=cx, 
              mixscape_class=gx, 
              ntc=ntc, 
              p=fish$p.value, 
              OR=fish$estimate, 
              sample=sx, 
              total.cells=sum(mx),
              guide.cells=nrow(pDT2[gene  == gx])
            ))
          }
        }
      }
    }
  }
  
  # save file
  if(nrow(res) < 3) next
  res[,padj := p.adjust(p, method="BH")]
  res[, log2OR := pmax(-5, pmin(5, log2(OR)))]
  res[,grp := paste(mixscape_class, sample)]
  write.tsv(res[,-c("grp"), with=F], out("Guides_Fisher_Mixscape_",fish.test.x,".tsv"))
  
  # plot
  res <- hierarch.ordering(res, toOrder = "grp", orderBy = "Clusters", value.var = "log2OR", aggregate = TRUE)
  ggplot(res, aes(
    x=Clusters,
    y=ntc,
    color=log2OR,
    size=pmin(-log10(padj), 5))) +
    geom_point(shape=16) +
    scale_color_gradient2(name="log2OR", low="blue", high="red") +
    scale_size_continuous(name="padj") +
    facet_grid(mixscape_class + sample ~ ., space = "free", scales = "free") +
    theme_bw(12) +
    theme(strip.text.y = element_text(angle=0)) +
    xRot()
  ggsave(out("Guides_Fisher_Mixscape_",fish.test.x,".pdf"), w=10, h=length(unique(res$grp)) * 0.50 + 1, limitsize = FALSE)
}



# specific gene
gg <- "Kmt2d"
fish.test.x <- "noBcells_in.vivo_withMixscape"
sx <- "in.vivo_OP2_14d"
gx <- "Rbbp4_BR_14486"
cx <- "Gran. P"
pDT <- fish.test.sets[[fish.test.x]]
pDT[, perturbed := !(mixscape_class.global %in% c("NP", "NTC") | is.na(mixscape_class.global))]
pDT[, gene := gsub("_.+", "", guide)]
pDT <- pDT[(gene == gg & perturbed == TRUE) | mixscape_class == "NTC"]
pDT <- pDT[, .N, by=c("gene", "Clusters", "sample")]
pDT[, sum := sum(N), by=c("Clusters", "sample")]
pDT[, fraction := N/sum]
#pDT[rel2NTCs == Inf, rel2NTCs := 1]
pDT <- pDT[gene != "NTC"]
ggplot(pDT, aes(x=cleanCelltypes(Clusters, reverse = FALSE), y=fraction, fill=sample)) + 
  themeNF(rotate=TRUE) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label=paste0(N,"/", sum))) + 
  ggtitle(gg)


annListOriginal[tissue == "in.vivo"][grepl(gg, guide)][,.N, by=c("Clusters", "tissue", "sample")]
