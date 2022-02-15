

print("----------------")
print(tissue.name)
print(analysis.name)
print("----------------")

# Check data --------------------------------------------------------------
stopifnot(all(c("monocle.obj") %in% ls()))
stopifnot("clusters.final" %in% colnames(colData(monocle.obj)))

fData(monocle.obj)$gene_short_name <- row.names(fData(monocle.obj))

# Load packages and functions ---------------------------------------------
require(umap)
require(igraph)
require(nebula)
require(fgsea)
library(SingleR)
require(doMC)
source("src/FUNC_Monocle_PLUS.R")



# Load additional data ----------------------------------------------------

# Markers
marker.genes <- fread("metadata/markers.csv")

# Marker signautres
ff <- list.files(dirout_load("SCRNA_06_01_Markers")(tissue.name), pattern="Signatures_", full.names = TRUE)
names(ff) <- gsub("^Signatures_(.+?).csv$", "\\1", basename(ff))
marker.signatures <- lapply(ff, function(fx) as.matrix(read.csv(fx)))

# SingleR
sx <- monocle.obj$sample[1]
singleR.res <- list()
for(sx in unique(monocle.obj$sample)){
  ff <- list.files(dirout_load("SCRNA_05_01_SingleR")(sx), pattern = "cell_types_.*.csv", full.names = TRUE)
  singleR.res[[sx]] <- setNames(lapply(ff, fread), gsub("cell_types_(.+).csv", "\\1", basename(ff)))
}
singleR.res <- rbindlist(lapply(singleR.res, function(xx) rbindlist(xx, fill=TRUE,idcol = "db")), fill=TRUE, idcol="sample")
singleR.res[, cellname := paste0(cell, "_", sample)]
izzo.proj.cts <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjIzzo_celltypes.RDS"))
singleR.merge <- readRDS(dirout_load("SCRNA_06_02_MergeMarkers")("CellTypes_", tissue.name,".RDS"))
singleR.res <- rbind(
  singleR.res[,c("cellname", "labels", "db"),with=F],
  setNames(data.table(izzo.proj.cts[,c("rn", "functional.cluster"), with=F], db="IzzoProjected"), c("cellname", "labels", "db")),
  singleR.merge[,c("cellname", "labels", "db"),with=F]
)

# CITESEQ
(load(PATHS$SCRNA$Citeseq)) # citeseq.MT

# ANNOTATION ------------------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN, sep="\t")

# Collect ANNOTATION --------------------------------------------------------------
ann <- data.table(data.frame(colData(monocle.obj)@listData), keep.rownames = TRUE)
umap <- setNames(
  data.table(reducedDims(monocle.obj)$UMAP, keep.rownames = TRUE), 
  c("rn", "UMAP1", "UMAP2"))
ann <- merge(ann, umap, by="rn", all=TRUE)
ann[,Clusters := clusters.final]

# fix guide annotations ---------------------------------------------------
for(sx in unique(ann$sample)){
  sx <- "DM_Test1_NM_6d_1"
  if(sum(!is.na(ann[sample == sx]$CRISPR_Cellranger)) == 0 & sum(!is.na(ann[sample == sx]$CRISPR.Guide.Capture)) > 0){
    print(sx)
    # Most samples use the guides from CRIPSR_Cellranger, but a couple of samples use other pipeliens stored in CRISPR.Guide.Capture
    ann[sample == sx, CRISPR_Cellranger := CRISPR.Guide.Capture] 
  }
}
for(colx in c(
  grep("CRISPR", colnames(ann), value=T),
  grep("mixscape", colnames(ann), value=T),
  grep("guide", colnames(ann), value=T))){
  message("Setting all NAs in ", colx, " from sample WT-LSK_OP0_NM_7d_1")
  ann[sample == "WT-LSK_OP0_NM_7d_1", (colx) := NA]
}


# SETUP ENDS HERE ---------------------------------------------------------





# Output annotation -------------------------------------------------------
write.tsv(ann, out("Annotation.tsv"))


# ADDITIONAL QC  --------------------------------------------------------
qcm <- "nCount_RNA"
for(qcm in c("percent.mt", "nFeature_RNA", "nCount_RNA")){
  print(qcm)
  ann$measure <- ann[[qcm]]
  p <- ggplot(ann, aes(y=measure + 0.1, x=factor(Clusters))) +
    geom_violin(color=NA, fill="lightblue") +
    geom_boxplot(fill=NA, coef=Inf) +
    scale_y_log10() +
    theme_bw(12) +
    ylab(qcm) +
    ggtitle(qcm) +
    xRot()
  ggsave(out("QC_", qcm, "_Clusters.pdf"), w=5,h=4, plot=p)

  ggplot(ann, aes(x=UMAP1, y=UMAP2)) +
    theme_bw(12) +
    stat_summary_hex(bins=100, aes(z=measure),fun=mean) +
    scale_fill_hexbin() +
    #scale_fill_gradient(low="lightgrey", high="blue") +
    ggtitle(qcm)
  ggsave(out("QC_", qcm, "_UMAP.pdf"), w=5,h=4)
  }
ann$measure <- NULL


# SAMPLES -----------------------------------------------------------------
ggplot(ann, aes(x=sample)) +
  theme_bw() +
  xRot() +
  geom_bar()
ggsave(out("Samples_Numbers.pdf"), w=0.5 * length(unique(ann$sample)) + 1, h=6)

ggplot(ann, aes(x=UMAP1, y=UMAP2)) +
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~sample, ncol=5)
ggsave(out("Samples_UMAP.pdf"), w=5*2+2,h=ceiling(length(unique(ann$sample))/5) * 2 + 1)

ggplot(ann[mixscape_class.global == "NTC"], aes(x=UMAP1, y=UMAP2)) +
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~sample, ncol=5)
ggsave(out("Samples_UMAP_NTCs.pdf"), w=5*2+2,h=ceiling(length(unique(ann$sample))/5) * 2 + 1)

pDT <- ann[,.N, by=c("sample", "Clusters")]
pDT[,sumS := sum(N), by="sample"]
pDT[,sumC := sum(N), by="Clusters"]
pDT[,percentS := N/sumS*100]
pDT[,percentC := N/sumC*100]
ggplot(pDT, aes(y=sample,x=factor(as.numeric(Clusters)), size=percentS, color=percentC)) +
  scale_size_continuous(name="% of sample") +
  scale_color_gradient(name="% of cluster", low="black", high="red") +
  theme_bw() +
  geom_point()
ggsave(out("Samples_Clusters.pdf"), w=8,h=length(unique(pDT$sample)) * 0.3+1)


# broad samples -----------------------------------------------------------------
ggplot(ann, aes(x=UMAP1, y=UMAP2)) +
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~sample_broad, ncol=3)
ggsave(out("SamplesBroad_UMAP.pdf"), w=3*2+2,h=3 + 1)

pDT <- ann[,.N, by=c("sample_broad", "Clusters")]
pDT[,sumS := sum(N), by="sample_broad"]
pDT[,sumC := sum(N), by="Clusters"]
pDT[,percentS := N/sumS*100]
pDT[,percentC := N/sumC*100]
ggplot(pDT, aes(y=sample_broad,x=factor(as.numeric(Clusters)), size=percentS, color=percentC)) +
  scale_size_continuous(name="% of sample_broad") +
  scale_color_gradient(name="% of cluster", low="black", high="red") +
  theme_bw() +
  geom_point()
ggsave(out("SamplesBroad_Clusters.pdf"), w=8,h=length(unique(pDT$sample_broad)) * 1+1)


# Cell cycle --------------------------------------------------------------
ggplot(ann, aes(x=UMAP1, y=UMAP2)) +
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_hexbin() +
  facet_grid(. ~ Phase)
ggsave(out("CellCycle_UMAP.pdf"), w=16,h=5)



# CLUSTERS ----------------------------------------------------
# UMAP
ggplot(ann, aes(x=UMAP1, y=UMAP2)) +
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  geom_label(data=ann[,.(UMAP1=median(UMAP1), UMAP2=median(UMAP2)), by="Clusters"], aes(label=Clusters), fill="#ffffffaa")
ggsave(out("Clusters_UMAP.pdf"), w=6,h=5)

p <- ggplot(ann, aes(x=UMAP1, y=UMAP2, color=Clusters)) +
  theme_bw(12) +
  geom_point(size=0.5) +
  geom_text(data=ann[,.(UMAP1=median(UMAP1), UMAP2=median(UMAP2)), by="Clusters"], aes(label=Clusters), color="black")
ggsave(out("Clusters_UMAP.points.jpg"), w=7,h=5, plot=p)


# Cell type MARKERS  ------------------------------------------------------

# Markers
plot_genes_by_group(monocle.obj, markers = marker.genes$Name, group_cells_by = "clusters.final") + scale_size_continuous(range=c(0,5))
ggsave(out("Markers_Clusters.pdf"), w=6,h=8)

p <- plot_cells_umap_hex_NF(monocle.obj,genes = sort(marker.genes$Name))
ggsave(out("Markers_UMAP_hex.pdf"), w=30,h=30, plot=p)

p <- plot_cells_umap_hex_NF(monocle.obj, scale=TRUE, genes = sort(marker.genes$Name))
ggsave(out("Markers_UMAP_hex_scale.pdf"), w=30,h=30, plot=p)


# CELLTYPES SingleR -------------------------------------------------------
singleR.sum <- data.table()
srx <- singleR.res$db[1]
for(srx in unique(singleR.res$db)){
  print(srx)
  singleR.resX <- singleR.res[db == srx]

  # Predicted cells on UMAP
  pDT.pc <- merge(
    singleR.resX[,c("cellname", "labels"), with=F],
    ann[,c("UMAP1", "UMAP2", "rn", "Clusters"),with=F],
    by.x="cellname", by.y="rn")
  pDT.pc <- pDT.pc[labels %in% pDT.pc[,.N, by="labels"][N > 10]$labels]
  p <- ggplot(pDT.pc, aes(x=UMAP1, y=UMAP2)) +
    theme_bw(12) +
    geom_hex(bins=100) +
    scale_fill_hexbin() +
    facet_wrap(~labels, ncol=5) +
    theme_bw(12) +
    ggtitle(srx)
  ggsave(
    out("SingleR_", srx, "_UMAP_Predicted",".pdf"),
    w=5*3+2,
    h=ceiling(length(unique(pDT.pc$labels))/5)*3+1,
    limitsize = FALSE,
    plot=p)

  # Clusters - Predictions
  pDT.ann <- merge(singleR.resX[,c("labels", "cellname")],ann, by.x="cellname", by.y="rn")
  pDT.ann <- pDT.ann[,.N, by=c("Clusters", "labels", "tissue")]
  pDT.ann[,sum := sum(N), by=c("Clusters", "tissue")]
  pDT.ann[,percent := N/sum*100]
  if(!any(is.na(as.numeric(pDT.ann$Clusters)))) pDT.ann[, Clusters := as.numeric(Clusters)]
  ggplot(pDT.ann, aes(x=factor(Clusters), y=labels, fill=percent)) +
    theme_bw(12) +
    geom_tile() +
    facet_grid(. ~ tissue) +
    scale_fill_gradient(limits=c(0,100), low="white", high="red") +
    ggtitle(srx)
  ggsave(out("SingleR_", srx, "_Clusters_", "PercPredicted", ".pdf"),
         h=length(unique(pDT.ann$labels)) * 0.3+1,
         w=length(unique(pDT.ann$Clusters)) * 0.3 * 3+2,
         limitsize = FALSE)

  pDT.ann$dataset <- srx
  singleR.sum <- rbind(singleR.sum, pDT.ann)
}
singleR.sum[,id := paste(dataset, labels, tissue)]
pDT <- merge(singleR.sum[,max(percent), by=c("id")][V1 > 1][,c("id")], singleR.sum, by=c("id"))
pDT <- hierarch.ordering(pDT, toOrder = "Clusters", orderBy = "labels", value.var = "percent", aggregate = TRUE)
pDT <- hierarch.ordering(pDT, toOrder = "labels", orderBy = "Clusters", value.var = "percent", aggregate = TRUE)
ggplot(pDT, aes(y=Clusters, x=labels, fill=percent)) +
  theme_bw(12) +
  geom_tile() +
  facet_grid(tissue ~ gsub("_", "\n", dataset), scales = "free", space = "free") +
  scale_fill_gradient(limits=c(0,100), low="white", high="red") +
  xRot()
ggsave(out("SingleR_0_Clusters_", "PercPredicted", ".pdf"),
       w=nrow(pDT[,.N, by=c("dataset", "labels")]) * 0.2+2,
       h=length(unique(pDT$Clusters)) * 0.2 * 3+1,
       limitsize = FALSE)



# CellTypes from Marker signatures --------------------------------------------------
mnam <- names(marker.signatures)[1]
for(mnam in names(marker.signatures)){
  mx <- marker.signatures[[mnam]]

  stopifnot(all(ann$rn %in% row.names(mx)))

  pDT <- merge(ann[,c("rn", "UMAP1", "UMAP2")], melt(data.table(mx, keep.rownames = TRUE), id.vars = "rn"), by="rn")
  pDT[, value.norm := scale(value), by="variable"]
  # filter what to show?
  # pDT[value.norm > 2][,.(UMAP1 = sd(UMAP1), UMAP2 = sd(UMAP2)), by="variable"][,.(mean(UMAP1, UMAP2)), by="variable"][order(V1)]

  ggplot(pDT, aes(x=UMAP1, y=UMAP2)) +
    stat_summary_hex(aes(z=value),fun=mean, bins=100) +
    #scale_fill_gradient2(low="blue", midpoint = 0, high="red") +
    scale_fill_hexbin() +
    theme_bw(12) +
    facet_wrap(~variable) +
    ggtitle(paste("Marker Signatures", mnam))
  ggsave(out("Markers_Signatures_",mnam,"_UMAP_raw.pdf"), w=12,h=12)

  ggplot(pDT, aes(x=UMAP1, y=UMAP2)) +
    stat_summary_hex(aes(z=pmin(3, value.norm)),fun=mean, bins=100) +
    scale_fill_gradient2(low="blue", midpoint = 0, high="red") +
    #scale_fill_hexbin() +
    theme_bw(12) +
    facet_wrap(~variable) +
    ggtitle(paste("Marker Signatures", mnam))
  ggsave(out("Markers_Signatures_",mnam,"_UMAP_scaled.pdf"), w=12,h=12)

  cleanDev(); pdf(out("Markers_Signatures_",mnam,"_Clusters.pdf"),w=8,h=6)
  pheatmap(sapply(with(ann, split(rn, Clusters)), function(cx) colMeans(mx[cx,,drop=F])))
  dev.off()
}


# ANTIBODIES --------------------------------------------------------------
if("DM_CITEseq-2_NA_NM_1" %in% ann$sample){

  abMT <- citeseq.MT$`Antibody Capture`
  abMT <- SCRNA.TPXToLog(SCRNA.RawToTPX(abMT, scale.factor = 1e6))
  ann.c1 <- ann[sample == "DM_CITEseq-2_NA_NM_1"]

  # UMAP
  res <- data.table()
  abx <- row.names(abMT)[1]
  for(abx in row.names(abMT)){
    pDT <- copy(ann.c1)
    pDT$Signal <- abMT[abx, gsub("\\-1.+", "-1", ann.c1$rn)]
    pDT$Antibody <- abx
    res <- rbind(res, pDT)
  }
  res[,Signal.norm := scale(Signal), by="Antibody"]
  ggplot(res, aes(x=UMAP1, y=UMAP2)) +
    stat_summary_hex(bins = 100, aes(z=Signal.norm),fun=mean) +
    scale_fill_gradient2(low="blue", high="red") +
    facet_wrap(~Antibody) +
    theme_bw(12)
  ggsave(out("Antibodies_UMAP.pdf"), w=12+2, h=9+1)

  write.tsv(res, out("Antibodies.tsv"))

  # Correlations
  cMT <- corS(t(as.matrix(abMT)), use="pairwise.complete.obs")
  diag(cMT) <- NA
  dist <- as.dist(1-cMT)
  cleanDev(); pdf(out("Antibodies_Correlation.pdf"), w=5,h=4)
  pheatmap(cMT,
           clustering_distance_rows = dist,
           clustering_distance_cols = dist,
           breaks=seq(-1,1, 0.01), color=COLORS.HM.FUNC(200),
  )
  dev.off()

  # Percentile plots
  probx=0.9
  res[,percentile := quantile(Signal.norm, probs=probx, na.rm=TRUE), by="Antibody"]
  resN <- res[, sum(Signal.norm > percentile), by=c("Clusters", "Antibody")]
  resN[,clSize := sum(V1), by="Clusters"]
  stopifnot(all(resN[, length(unique(clSize)), by="Clusters"]$V1 == 1))
  resN[,percentage := V1/clSize*100]
  resN <- hierarch.ordering(resN, toOrder = "Clusters", orderBy = "Antibody", value.var = "percentage")
  resN <- hierarch.ordering(resN, toOrder = "Antibody", orderBy = "Clusters", value.var = "percentage")
  resN[is.na(percentage), percentage := 0]
  ggplot(resN, aes(x=Clusters, y=Antibody, fill=percentage)) +
    geom_tile() +
    ggtitle("Percent of cell in Clusters\nthat are above 90th percentile of antibody signal") +
    scale_fill_gradient(low="white", high="blue")
  ggsave(out("Antibodies_Percentile.pdf"), w=5, h=4)
}


# GUIDES ------------------------------------------------------

# . Number of guides --------------------------------------------------------
# Number of cells per guide / mixscape
ggplot(ann[!is.na(mixscape_class.global)], aes(x=guide, fill=mixscape_class.global)) +
  geom_bar(position="dodge") +
  facet_wrap(~sample, scales = "free") +
  theme_bw(12) + xRot() +
  ylab("Cells")
ggsave(out("Guides_Counts.pdf"), w=12,h=10)
# % of cells assigned
ggplot(ann[,sum(!is.na(guide))/.N*100, by="sample"], aes(x=sample, y=V1)) +
  geom_bar(stat="identity") +
  theme_bw(12) + xRot() +
  ylab("Percent of cells with assigned guide (NP or KO)")
ggsave(out("Guides_Counts_PercAssigned.pdf"), w=6,h=4)

# . Guides in UMAP ----------------------------------------------------------
pDT <- ann[!is.na(guide)]
pDT[, gene := gsub("_.+","", guide)]
(gx <- pDT[gene != "NTC"]$gene[1])
for(gx in unique(pDT[gene != "NTC"]$gene)){
  xDT <- pDT[gene %in% c(gx, "NTC")]
  cols <- length(unique(paste(xDT$CRISPR_Cellranger, xDT$mixscape_class.global)))
  rows <- length(unique(xDT$sample_broad))
  ggplot(xDT, aes(x=UMAP1, y=UMAP2)) + 
    geom_hex(bins=50) +
    scale_fill_gradient(low="lightgrey", high="blue") +
    facet_grid(sample_broad~CRISPR_Cellranger + mixscape_class.global) +
    geom_text(data=xDT[,.N, by=c("sample_broad", "CRISPR_Cellranger", "mixscape_class.global")],
              aes(label=N), x=0, y=0) +
    theme_bw(12)
  ggsave(out("Guides_UMAP_", gx, "_samplebroad.pdf"), w=cols * 2 + 2,h=rows * 2 + 1, limitsize = FALSE)
  
  cols <- length(unique(paste(xDT$CRISPR_Cellranger, xDT$mixscape_class.global)))
  rows <- length(unique(xDT$sample))
  ggplot(xDT, aes(x=UMAP1, y=UMAP2)) + 
    geom_hex(bins=50) +
    scale_fill_gradient(low="lightgrey", high="blue") +
    facet_grid(sample~CRISPR_Cellranger + mixscape_class.global) +
    geom_text(data=xDT[,.N, by=c("sample", "CRISPR_Cellranger", "mixscape_class.global")],
              aes(label=N), x=0, y=0) +
    theme_bw(12)
  ggsave(out("Guides_UMAP_", gx, "_sample.pdf"), w=cols * 2 + 2,h=rows * 2 + 1, limitsize = FALSE)
}



# . Guides on UMAP - NTCs ---------------------------------------------
ggplot(ann[mixscape_class == "NTC"], aes(x=UMAP1, y=UMAP2)) +
  geom_hex(bins=50) +
  theme_bw(12) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~ sample + CRISPR_Cellranger)
ggsave(out("Guides_UMAP_", "NTCs", ".pdf"), w=20,h=20)


# . Guides in cluster - NTCs ----------------------------------------------
res <- data.table()
pDT1 <- ann[mixscape_class == "NTC"]
sx <- pDT1$sample_broad[1]
for(sx in unique(pDT1$sample_broad)){
  pDT2 <- pDT1[sample_broad == sx]
  for(gx in unique(pDT2$CRISPR_Cellranger)){
    for(cx in unique(pDT2$Clusters)){
      #message(gx, "-", cx)
      mx <- as.matrix(with(pDT2, table(Clusters == cx, CRISPR_Cellranger == gx)))
      #print(mx)
      if(all(dim(mx) == c(2,2))){
        fish <- fisher.test(mx)
        res <- rbind(res, data.table(Clusters=cx, mixscape_class=gx, p=fish$p.value, OR=fish$estimate, sample=sx, total.cells=sum(mx)))
      }
    }
  }
}
res[,padj := p.adjust(p, method="BH")]
res[, log2OR := log2(pmin(5, OR + min(res[OR != 0]$OR)))]
res[,grp := paste(mixscape_class, sample)]
res <- hierarch.ordering(res, toOrder = "grp", orderBy = "Clusters", value.var = "log2OR")
#res <- hierarch.ordering(res, toOrder = "Clusters", orderBy = "grp", value.var = "log2OR")
ggplot(res, aes(
  x=Clusters,
  y=mixscape_class,
  color=log2OR,
  size=pmin(-log10(padj), 5))) +
  geom_point(shape=16) +
  scale_color_gradient2(name="log2OR", low="blue", high="red") +
  scale_size_continuous(name="padj") +
  facet_grid(sample ~ ., space = "free", scales = "free") +
  theme_bw(12) +
  theme(strip.text.y = element_text(angle=0)) +
  xRot()
ggsave(out("Guides_Fisher_NTCs.pdf"), w=10, h=length(unique(res$grp)) * 0.25 + 1, limitsize = FALSE)


# . Guides per cluster - fisher test - MIXSCAPE ------------------------------------------------------

# FIRST PREPARE DIFFERENT ANALYSIS (in different input tables)
fish.test.sets <- list()
x <- ann[!is.na(mixscape_class)][mixscape_class.global != "NP"]
# basic analysis of everything
fish.test.sets[[paste0("basic")]] <- x
# without b cells
fish.test.sets[[paste0("noBcells")]] <- x[!grepl("B.cell", Clusters)]
# Early branching analysis
eba <- list(MEP=c("MEP (early)", "MEP"),GMP=c("GMP (early)", "Mono", "GMP"))
xxx <- x[Clusters %in% do.call(c, eba)]
for(xnam in names(eba)){xxx[Clusters %in% eba[[xnam]], Clusters := xnam]}
fish.test.sets[["earlyBranches"]] <- xxx
# remove those with 0 rows
fish.test.sets <- fish.test.sets[sapply(fish.test.sets, nrow) > 0]

# START ANALYSIS
(fish.test.x <- names(fish.test.sets)[1])
for(fish.test.x in names(fish.test.sets)){
  res <- data.table()
  pDT1 <- fish.test.sets[[fish.test.x]]
  #pDT1 <- pDT1[!grepl("B.cell", Clusters) & !grepl("CLP", Clusters)]
  stopifnot(sum(is.na(pDT1$CRISPR_Cellranger)) == 0)
  sx <- pDT1$sample_broad[1]
  for(sx in unique(pDT1$sample_broad)){
    pDT2 <- pDT1[sample_broad == sx]
    gx <- pDT2$CRISPR_Cellranger[1]
    for(gx in unique(pDT2[mixscape_class != "NTC"]$CRISPR_Cellranger)){
      for(cx in unique(pDT2$Clusters)){
        for(ntc in unique(pDT2[mixscape_class == "NTC"][,.N, by="CRISPR_Cellranger"][N > 20]$CRISPR_Cellranger)){
          mx <- as.matrix(with(pDT2[CRISPR_Cellranger %in% c(gx, ntc)], table(Clusters == cx, CRISPR_Cellranger == gx)))
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
              guide.cells=nrow(pDT2[CRISPR_Cellranger  == gx])
            ))
          }
        }
      }
    }
  }
  
  # modify and save results
  res[,padj := p.adjust(p, method="BH")]
  res[, log2OR := log2(pmin(5, OR + min(res[OR != 0]$OR)))]
  res[,grp := paste(mixscape_class, sample)]
  write.tsv(res[,-c("grp"), with=F], out("Guides_Fisher_Mixscape_",fish.test.x,".tsv"))
  
  # Plot
  if(nrow(res) < 3) next
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


# Guides Signature differential analysis ----------------------------------

# Calculate stats
mMT <- lapply(names(marker.signatures), function(xnam){x <- marker.signatures[[xnam]]; colnames(x) <- paste(xnam, colnames(x)); x})
mMT <- do.call(cbind, mMT)
mMT <- scale(mMT)
#mMT <- cbind(marker.signatures$Larry, marker.signatures$PanglaoDB[,c("Basophils", "Megakaryocytes")])

sx <- ann$sample_broad[13]
registerDoMC(cores=8)
resSDA <- foreach(sx = unique(ann$sample_broad)) %dopar% {
  resSDA <- data.table()
  annS <- ann[sample_broad == sx]
  if(nrow(annS[mixscape_class.global == "NTC"]) < 10) return(data.table())
  guides <- unique(annS[mixscape_class.global == "KO"][,.N, by="guide"][N > 10]$guide)
  sigx <- colnames(mMT)[1]
  for(sigx in colnames(mMT)){
    xNTC <- mMT[annS[mixscape_class.global == "NTC"]$rn, sigx]
    (guidex <- guides[1])
    for(guidex in guides){
      x <- mMT[annS[guide == guidex & mixscape_class.global == "KO"]$rn, sigx]
      resSDA <- rbind(resSDA, data.table(
        sample=sx,
        guide=guidex,
        sig=sigx,
        p=wilcox.test(x, xNTC)$p.value,
        d=median(x) - median(xNTC)
      ))
    }
  }
  return(resSDA)
}
resSDA <- do.call(rbind, resSDA)
resSDA[, padj := p.adjust(p, method="BH")]
resSDA[, gene := gsub("_.+", "", guide)]
write.tsv(resSDA, out("SigDA.tsv"))

# Plot stats
ggplot(resSDA, aes(x=paste(sample, guide),y=sig, color=d, size=pmin(5, -log10(padj)))) +
  theme_bw(12) +
  geom_point() +
  scale_color_gradient2(low="blue", high="red") +
  facet_grid(. ~ gene, scales = "free", space = "free") +
  xRot()
ggsave(out("SigDA_Stats.pdf"), w=29,h=20)


# Signatures on UMAP
pDT <- ann[mixscape_class.global %in% c("KO", "NTC") & guide %in% c(guides, "NTC")]
grps <- length(unique(pDT$guide))
ggplot(pDT, aes(x=UMAP1, y=UMAP2)) +
  geom_hex(bins=50) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~guide, ncol = 7) +
  theme_bw(12)
ggsave(out("SigDA_UMAP.pdf"), w=7*4,h=ceiling(grps/7)*4+1, limitsize = FALSE)

