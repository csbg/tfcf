source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("ECCITE4_04_Integration/")

require(Seurat)
require(limma)
require(edgeR)

# Read cellranger analysis results --------------------------------------------


AGG.CSV <- fread(paste(Sys.getenv("DATA"), "ECCITE4_INT", "outs", "aggregation.csv", sep="/"))
AGG.CSV$i <- 1:nrow(AGG.CSV)
AGG.CSV[,sample_id := paste0("ECCITE4_", sample_id)]


# Read data  and Seurat integration --------------------------------------------
seurat.file <- out("SeuratObject.RData")
if(!file.exists(seurat.file)){
  
  seurat.list <- list()
  additional.info <- list()
  
  for(dsx in rev(c("ECCITE4_WT", "ECCITE4_Cas9"))){
    print("---------------")
    print(dsx)
    print("---------------")
    path <- paste(Sys.getenv("DATA"), dsx, "outs", "filtered_feature_bc_matrix.h5", sep = "/")
    if(!file.exists(path)){
      message("File for ", dsx, " not found")
      next
    }
    data <- Read10X_h5(path)
    data.gx <- NA
    if(class(data) == "dgCMatrix"){
      data.gx <- data
    }
    if(is.list(data) && "Gene Expression" %in% names(data)){
      data.gx <- data[["Gene Expression"]]
      additional.info[[dsx]] <- data[names(data) != "Gene Expression"]
    }
    if(class(data.gx) != "dgCMatrix"){
      message("Data for ", dsx, " not in the right format")
    }
    seurat.obj <- CreateSeuratObject(counts = data.gx, project = dsx, min.cells = 5)
    seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
    seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
    seurat.obj <- NormalizeData(seurat.obj, verbose = FALSE)
    seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
    seurat.list[[dsx]] <- seurat.obj
  }
  
  # Integrate
  anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:20)
  rm(list=c("seurat.list", "data.gx", "data"))
  sobj <- IntegrateData(anchorset = anchors, dims = 1:20)
  rm(list="anchors")
  
  DefaultAssay(sobj) <- "integrated"
  sobj@meta.data$dataset <- sobj@meta.data$orig.ident
  
  # Cell Cycle scoring
  sobj <- CellCycleScoring(sobj, 
                           s.features = cc.genes$s.genes, 
                           g2m.features = cc.genes$g2m.genes, 
                           set.ident = TRUE)
  
  # Process full dataset
  sobj <- ScaleData(sobj, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 30, verbose = FALSE)
  sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:20)
  sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:20)
  sobj <- FindClusters(sobj, resolution = 0.5)
  
  # Store full dataset
  sobj <- SCRNA.DietSeurat(sobj)
  save(sobj, additional.info, file=seurat.file)
} else {
  print("Loading Seurat file")
  load(seurat.file)
}
sobj <- SCRNA.UndietSeurat(sobj)



# Cluster markers ---------------------------------------------------------
file.markers <- out("ClusterMarkers.tsv")
if(!file.exists(file.markers)){
  sobj <- ScaleData(sobj, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
  cl.markers <- FindAllMarkers(sobj, assay = "RNA", only.pos = TRUE, min.cells.feature=10, min.pct = 0.2)
  cl.markers <- data.table(cl.markers, keep.rownames = TRUE)
  write.table(cl.markers, file.markers)
  
  pDT <- cl.markers[gene %in% cl.markers[p_val_adj < 0.05][order(avg_log2FC, decreasing = TRUE)][,head(.SD, n=10), by="cluster"]$gene]
  pDT <- hierarch.ordering(pDT, toOrder = "gene", orderBy = "cluster", value.var = "avg_log2FC")
  #pDT <- hierarch.ordering(pDT, toOrder = "cluster", orderBy = "gene", value.var = "avg_log2FC")
  ggplot(pDT, aes(y=gene, x=factor(cluster), color=avg_log2FC, size=-log10(p_val_adj + 1e-5))) + 
    geom_point() +
    scale_color_gradient2(low="blue", high="red") +
    theme_bw(12)
  ggsave(out("ClusterMarkers.pdf"), w=12, h=20)
}




# Collect metadata --------------------------------------------------------------
ann <- data.table(sobj[["umap"]]@cell.embeddings, keep.rownames = TRUE)
ann <- merge(ann, data.table(sobj@meta.data, keep.rownames = TRUE), by="rn")
ann[, Barcode := gsub("_\\d+$", "", rn)]
#ann <- merge(ann, clusters, by=c("Barcode", "dataset"))
colnames(ann) <- gsub("_", ".", colnames(ann))
ann[,cluster := as.numeric(as.character(seurat.clusters))]
ann$seurat.clusters <- NULL
ann$integrated.snn.res.0.5 <- NULL
ann$orig.ident <- NULL
ann$old.ident <- NULL


# Add extra meta data (antibodies, guides) --------------------------------
ann.orig <- copy(ann)
cutoff <- 5
dsx <- names(additional.info)[1]
for(dsx in names(additional.info)){
  print(dsx)
  tx <- names(additional.info[[dsx]])[1]
  for(tx in names(additional.info[[dsx]])){
    txn <- make.names(tx)
    print(tx)
    x <- additional.info[[dsx]][[tx]]
    x.max <- apply(x, 2, max)
    x.sum <- apply(x, 2, sum)
    x <- as.matrix(x[,x.max/x.sum > 0.75 & x.sum > 20])
    ii <- apply(x, 2, function(col) which(col >= max(col))) # for each cell (column) get the rows that are above the cutoff
    #labelsx <- sapply(ii, function(x) paste(sort(names(x)), collapse=","))
    labelsx <- sapply(ii, function(ix) row.names(x)[ix])
    labelsx <- setNames(data.table(data.frame(labelsx), keep.rownames = TRUE), c("Barcode", txn))
    prep.ann <- if(txn %in% colnames(ann)) ann[, -txn,with=F] else ann
    ann <- rbind(ann[dataset != dsx], merge(prep.ann[dataset == dsx], labelsx, by="Barcode", all.x=TRUE), fill=TRUE)
  }
}
ann[, Guide := CRISPR.Guide.Capture]
ann[, Guide := gsub("_.+$", "", Guide)]
ann[grepl(",", Guide), Guide := NA]


# for(dsx in unique(ann$dataset)){
#   if(grepl("CITESEQ", dsx)){
#     stopifnot(nrow(ann[dataset == dsx & !is.na(CRISPR.Guide.Capture)]) == 0)
#     stopifnot(nrow(ann[dataset == dsx & !is.na(Antibody.Capture)]) != 0)
#   } else {
#     stopifnot(nrow(ann[dataset == dsx & !is.na(CRISPR.Guide.Capture)]) != 0)
#     stopifnot(nrow(ann[dataset == dsx & !is.na(Antibody.Capture)]) == 0)
#   }
# }

write.tsv(ann, out("Metadata_1_Original.tsv"))


# SETUP ENDS HERE ---------------------------------------------------------



# PLOT DATASET AND CELL CYCLE -------------------------------------------------------------------
ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_point(aes(color=dataset), alpha=1, size=0.1) +
  theme_bw(12)
ggsave(out("UMAP_datasets.pdf"), w=6, h=5)

ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_point(aes(color=Phase), alpha=1, size=0.1) +
  theme_bw(12)
ggsave(out("UMAP_CellCycle.pdf"), w=6, h=5)

# PLOT INTEGRATED CLUSTERS ------------------------------------------------
ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_point(aes(color=factor(cluster))) +
  facet_wrap( ~ dataset, ncol = 2) +
  theme_bw(12) +
  geom_label(data=ann[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by=c("cluster", "dataset")], aes(label=cluster))
ggsave(out("UMAP_clusters_points.pdf"), w=11, h=10)

ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap( ~ dataset, ncol = 2) +
  theme_bw(12) +
  geom_label(data=ann[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by=c("cluster", "dataset")], aes(label=cluster))
ggsave(out("UMAP_cluster.pdf"), w=11, h=10)



# Plot marker genes -------------------------------------------------------
str(bulk.data <- as.matrix(read.table("metadata/INVITRO_BulkExpression.csv", row.names = 1, sep=",", header = TRUE)))
# datG <- calcNormFactors(DGEList(bulk.data))
# desMat <- data.frame(row.names = colnames(bulk.data), group=gsub(".\\d$", "", colnames(bulk.data)))
# desMat$group <- factor(desMat$group, levels=c("LSK", "MEP", "GMP"))
# desMat <- model.matrix(~group, data = desMat)
# voomRes <- voom(datG, design = desMat,plot=FALSE)
# fit <- lmFit(voomRes, design=desMat)
# fit <- eBayes(fit)
# res <- data.table()
# for(coefx in colnames(coef(fit))){res <- rbind(res, data.table(topTable(fit, number = nrow(voomRes), coef=coefx), coef=coefx, keep.rownames = TRUE))}
# ll <- with(res[coef != "(Intercept)"][adj.P.Val < 0.2], split(rn, paste(coef, ifelse(logFC > 0, "up", "down"))))
# ll <- gplots::venn(ll, intersections = TRUE)
# ii <- attr(ll, "intersections")
# markers <- list(
#   GMP=ii[["groupGMP up"]],
#   MEP=ii[["groupMEP up"]],
#   LSK=ii[["groupGMP down:groupMEP down"]])
str(bulk.data <- as.matrix(read.table("metadata/INVITRO_BulkExpression.csv", row.names = 1, sep=",", header = TRUE)))
bulk.data <- sapply(split(colnames(bulk.data), gsub(".\\d$", "", colnames(bulk.data))), function(sx){rowMeans(bulk.data[,sx,drop=F])})
bulk.data <- bulk.data[apply(bulk.data, 1, max) > median(apply(bulk.data, 1, max)),]
bulk.data.fc <- sapply(colnames(bulk.data), function(cx){
  bulk.data[,cx] / apply(bulk.data[,colnames(bulk.data) != cx], 1, max)
})
pheatmap(corS(bulk.data))
markers <- lapply(apply(bulk.data.fc, 2, function(cx) which(cx > 3)), names)

ann.markers <- copy(ann)
mx <- names(markers)[1]
for(mx in names(markers)){
  message(mx)
  m <- sobj@assays$RNA@data[,ann.markers$rn]
  m <- m[intersect(markers[[mx]], row.names(m)),]
  m <- m[apply(m != 0, 1, sum) > ncol(m)/50,]
  print(dim(m))
  ann.markers[[paste0("Signature_", mx)]] <- rowMeans(scale(t(m)))
}

ggplot(melt(ann.markers, id.vars = c("rn", "UMAP.1", "UMAP.2"), measure.vars = grep("Signature_", colnames(ann.markers), value = TRUE)),
       aes(x=UMAP.1, UMAP.2, z=value)) +
  stat_summary_hex(fun = mean) +
  scale_fill_gradient2(low="blue", high="red") +
  facet_grid(. ~ variable) +
  theme_bw(12)
ggsave(out("UMAP_celltype_signatures.pdf"), w=13, h=4)


# Remove bad cells --------------------------------------------------------
for(qcm in c("percent.mt", "nFeature.RNA", "nCount.RNA")){
  print(qcm)
  ann$measure <- ann[[qcm]]
  p <- ggplot(ann, aes(y=measure + 0.1, x=factor(cluster))) + 
    geom_violin(color=NA, fill="lightblue") + 
    geom_boxplot(fill=NA, coef=Inf) +
    scale_y_log10() +
    theme_bw(12) +
    ylab(qcm)
    ggtitle(qcm)
  ggsave(out("QC_", qcm, ".pdf"), w=5,h=4, plot=p)
}

ann[,cluster.qual.keep :=TRUE]
ann[cluster %in% ann[,median(percent.mt), by="cluster"][V1 < 1]$cluster, cluster.qual.keep := FALSE]
write.tsv(ann, out("Metadata_2_Cleaned.tsv"))




# Now keep working with cleaned dataset -----------------------------------
ann <- ann[cluster.qual.keep == TRUE]


# 
# # ANTIBODIES --------------------------------------------------------------
# abMT <- additional.info$CITESEQ2$`Antibody Capture`
# abMT <- SCRNA.TPXToLog(SCRNA.RawToTPX(abMT, scale.factor = 1e6))
# ann.c1 <- ann[dataset == "CITESEQ2"]
# 
# # Plot signal
# res <- data.table()
# abx <- row.names(abMT)[1]
# for(abx in row.names(abMT)){
#   pDT <- copy(ann.c1)
#   pDT$Signal <- abMT[abx,ann.c1$Barcode]
#   pDT$Antibody <- abx
#   res <- rbind(res, pDT)
# }
# res[,Signal.norm := scale(Signal), by="Antibody"]
# ggplot(res, aes(x=UMAP.1, y=UMAP.2)) + 
#   stat_summary_hex(aes(z=Signal.norm),fun=mean) +
#   scale_fill_gradient2(low="blue", high="red") +
#   facet_wrap(~Antibody) +
#   theme_bw(12)
# ggsave(out("Antibodies_UMAP.pdf"), w=12+2, h=9+1)
# 
# cMT <- corS(t(as.matrix(abMT)), use="pairwise.complete.obs")
# diag(cMT) <- NA
# dist <- as.dist(1-cMT)
# cleanDev(); pdf(out("Antibodies_Correlation.pdf"), w=5,h=4)
# pheatmap(cMT,
#          clustering_distance_rows = dist,
#          clustering_distance_cols = dist,
#          breaks=seq(-1,1, 0.01), color=COLORS.HM.FUNC(200),
# )
# dev.off()
# 
# probx=0.9
# res[,percentile := quantile(Signal.norm, probs=probx, na.rm=TRUE), by="Antibody"]
# resN <- res[, sum(Signal.norm > percentile), by=c("cluster", "Antibody")]
# resN[,clSize := sum(V1), by="cluster"]
# stopifnot(all(resN[, length(unique(clSize)), by="cluster"]$V1 == 1))
# resN[,percentage := V1/clSize*100]
# resN <- hierarch.ordering(resN, toOrder = "cluster", orderBy = "Antibody", value.var = "percentage")
# resN <- hierarch.ordering(resN, toOrder = "Antibody", orderBy = "cluster", value.var = "percentage")
# resN[is.na(percentage), percentage := 0]
# ggplot(resN, aes(x=cluster, y=Antibody, fill=percentage)) + 
#   geom_tile() + 
#   ggtitle("Percent of cell in cluster\nthat are above 90th percentile of antibody signal") + 
#   scale_fill_gradient(low="white", high="blue")
# ggsave(out("Antibodies_Percentile.pdf"), w=5, h=4)



# MIXSCAPE Analysis -----------------------------------------------------
eccite.file <- out("Mixscape_sobj.RData")
if(!file.exists(eccite.file)){
  sobj2 <- sobj
  sobj2@meta.data$guide <- ann[match(row.names(sobj2@meta.data), rn)]$Guide
  with(data.table(sobj2@meta.data), table(guide, dataset))
  #with(g2, table(Guide))
  sobj2 <- sobj2[,!is.na(sobj2@meta.data$guide)]
  with(data.table(sobj2@meta.data), table(guide, dataset))
  eccite <- CalcPerturbSig(
    object = sobj2, assay = "RNA", slot = "data", 
    gd.class ="guide",nt.cell.class = "NTC",
    reduction = "pca",ndims = 20,num.neighbors = 20,
    new.assay.name = "PRTB")
  eccite <- ScaleData(object = eccite, assay = "PRTB", do.scale = F, do.center = T)
  eccite <- RunMixscape(
    object = eccite,assay = "PRTB",slot = "scale.data",
    labels = "guide",nt.class.name = "NTC",
    min.de.genes = 5,
    iter.num = 10,
    de.assay = "RNA",
    verbose = TRUE,
    prtb.type = "KO")
  
  table(eccite$mixscape_class.global, eccite$guide)
  prop.table(table(eccite$mixscape_class.global, eccite$guide),2)
  
  eccite.diet <- SCRNA.DietSeurat(sobj = eccite)
  save(eccite.diet, file=eccite.file)
} else {
  print("Loading Mixscape results")
  load(eccite.file)
}
eccite <- SCRNA.UndietSeurat(eccite.diet)

ann.orig <- fread(out("Metadata_2_Cleaned.tsv"))
ann <- merge(ann.orig, data.table(eccite@meta.data, keep.rownames = TRUE)[,c("rn", "mixscape_class", "mixscape_class_p_ko", "mixscape_class.global"),with=F], by="rn", all.x=TRUE)
write.tsv(ann, out("Metadata_3_Mixscape.tsv"))



# Export results for cellranger -------------------------------------------
if("AGG.CSV" %in% ls()){
  ann.exp <- fread(out("Metadata_3_Mixscape.tsv"))
  ann.exp <- merge(ann.exp, AGG.CSV[,c("sample_id", "i"),with=F], by.x="dataset", by.y="sample_id", all.x=TRUE)
  stopifnot(!any(is.na(ann.exp$i)))
  ann.exp[,Barcode := paste0(gsub("-.+$", "", Barcode), "-", i)]
  write.table(ann.exp[cluster.qual.keep == TRUE][,c("Barcode", "UMAP.1", "UMAP.2"),with=F], file=out("Cellranger_UMAP.csv"), sep=",", col.names = c("Barcode", "UMAP-1", "UMAP-2"), quote=F, row.names = F)
  write.table(ann.exp[cluster.qual.keep == TRUE][,c("Barcode", "mixscape_class"),with=F], file=out("Cellranger_MIXSCAPE.csv"), sep=",", col.names = c("Barcode", "MIXSCAPE"), quote=F, row.names = F)
  write.table(ann.exp[cluster.qual.keep == TRUE][,c("Barcode", "cluster"),with=F], file=out("Cellranger_Clusters.csv"), sep=",", col.names = c("Barcode", "Clusters_Seurat"), quote=F, row.names = F)
}



# PLOT GUIDES ------------------------------------------------
ann <- fread(out("Metadata_3_Mixscape.tsv"))
ggplot(ann[mixscape_class.global != "NP"], aes(x=UMAP.1, y=UMAP.2)) + 
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  #geom_point(data=, aes(color=mixscape_class)) +
  facet_wrap(~ gsub("_.+", "", mixscape_class)) +
  theme_bw(12)
ggsave(out("UMAP_Guides.pdf"), w=7, h=5)

ggplot(ann[mixscape_class.global == "NP"], aes(x=UMAP.1, y=UMAP.2)) + 
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  #geom_point(data=, aes(color=mixscape_class)) +
  facet_wrap(~ gsub("_.+", "", mixscape_class)) +
  theme_bw(12)
ggsave(out("UMAP_Guides_NP.pdf"), w=7, h=5)

ggplot(ann[mixscape_class.global != "NTC"], aes(x=mixscape_class, fill=dataset)) + geom_bar() + theme_bw(12) + xRot() 
ggsave(out("Guides_Barplot.pdf"), w=7,h=5)

# GUIDE ENRICHMENT
res <- data.table()
pDT <- ann[!is.na(mixscape_class)][mixscape_class.global != "NP"]
for(gx in unique(pDT[mixscape_class != "NTC"]$mixscape_class)){
  for(cx in unique(pDT$cluster)){
    message(gx, "-", cx)
    mx <- as.matrix(with(ann[mixscape_class %in% c(gx, "NTC")], table(cluster == cx, mixscape_class == gx)))
    print(mx)
    if(dim(mx) == c(2,2)){
      fish <- fisher.test(mx)
      res <- rbind(res, data.table(Cluster=cx, mixscape_class=gx, p=fish$p.value, OR=fish$estimate))
    }
  }
}
res[,padj := p.adjust(p, method="BH")]
res[padj < 0.05]
res[, log2OR := log2(pmin(5, OR + min(res[OR != 0]$OR)))]
res <- hierarch.ordering(res, toOrder = "mixscape_class", orderBy = "Cluster", value.var = "log2OR")
res <- hierarch.ordering(res, toOrder = "Cluster", orderBy = "mixscape_class", value.var = "log2OR")
ggplot(res, aes(
  x=Cluster,
  y=mixscape_class, 
  color=log2OR, 
  size=pmin(-log10(padj), 5))) + 
  geom_point(shape=16) +
  scale_color_gradient2(name="log2OR", low="blue", high="red") +
  scale_size_continuous(name="padj") + 
  theme_bw(12)
ggsave(out("Guides_Fisher.pdf"), w=5, h=4)




# 
# # Perform DE based on Mixscape --------------------------------------------
# 
# # # Global approach fast analysis --------------------------------------------
# # Idents(eccite) <- "mixscape_class"
# # res <- data.table()
# # for(guidex in unique(grep(" KO$", eccite@meta.data$mixscape_class, value = TRUE))){
# #   print(guidex)
# #   markers <- FindMarkers(
# #     eccite, ident.1 = guidex, ident.2 = "NTC", logfc.threshold = log(2),
# #     assay = "RNA",only.pos = FALSE, test.use="roc")
# #   res <- rbind(res, data.table(markers, guide=guidex, keep.rownames = TRUE))
# # }
# # markers.global.fast <- res
# # write.tsv(markers.global.fast, out("CRISPR.DE.global.fast.tsv"))
# # 
# # # Cluster-specific approach fast analysis --------------------------------------------
# # gx <- "Pu.1 KO"
# # Idents(eccite) <- "mixscape_class"
# # res <- data.table()
# # for(gx in unique(grep(" KO$", eccite@meta.data$mixscape_class, value = TRUE))){
# #   if(gx == "NTC") next
# #   cons.markers <- FindConservedMarkers(min.cells.group=0,
# #                                        object = eccite, ident.1=gx, ident.2 = "NTC", 
# #                                        assay = "RNA", grouping.var = "seurat_clusters",
# #                                        test.use="roc")
# #   
# #   x <- data.table(cons.markers, keep.rownames = TRUE)
# #   
# #   i <- 0
# #   for(i in unique(as.character(eccite@meta.data$seurat_clusters))){
# #     ret <- data.table(crispr=gx, x[,c("rn", grep(paste0("^", i, "_"), colnames(x), value=TRUE)),with=F], cluster=i)
# #     colnames(ret) <- gsub("^\\d+_", "", colnames(ret))
# #     res <- rbind(res, ret, fill=TRUE)
# #   }
# # }
# # markers.cluster.fast <- res[!is.na(avg_log2FC)]
# # write.tsv(markers.cluster.fast, out("CRISPR.DE.cluster.fast.tsv"))
# # 
# # 
# # # Collect genes that appear across analysis -------------------------------
# # markers.global.fast[,length(unique(rn)), by=c("guide")]
# # markers.cluster.fast[,length(unique(rn)), by=c("crispr")]
# # 
# # genes.to.test <- unique(markers.global.fast$rn, markers.cluster.fast$rn, )
# # # CANNOT USE THIS BECAUSE SEURAT IS BUGGY
# 
# # Global approach negbinom analysis --------------------------------------------
# Idents(eccite) <- "mixscape_class"
# res <- data.table()
# guidex <- "Pu.1 KO"
# for(guidex in unique(grep(" KO$", eccite@meta.data$mixscape_class, value = TRUE))){
#   print(guidex)
#   markers <- FindMarkers(
#     eccite, ident.1 = guidex, ident.2 = "NTC", 
#     #min.cells.group = 0,min.cells.feature = 0,
#     #logfc.threshold = 0, min.pct = 0, features = genes.to.test,
#     assay = "RNA",only.pos = FALSE, test.use="negbinom")
#   res <- rbind(res, data.table(markers, guide=guidex, keep.rownames = TRUE))
# }
# markers.global.negbinom <- res
# write.tsv(markers.global.negbinom, out("CRISPR.DE.global.negbinom.tsv"))
# 
# 
# # Confirm with global calculation of log fold changes - NOT WORKING BECAUSE SEURAT IS BUGGY ------------------------------------
# # me <- data.table(eccite@meta.data, keep.rownames = TRUE)
# # dt <- SCRNA.LogToTPX(eccite@assays$RNA@data)
# # 
# # gg <- markers.global.negbinom[guide == guidex]$rn
# # m1 <- rowMeans(dt[gg,me[mixscape_class == guidex]$rn])
# # m2 <- rowMeans(dt[gg,me[mixscape_class == "NTC"]$rn])
# # fc <- log2(m1+1) - log2(m2+1)
# # data.table(fc[1:5])
# # markers.global.negbinom[guide == guidex][1:5]
# # fc <- FoldChange(eccite,ident.1 = guidex, ident.2 = "NTC")
# # ii <- intersect(row.names(fc), markers.global.negbinom[guide == guidex]$rn)
# # plot(fc[ii,1], markers.global.negbinom[match(ii, rn)]$avg_log2FC)
# 
# 
# # Cluster-specific approach negbinom analysis --------------------------------------------
# gx <- "Pu.1 KO"
# Idents(eccite) <- "mixscape_class"
# res <- data.table()
# for(gx in unique(grep(" KO$", eccite@meta.data$mixscape_class, value = TRUE))){
#   if(gx == "NTC") next
#   cons.markers <- FindConservedMarkers(min.cells.group=0,
#                                        object = eccite, ident.1=gx, ident.2 = "NTC", 
#                                        assay = "RNA", grouping.var = "seurat_clusters",
#                                        test.use="negbinom")
#   
#   x <- data.table(cons.markers, keep.rownames = TRUE)
#   
#   i <- 0
#   for(i in unique(as.character(eccite@meta.data$seurat_clusters))){
#     ret <- data.table(crispr=gx, x[,c("rn", grep(paste0("^", i, "_"), colnames(x), value=TRUE)),with=F], cluster=i)
#     colnames(ret) <- gsub("^\\d+_", "", colnames(ret))
#     res <- rbind(res, ret, fill=TRUE)
#   }
# }
# 
# markers.cluster.negbinom <- res[!is.na(avg_log2FC)]
# write.tsv(markers.cluster.negbinom, out("CRISPR.DE.cluster.negbinom.tsv"))
# 
# 
# 
# # Confirm with global calculation of log fold changes - NOT WORKING BECAUSE SEURAT IS BUGGY ------------------------------------
# # me <- data.table(eccite@meta.data, keep.rownames = TRUE)
# # dt <- SCRNA.LogToTPX(eccite@assays$RNA@data)
# # 
# # gg <- markers.global.negbinom[guide == guidex]$rn
# # m1 <- rowMeans(dt[gg,me[mixscape_class == guidex]$rn])
# # m2 <- rowMeans(dt[gg,me[mixscape_class == "NTC"]$rn])
# # fc <- log2(m1+1) - log2(m2+1)
# # data.table(fc[1:5])
# # markers.global.negbinom[guide == guidex][1:5]
# 
# # ii <- intersect(row.names(fc), markers.global.negbinom[guide == guidex]$rn)
# # plot(fc[ii,1], markers.global.negbinom[match(ii, rn)]$avg_log2FC)
# 
# # gx <- "Kmt2d KO"
# # with(markers.cluster.negbinom[p_val_adj < 0.05 & crispr == gx][order(rn)], split(rn, cluster))
# # markers.global.negbinom[guide == gx & p_val_adj < 0.05][rn == "Coro1a"]
# # 
# # VlnPlot(eccite, features = "Coro1a", idents = c(gx, "NTC"), assay="RNA")
# # VlnPlot(eccite, features = "Coro1a", idents = c(gx, "NTC"), assay="RNA", group.by = "mix_cluster")
# # 
# # fc <- FoldChange(eccite,ident.1 = gx, ident.2 = "NTC", assay="RNA")
# # cor(fc[markers.global.negbinom[guide == gx]$rn,1], markers.global.negbinom[guide == gx]$avg_log2FC)
# # 
# # Idents(eccite) <- "mix_cluster"
# # x <- markers.cluster.negbinom[cluster == 1 & crispr == gx]
# # 
# # cor(fc[x$rn,1], x$avg_log2FC)
# # plot(fc[x$rn,1], x$avg_log2FC)
# 
# 
# 
# # Combine results with again calculated logFCs ----------------------------
# res <- data.table()
# gx <- "Pu.1 KO"
# cx <- 1
# for(gx in unique(markers.global.negbinom$guide)){
#   print(gx)
#   Idents(eccite) <- "mixscape_class"
#   fc <- data.table(Cluster = "All", FoldChange(eccite,ident.1 = gx, ident.2 = "NTC", assay="RNA", fc.name="fc"), keep.rownames = TRUE)
#   fc <- merge(fc, markers.global.negbinom[guide == gx][,c("rn", "p_val", "avg_log2FC", "p_val_adj")], by="rn", all=TRUE)
#   res <- rbind(res, data.table(fc, guide=gx))
#   
#   Idents(eccite) <- "mix_cluster"
#   for(cx in unique(markers.cluster.negbinom$cluster)){
#     print(cx)
#     gxn <- paste(gsub(" KO", "_KO", gx), cx)
#     if(sum(eccite@meta.data$mix_cluster == gxn) < 5) next
#     fc <- data.table(Cluster = cx, FoldChange(eccite,ident.1 = gxn, ident.2 = paste("NTC", cx), assay="RNA", fc.name="fc"), keep.rownames = TRUE)
#     fc <- merge(fc, markers.cluster.negbinom[crispr == gx & cluster == cx][,c("rn", "p_val", "avg_log2FC", "p_val_adj")], by="rn", all=TRUE)
#     res <- rbind(res, data.table(fc, guide=gx))
#   }
# }
# 
# res[,sig := sum(p_val_adj < 0.05, na.rm=TRUE), by=c("guide", "rn")]
# table(is.na(res$sig))
# 
# 
# res[,id := paste(guide, Cluster)]
# resMT <- toMT(res, row="rn", col="id", val = "fc")
# cMT <- corS(resMT, use="pairwise.complete.obs")
# diag(cMT) <- NA
# dd <- as.dist(1-cMT)
# # cleanDev(); pdf(out("Guides_DE_corHM.pdf"),w=12,h=12)
# cleanDev(); pdf("~/Guides_DE_corHM.pdf",w=12,h=12)
# pheatmap(cMT,
#          clustering_distance_rows = dd, clustering_distance_cols = dd,
#          breaks=seq(-1,1, 0.01), color=COLORS.HM.FUNC(200))
# dev.off()
# 
# pDT <- data.table()
# for(gx in unique(res$guide)){
#   #gx <- "Pu.1 KO"
#   resx <- res[guide == gx]
#   resMT <- toMT(resx, row="rn", col="Cluster", val = "fc")
#   delta <- resMT[, "All"] - apply(resMT[,colnames(resMT) != "All", drop=F], 1, max, na.rm=TRUE)
#   delta <- delta[2 * delta > resMT[, "All"]]
#   gg <- c(
#     head(sort(delta[names(delta) %in% resx[Cluster == "All"][fc > 0][p_val_adj < 0.05]$rn], decreasing = TRUE), 5)
#     #head(sort(delta[names(delta) %in% resx[Cluster == "All"][fc > 0][p_val_adj < 0.05]$rn], decreasing = FALSE), 5)
#   )
#   pDT <- rbind(pDT, resx[rn %in% names(gg)])
# }
# 
# 
# # pDT <- res[sig >= 1]
# # pDT <- pDT[order(abs(fc), decreasing = TRUE)][,head(.SD, n=10), by=c("Cluster", "guide")]
# 
# ggplot(pDT, aes(x=Cluster, y=rn, color=fc, size=abs(fc))) + 
#   geom_point() +
#   scale_size_continuous(range = c(0,5)) +
#   geom_point(data=pDT[p_val_adj < 0.05], shape=1, color="black") +
#   scale_color_gradient2(low="blue", mid="white", high="red") +
#   facet_grid(guide ~ ., space = "free", scales = "free") +
#   theme_bw(12) +
#   xlab("Cluster where logFC was calculated")
# ggsave("~/DE.res.pdf",w=4, h=12)
# 
# xDT <- unique(pDT[,c("rn", "guide"),with=F])
# res <- data.table()
# i <- 1
# for(i in 1:nrow(xDT)){
#   me <- data.table(eccite@meta.data, keep.rownames = TRUE)
#   me <- me[mixscape_class %in% c(xDT[i]$guide, "NTC")]
#   me$Expression <- eccite@assays$RNA@data[xDT[i]$rn, me$rn]
#   me$Gene <- xDT[i]$rn
#   me$Guide_violin <- xDT[i]$guide
#   res <- rbind(res, me)
# }
# ggplot(res, aes(x=seurat_clusters, y=Expression, fill= guide == "NTC", color= guide == "NTC")) + 
#   #geom_violin(color=NA) +
#   geom_boxplot(fill=NA, coef=1e20) + 
#   scale_color_manual(values=c("#33a02c", "#6a3d9a")) + 
#   scale_fill_manual(values=c("#33a02c", "#6a3d9a")) + 
#   theme(panel.grid = element_blank()) +
#   theme_bw(12) +
#   facet_wrap(~ Gene + Guide_violin, scales = "free_y")
# ggsave("~/DE.violin.res.pdf",w=16, h=12)
# 
# # # Create Dotplot of Mixscape results ----------------------------------------------------------
# # gg <- unique(c(
# #   res[order(avg_log2FC, decreasing = TRUE)][,head(.SD,n=5), by="guide"]$rn,
# #   res[order(avg_log2FC, decreasing = FALSE)][,head(.SD,n=5), by="guide"]$rn))
# # 
# # eccite@meta.data$mix_cluster <- with(eccite@meta.data, paste(gsub(" ", "_", mixscape_class), seurat_clusters))
# # Idents(eccite) <- "mixscape_class"
# # avExp <- AverageExpression(eccite, assays="RNA")
# # perExp <- sapply(
# #   with(data.table(eccite@meta.data, keep.rownames = TRUE), split(rn, mixscape_class)),
# #   function(cells){
# #     rowSums(eccite@assays$RNA@counts[,cells, drop = FALSE] != 0) / length(cells)
# #   }) * 100
# # 
# # dotplotDT <- merge(
# #   melt(data.table(avExp$RNA, keep.rownames = TRUE), id.vars = "rn", value.name = "average",variable.factor = FALSE),
# #   melt(data.table(perExp, keep.rownames = TRUE), id.vars = "rn", value.name = "percentage",variable.factor = FALSE),
# #   by=c("rn", "variable"))
# # 
# # dotplotDT <- merge(dotplotDT, dotplotDT[variable == "NTC"][,c("rn", "average")], by="rn", suffixes = c("", "_NTC"))
# # dotplotDT[, scaleExp := average - average_NTC]
# # dotplotDT[, scaleExp := scale(scaleExp, center = FALSE), by="rn"]
# # 
# # #dotplotDT <- cbind(dotplotDT, setNames(data.table(do.call(rbind, strsplit(dotplotDT$variable, " "))), c("guide", "cluster")))
# # markers.all.plot <- markers.all
# # markers.all.plot[, percentage := pct.1 * 100]
# # markers.all.plot[, variable := guide]
# # pDT <- dotplotDT[rn %in% gg]
# # pDT <- hierarch.ordering(pDT, toOrder = "variable", orderBy = "rn", value.var = "scaleExp")
# # pDT <- hierarch.ordering(pDT, toOrder = "rn", orderBy = "variable", value.var = "scaleExp")
# # ggplot(pDT, aes(x=variable, y=rn, size=percentage, color=scaleExp)) + 
# #   geom_point() + 
# #   geom_point(data=markers.all.plot[rn %in% gg & p_val_adj < 0.05], shape=1, color="black") + 
# #   scale_color_gradient2(low='blue', high = "red") +
# #   scale_size_continuous(limits=c(0, 100)) +
# #   theme_bw() +
# #   xRot()
# # ggsave(out("Mixscape_Guides.pdf"), w=10,h=12)
# # 
# # DoHeatmap(size = 3,
# #   eccite, features = gg, group.by = "mix_cluster", assay = "RNA", slot="data",
# #   cells=data.table(eccite@meta.data, keep.rownames = TRUE)[guide %in% c("NTC", "Pu.1")]$rn)
# # ggsave(out("Mixscape_HM.jpeg"), w=20,h=6, dpi = 100)
# # 
# # DotPlot(eccite, features = gg, assay="RNA") + xRot() + coord_flip()
# # ggsave(out("Mixscape_Guides_Seurat.pdf"), w=10,h=12)
# # 
# # markers.all.MT <- toMT(markers.all, row = "rn", col = "guide", val = "avg_log2FC")
# # cMT <- corS(markers.all.MT, use="pairwise.complete.obs")
# # diag(cMT) <- NA
# # dd <- as.dist(1-cMT)
# # cleanDev(); pdf(out("Mixscape_corHM.pdf"),w=4.5,h=4)
# # pheatmap(cMT,
# #          clustering_distance_rows = dd, clustering_distance_cols = dd,
# #          breaks=seq(-1,1, 0.01), color=COLORS.HM.FUNC(200))
# # dev.off()
# # 
# # 
# # # Mixscape - Second DE approach, by cluster ---------------------------------------------
# # stop("Would be nice to get the below consistent with the above (test same genes)")
# # 
# # write.tsv(res, out("Mixscape_byCluster.tsv"))
# # 
# # markers.broad <- fread(out("Mixscape_Guides.tsv"))
# # markers.byClu <- fread(out("Mixscape_byCluster_Guides.tsv"))
# # 
# # markers.broad[,crispr := guide]
# # markers.broad <- markers.broad[,-c("guide", "percentage", "variable"), with=F]
# # 
# # markers.comb <- rbind(markers.broad, markers.byClu, fill=TRUE)
# # markers.comb[, cluster := as.character(cluster)]
# # markers.comb[is.na(cluster), cluster := "Broad"]
# # markers.comb$cluster <- factor(markers.comb$cluster, levels=c("Broad", as.character(0:max(markers.byClu$cluster))))
# # 
# # for(gx in unique(markers.comb$crispr)){
# #   pDT <- markers.comb[crispr == gx]
# #   pDT <- pDT[rn %in% pDT[order(abs(avg_log2FC), decreasing = TRUE)][,head(.SD, n=5),by=c("cluster")]$rn]
# #   ggplot(pDT, aes(x=cluster, y=rn, color=avg_log2FC, size=pmin(5, -log10(p_val_adj)))) + 
# #     geom_point() + scale_color_gradient2(low="blue", high="red") +
# #     theme_bw(12) +
# #     ggtitle(gx)
# #   ggsave(out("Mixscape_byCluster_",gx,".pdf"), w=8,h=8)
# # }


