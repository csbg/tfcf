source("src/00_init.R")

base.dir <- "SCRNA_41_01_GSEA/"
out <- dirout(base.dir)

require(fgsea)


# EnrichR -----------------------------------------------------------------
(load(PATHS$RESOURCES$Enrichr.mouse))


# Load data ---------------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_33_DE_Nebula_testClustering")("in.vivo_14d_noClusters"), pattern="DEG_Results_all.tsv", full.names = TRUE, recursive = TRUE)
names(ff) <- gsub("^.+\\/", "", dirname(ff))
ff <- lapply(ff, fread)
DE.RES <- rbindlist(ff, idcol = "tissue")
stopifnot(all(grepl("^GuideDE", DE.RES$term)))
DE.RES[abs(estimate) > 15, estimate := min(15,  abs(estimate)) * sign(estimate)]



# fgsea -------------------------------------------------------------------
gsea.res <- data.table()
de.grp <- DE.RES$guide[1]
for(tx in unique(DE.RES$tissue)){
  for(de.grp in unique(DE.RES[tissue == tx]$guide)){
    for(dbx in names(enr.terms)){
      gsea.res <- rbind(gsea.res, data.table(fgsea(
        pathways=enr.terms[[dbx]], 
        stats=with(DE.RES[tissue == tx][guide == de.grp], setNames(estimate, nm=gene_id))), 
      grp=de.grp,
      tissue=tx,
      db=dbx))
    }
  }
}
saveRDS(gsea.res, file=out("FGSEA.RDS"))



# Plots -------------------------------------------------------------------
if(!"gsea.res" %in% ls()) gsea.res <- readRDS(out("FGSEA.RDS"))

# cleanup / export results
gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge, function(vec) paste(vec[1:10], collapse = ","))
write.tsv(gsea.res.export, out("GSEA_significant",".tsv"))

for(dbx in unique(gsea.res$db)){
  write.tsv(gsea.res.export[db == dbx], out("GSEA_significant_",dbx,".tsv"))
}


# Prepare for plotting
for(dbx in unique(gsea.res$db)){
  pDT <- gsea.res[db == dbx]
  pw.display <- unique(pDT[padj < 0.05][order(NES)][, tail(.SD, n=5), by=c("grp", "tissue")]$pathway)
  pDT <- pDT[pathway %in% pw.display]
  pDT <- hierarch.ordering(pDT, "pathway", "grp", "NES", TRUE)
  pDT <- hierarch.ordering(pDT, "grp", "pathway", "NES", TRUE)
  ggplot(pDT, aes(x=grp, y=pathway, color=NES, size=pmin(5, -log10(padj)))) +
    geom_point() + scale_color_gradient2(low="blue", mid="white", high="red") +
    geom_point(data=pDT[padj < 0.05], shape=1, color="black") +
    scale_size_continuous(range=c(0,5), limits = c(0,5)) +
    theme_bw(12) +
    xRot() +
    facet_grid(. ~ tissue, space="free", scales="free") +
    theme(strip.text.y=element_text(angle=0))
  ggsave(out("GSEA_plot_",dbx,".pdf"), w=30,h=length(unique(pDT$pathway)) * 0.2 + 3, limitsize = FALSE)
}


