source("src/00_init.R")

base.dir <- "SCRNA_41_01_GSEA/"
out <- dirout(base.dir)

require(fgsea)


# EnrichR -----------------------------------------------------------------
(load(PATHS$RESOURCES$Enrichr.mouse))


# Load data ---------------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_32_DE_Nebula_simple")(""), pattern="DEG_Results_all.tsv", full.names = TRUE, recursive = TRUE)
names(ff) <- gsub("^.+\\/", "", dirname(ff))
ff <- lapply(ff, fread)
DE.RES <- rbindlist(ff, idcol = "tissue")
stopifnot(all(grepl("^GuideDE", DE.RES$term)))
DE.RES[abs(estimate) > 15, estimate := min(15,  abs(estimate)) * sign(estimate)]



# fgsea -------------------------------------------------------------------
gsea.res <- data.table()
de.grp <- DE.RES$guide[1]
for(de.grp in unique(DE.RES$guide)){
  for(dbx in names(enr.terms)){
    gsea.res <- rbind(gsea.res, data.table(fgsea(
      pathways=enr.terms[[dbx]], 
      stats=with(DE.RES[guide == de.grp], setNames(estimate, nm=gene_id))), 
    grp=de.grp,
    db=dbx))
  }
}
saveRDS(gsea.res, file=out("FGSEA.RDS"))