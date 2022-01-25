source("src/00_init.R")

base.dir <- "FULLINT_10_04_ProjectUMAPs/"
out <- dirout(base.dir)


# Folders -----------------------------------------------------------------
list.files(dirout_load("")(""))
inDir.funcs <- list(
  "in.vivo"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vivo"),
  "in.vitro"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vitro"),
  "leukemia"=dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")
)


# Read in vivo data and perform differnetial expression -------------------
mobjs <- list()
for(tissuex in c("in.vitro", "leukemia", "in.vivo")){
  (load(inDir.funcs[[tissuex]]("MonocleObject.RData")))
  mobjs[[tissuex]] <- monocle.obj
}


x <- mobjs$in.vivo

as.Seurat.NF <- function(x){
  logcounts(x) <- counts(x)
  as.Seurat(x)
}

ref <- as.Seurat.NF(mobjs$in.vivo)
query <- as.Seurat.NF(mobjs$in.vitro)

proj <- make.projection(
  query = query,
  ref = ref,
  filter.cells = FALSE
)