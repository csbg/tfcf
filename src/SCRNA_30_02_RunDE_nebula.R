source("src/00_init.R")

require(nebula)
require(doMC)
source("src/FUNC_Monocle_PLUS.R")

basedir <- "SCRNA_30_DE_Nebula/"
outB <- dirout(basedir)


# Sample annotation -------------------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)
TISSUES <- PATHS$SCRNA$MONOCLE.NAMES

# load datasets -----------------------------------------------------------
mobjs <- list()
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}


# Which analysis to take celltypes from? -----------------------------------
ANALYSIS <- "monocle.singleR"




# # Analysis of Myeloid cells ---------------------------------------------------------------
# ct.use <- "myeloid"
# (tissue.name <- names(mobjs)[1])
# for(tissue.name in names(mobjs)){
# 
#   # Monocle object
#   monocle.obj <- mobjs[[tissue.name]]
#   monocle.obj <- monocle.obj[, monocle.obj$timepoint != "28d"]
#   monocle.obj <- monocle.obj[, monocle.obj$sample != "WT-LSK_OP0_NM_7d_1"]
# 
#   # Filter only celltypes of interest
#   ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
#   ann <- ann[rn %in% colnames(monocle.obj)]
#   sort(unique(ann$Clusters))
#   for(x in c("Ery", "B.cell", "CLP", "MEP", "T.cell")){ann <- ann[!grepl(x, Clusters)]}
#   sort(unique(ann$Clusters))
# 
#   # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
#   monocle.obj <- monocle.obj[, ann$rn]
#   monocle.obj$clusterDE <- ann$Clusters
# 
#   out <- dirout(paste0(basedir, tissue.name, "_", ct.use))
# 
#   source("src/SCRNA_30_01_DE_FUNC_nebula.R")
# }
# 
# 
# # Analysis of Erythroid cells ---------------------------------------------------------------
# ct.use <- "erythroid"
# (tissue.name <- names(mobjs)[1])
# for(tissue.name in "in.vivo"){
# 
#   # Monocle object
#   monocle.obj <- mobjs[[tissue.name]]
#   monocle.obj <- monocle.obj[, monocle.obj$timepoint != "28d"]
#   monocle.obj <- monocle.obj[, monocle.obj$sample != "WT-LSK_OP0_NM_7d_1"]
# 
#   # Filter only celltypes of interest
#   ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
#   ann <- ann[rn %in% colnames(monocle.obj)]
#   ann.use <- data.table()
#   sort(unique(ann$Clusters))
#   for(x in c("Ery", "MEP")){ann.use <- rbind(ann.use, ann[!grepl(x, Clusters)])}
#   ann.use <- unique(ann.use)
#   ann <- ann.use
#   sort(unique(ann$Clusters))
# 
#   # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
#   monocle.obj <- monocle.obj[, ann$rn]
#   monocle.obj$clusterDE <- ann$Clusters
# 
#   out <- dirout(paste0(basedir, tissue.name, "_", ct.use))
# 
#   source("src/SCRNA_30_01_DE_FUNC_nebula.R")
# }
# 
# 
# 
# 
# # Analysis of All cells in vivo ---------------------------------------------------------------
# ct.use <- "everything"
# (tissue.name <- names(mobjs)[1])
# for(tissue.name in "in.vivo"){
# 
#   # Monocle object
#   monocle.obj <- mobjs[[tissue.name]]
#   monocle.obj <- monocle.obj[, monocle.obj$timepoint != "28d"]
#   monocle.obj <- monocle.obj[, monocle.obj$sample != "WT-LSK_OP0_NM_7d_1"]
# 
#   # Filter only celltypes of interest
#   ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
#   ann <- ann[rn %in% colnames(monocle.obj)]
#   sort(unique(ann$Clusters))
#   for(x in c("B.cell", "CLP", "T.cell")){ann <- ann[!grepl(x, Clusters)]}
#   sort(unique(ann$Clusters))
# 
#   # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
#   monocle.obj <- monocle.obj[, ann$rn]
#   monocle.obj$clusterDE <- ann$Clusters
# 
#   out <- dirout(paste0(basedir, tissue.name, "_", ct.use))
# 
#   source("src/SCRNA_30_01_DE_FUNC_nebula.R")
# }
# 
# 
# 
# # Analysis without Mixscape ---------------------------------------------------------------
# ct.use <- "noMixscape"
# (tissue.name <- names(mobjs)[1])
# for(tissue.name in c("ex.vivo", "leukemia")){
# 
#   # Monocle object
#   monocle.obj <- mobjs[[tissue.name]]
#   monocle.obj <- monocle.obj[, monocle.obj$timepoint != "28d"]
#   monocle.obj <- monocle.obj[, monocle.obj$sample != "WT-LSK_OP0_NM_7d_1"]
# 
#   # Filter only celltypes of interest
#   ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
#   ann <- ann[rn %in% colnames(monocle.obj)]
#   sort(unique(ann$Clusters))
#   for(x in c("Ery", "B.cell", "CLP", "MEP", "T.cell")){ann <- ann[!grepl(x, Clusters)]}
#   sort(unique(ann$Clusters))
# 
#   # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
#   monocle.obj <- monocle.obj[, ann$rn]
#   monocle.obj$clusterDE <- ann$Clusters
# 
#   # Remove Mixscape effects
#   monocle.obj$mixscape_class.global[monocle.obj$mixscape_class.global == "NP"] <- "KO"
# 
#   out <- dirout(paste0(basedir, tissue.name, "_", ct.use))
# 
#   source("src/SCRNA_30_01_DE_FUNC_nebula.R")
# }
# 
# 
# # Analysis of different time points ex.vivo ---------------------------------------------------------------
# tissue.name <- "ex.vivo"
# for(ct.use in unique(mobjs[[tissue.name]]$timepoint)){
#   
#   # Monocle object
#   monocle.obj <- mobjs[[tissue.name]]
#   monocle.obj <- monocle.obj[, monocle.obj$timepoint == ct.use]
#   monocle.obj <- monocle.obj[, monocle.obj$sample != "WT-LSK_OP0_NM_7d_1"]
#   
#   # Filter only celltypes of interest
#   ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
#   ann <- ann[rn %in% colnames(monocle.obj)]
#   sort(unique(ann$Clusters))
#   for(x in c("Ery", "B.cell", "CLP", "MEP", "T.cell")){ann <- ann[!grepl(x, Clusters)]}
#   sort(unique(ann$Clusters))
#   
#   # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
#   monocle.obj <- monocle.obj[, ann$rn]
#   monocle.obj$clusterDE <- ann$Clusters
#   
#   out <- dirout(paste0(basedir, tissue.name, "_", ct.use))
#   
#   source("src/SCRNA_30_01_DE_FUNC_nebula.R")
# }





# # Analysis of Lin- cells in vivo ---------------------------------------------------------------
# ct.use <- "linNeg"
# (tissue.name <- names(mobjs)[2])
# for(tissue.name in "in.vivo"){
#   
#   # Monocle object
#   monocle.obj <- mobjs[[tissue.name]]
#   monocle.obj <- monocle.obj[, monocle.obj$timepoint != "28d"]
#   monocle.obj <- monocle.obj[, monocle.obj$markers == "lin-"]
#   monocle.obj <- monocle.obj[, monocle.obj$sample != "WT-LSK_OP0_NM_7d_1"]
#   
#   # Filter only celltypes of interest
#   ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
#   ann <- ann[rn %in% colnames(monocle.obj)]
#   sort(unique(ann$Clusters))
#   for(x in c("B.cell", "CLP", "T.cell")){ann <- ann[!grepl(x, Clusters)]}
#   sort(unique(ann$Clusters))
#   
#   # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
#   monocle.obj <- monocle.obj[, ann$rn]
#   monocle.obj$clusterDE <- ann$Clusters
#   
#   out <- dirout(paste0(basedir, tissue.name, "_", ct.use))
#   
#   source("src/SCRNA_30_01_DE_FUNC_nebula.R")
# }


# # Analysis of most/all cells in vivo ---------------------------------------------------------------
# cells <- list(
#   in.vivo = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_in.vivo_noMixscape_Cells.tsv"))$rn,
#   ex.vivo = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_ex.vivo_noMixscape_Cells.tsv"))$rn,
#   leukemia = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_leukemia_noMixscape_Cells.tsv"))$rn
# )
# ct.use <- "everything"
# (tissue.name <- "leukemia")
# #for(tissue.name in "leukemia"){
# for(tissue.name in TISSUES){
#   
#   (timex <- mobjs[[tissue.name]]$timepoint[1])
#   for(timex in unique(mobjs[[tissue.name]]$timepoint)){
#     
#     # Monocle object
#     monocle.obj <- mobjs[[tissue.name]]
#     monocle.obj <- monocle.obj[, monocle.obj$timepoint == timex]
#     monocle.obj <- monocle.obj[, monocle.obj$sample != "WT-LSK_OP0_NM_7d_1"]
#     monocle.obj <- monocle.obj[, colnames(monocle.obj) %in% cells[[tissue.name]]]
#     
#     # Remove Mixscape effects
#     monocle.obj$mixscape_class.global[monocle.obj$mixscape_class.global == "NP"] <- "KO"
#     
#     if(ncol(monocle.obj) == 0) next
#     
#     # Filter only celltypes of interest
#     ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
#     ann <- ann[rn %in% colnames(monocle.obj)]
#     sort(unique(ann$Clusters))
#     # for(x in c("B.cell", "CLP", "T.cell")){ann <- ann[!grepl(x, Clusters)]}
#     # sort(unique(ann$Clusters))
#     
#     # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
#     monocle.obj <- monocle.obj[, ann$rn]
#     monocle.obj$clusterDE <- ann$Clusters
#     
#     out <- dirout(paste0(basedir, tissue.name, "_", timex, "_", ct.use))
#     
#     source("src/SCRNA_30_01_DE_FUNC_nebula.R")
#   }
# }


# # Analysis of most/all cells in vivo ---------------------------------------------------------------
cells <- list(
  in.vivo = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_in.vivo_noMixscape_Cells.tsv"))$rn,
  ex.vivo = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_ex.vivo_noMixscape_Cells.tsv"))$rn,
  leukemia = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_leukemia_noMixscape_Cells.tsv"))$rn
)
ct.use <- "everythingMonocle"
(tissue.name <- "leukemia")
for(tissue.name in TISSUES){

  (timex <- mobjs[[tissue.name]]$timepoint[1])
  for(timex in unique(mobjs[[tissue.name]]$timepoint)){

    # Monocle object
    monocle.obj <- mobjs[[tissue.name]]
    monocle.obj <- monocle.obj[, monocle.obj$timepoint == timex]
    monocle.obj <- monocle.obj[, monocle.obj$sample != "WT-LSK_OP0_NM_7d_1"]
    monocle.obj <- monocle.obj[, colnames(monocle.obj) %in% cells[[tissue.name]]]

    # Remove Mixscape effects
    #monocle.obj$mixscape_class.global[monocle.obj$mixscape_class.global == "NP"] <- "KO"

    if(ncol(monocle.obj) == 0) next

    # Filter only celltypes of interest
    ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
    ann <- ann[rn %in% colnames(monocle.obj)]
    sort(unique(ann$Clusters))
    # for(x in c("B.cell", "CLP", "T.cell")){ann <- ann[!grepl(x, Clusters)]}
    # sort(unique(ann$Clusters))

    # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
    monocle.obj <- monocle.obj[, ann$rn]
    monocle.obj$clusterDE <- ann$Clusters

    out <- dirout(paste0(basedir, tissue.name, "_", timex, "_", ct.use))

    source("src/SCRNA_30_01_DE_FUNC_nebula.R")
  }
}