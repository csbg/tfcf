

# Environmental variables for working on the CAME Cluster ----------------------------------------
if(Sys.getenv("GFS") == "" | Sys.getenv("GFS") == "/usr/local/AGFORTELNY/"){  # GFS is not set (we are on the CAME Cluster)
  print("Setting GFS PATH")
  if(dir.exists("/media/AGFORTELNY")){  # we are on the CAME Cluster
    Sys.setenv(CODEBASE=paste0(Sys.getenv("HOME"), "/code/"))
    Sys.setenv(GFS="/media/AGFORTELNY")
  } else {
    stop("Environmental variables missing")
  }
}

if(Sys.getenv("DATA") == ""){
	print("Setting DATA and others paths")
	Sys.setenv(
	  CODE=paste0(Sys.getenv("CODEBASE"), "/tfcf/"),
	  DATA=paste0(Sys.getenv("GFS"), "/PROJECTS/TfCf/Data/"),
	  RAWDATA=paste0(Sys.getenv("GFS"), "/PROJECTS/TfCf/RawData/"),
	  ANALYSIS=paste0(Sys.getenv("GFS"), "/PROJECTS/TfCf/Analysis/")
	)
}


# # Setwd -------------------------------------------------------------------
# setwd(paste0(Sys.getenv("CODE")))
# 

.libPaths()

# Packages and functions ----------------------------------------------------------------
require(data.table)
require(ggplot2)
require(pheatmap)
# require(renv)
require(hdf5r)
require(dplyr)
require(tidyr)
require(Seurat)
require(monocle3)
require(mixtools)


source(paste(Sys.getenv("CODEBASE"), "resources", "RFunctions", "Basics.R", sep="/"))
source(paste(Sys.getenv("CODEBASE"), "resources", "RFunctions", "scRNA_Basics.R", sep="/"))
source(paste(Sys.getenv("CODEBASE"), "resources", "RFunctions", "GSEA_hitlist.R", sep="/"))


# RENV LOCKFILE -----------------------------------------------------------
if(!dir.exists("lockfiles/")) dir.create("lockfiles/")
# renv::snapshot(lockfile = paste0("lockfiles/renv_lockfile_", make.names(Sys.time()), ".lockfile"), force = TRUE, prompt = FALSE)



# Enrichr DBs -------------------------------------------------------------
ENRICHR.DBS <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse",
                 "MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures",
                 "NCI-60_Cancer_Cell_Lines", "Cancer_Cell_Line_Encyclopedia",
                 "TRANSFAC_and_JASPAR_PWMs", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "TRRUST_Transcription_Factors_2019"
                 )

# FUNCTIONS ---------------------------------------------------------------
dirout <- function(out, ext="", init=TRUE){
  out.dir <- paste0(Sys.getenv("ANALYSIS"), "/", out, "/")
  if(init){
    dir.create(out.dir,showWarnings=FALSE); 
    message("Setting output directory: ", out.dir)
  }
  function(...){
    paste0(out.dir, paste0(...), ext)
  }
}

# The dirout_load function is used load data from another output folder. Does not create a directory
dirout_load <- function(out, ext=""){
  dirout(out=out, ext=ext, init=FALSE)
}

# Get main datasets
getMainDatasets <- function(){
  ff <- list.files(Sys.getenv("DATA"))
  ff <- ff[!grepl(".log$", ff)]
  ff <- ff[!grepl("_onlyRNA", ff)]
  ff <- ff[!grepl("RNAonly", ff)]
  ff <- ff[grepl("^ECCITE", ff) | grepl("^CITESEQ", ff)]
  ff <- ff[!grepl("ECCITE4_INT", ff)]
  ff <- ff[!grepl("ECCITE1_", ff)]
  ff <- ff[!grepl("LINES", ff)]
  #ff <- ff[!grepl("ECCITE8", ff)]
  list(folders=ff, dir=Sys.getenv("DATA"))
}




# Comparisons (POOLED) -------------------------------------------------------------
COMPARISONS <- list(
  LSK.CKIT=c("LSKd9", "cKit"),
  GMP.LSK=c("GMP", "LSKd7"),
  MEP.LSK=c("MEP", "LSKd7"),
  GMP.MEP=c("GMP", "MEP"),
  UND.MEP=c("Und", "MEP"),
  MYE.GMP=c("Mye", "GMP"),
  MYE.UND=c("Mye", "Und"),
  GMPcd11.DN=c("GMP.CD11bGr1", "GMP.DN"),
  DMCD34pos.neg=c("DM.CD34pos", "DM.CD34neg"),
  DMCFSEhigh.low=c("DM.CFSEhigh", "DM.CFSElow"),
  DMEry.LSC=c("DM.Ery", "DM.LSC"),
  DMMye.LSC=c("DM.Mye", "DM.LSC"),
  DMMye.Ery=c("DM.Mye", "DM.Ery"),
  DMCD11b.LSC=c("DM.LSC.CD11b", "DM.LSC")
)

COMPARISONS.healthy <- c(
  "LSK.CKIT",
  "GMP.LSK",
  "MEP.LSK",
  "GMP.MEP",
  "MYE.UND",
  "GMPcd11.DN"
)
cleanComparisons <- function(x, order=TRUE, ggtext=FALSE, dm="clean"){
  transformPretty <- function(i, gg=ggtext){
    if(gg){
      i <- gsub("^(.+)\\.(.+)$", "<strong style='color:#832424;'>\\1</strong> vs <strong style='color:#3A3A98;'>\\2</strong>", i)
    } else {
      i <- gsub("\\.", " vs ", i)
    }
    if(dm == "rm"){
      i <- sub("^DM", "", i)
    }
    if(dm == "clean"){
      i <- sub("^DM", "DM: ", i)
    }
    return(i)
  }
  
  x <- transformPretty(x)
  if(order) x <- factor(x, levels = intersect(transformPretty(names(COMPARISONS)), unique(x)))
  x
}


# Paths -------------------------------------------------------------------

PATHS <- list()

PATHS$RESOURCES <- list(
  HM.MAP = dirout_load("PPI_00_getData/")("BioMart_Human_Mouse_2021_07_27.txt"),
  Enrichr.mouse = dirout_load("EXT_02_EnrichR_Genesets")("Genesets_Mouse.RData")
)

PATHS$CHIP <- list()
PATHS$CHIP$Targets <- dirout_load("CHIP_20_01_Peaks_julen")("ChIP.Targets.RData")

PATHS$POOLED <- list()
PATHS$POOLED$DATA <- list(
  matrix=dirout_load("POOLED_01_CollectData")("Matrix.csv"),
  annotation=dirout_load("POOLED_01_CollectData")("Annotation.tsv"),
  matrix.aggregated=dirout_load("POOLED_09_CleanData")("Matrix_aggregated.csv"),
  annotation.aggregated=dirout_load("POOLED_09_CleanData")("Annotation_aggregated.tsv")
)
sapply(PATHS$POOLED$DATA, file.exists)


PATHS$FULLINT <- list()
PATHS$FULLINT$Monocle <- dirout_load("FULLINT_01_01_Integration")("MonocleObject.RData")
PATHS$FULLINT$DEG <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")("DEG_Results_nebula.RData")
PATHS$FULLINT$DEG.clean <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")("DEG_Results_all.tsv")
PATHS$FULLINT$DEG.clean.leukemia <- dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")("DEG_Results_all.tsv")
PATHS$FULLINT$DEG.clean.invitro <- dirout_load("FULLINT_10_01_BasicAnalysis_in.vitro")("DEG_Results_all.tsv")
PATHS$FULLINT$DEG.clean.invivo <- dirout_load("FULLINT_10_01_BasicAnalysis_in.vivo")("DEG_Results_all.tsv")
PATHS$FULLINT$DEG.ann <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")("DEG_Annnotation.tsv")
PATHS$FULLINT$DEG.logFCMT <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")("DEG_Results_logFCMT.csv")
PATHS$FULLINT$DEG.UMAP <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")("RegulatoryMap_UMAP_top.genes.tsv")


PATHS$CITESEQ1 <- list()
PATHS$CITESEQ1$DATA <- list(
  matrix=paste(Sys.getenv("DATA"), "CITESEQ1", "outs", "filtered_feature_bc_matrix.h5", sep="/"),
  umap=paste(Sys.getenv("DATA"), "CITESEQ1", "outs", "analysis", "umap", "2_components", "projection.csv", sep="/"),
  clusters=paste(Sys.getenv("DATA"), "CITESEQ1", "outs", "analysis", "clustering", "graphclust", "clusters.csv", sep="/"),
  de=paste(Sys.getenv("DATA"), "CITESEQ1", "outs", "analysis", "diffexp", "graphclust", "differential_expression.csv", sep="/")
)
sapply(PATHS$CITESEQ1$DATA, file.exists)

PATHS$CITESEQ1_CLEAN <- list()
PATHS$CITESEQ1_CLEAN$DATA <- list(
  matrix=paste(Sys.getenv("DATA"), "CITESEQ1_CLEAN", "outs", "filtered_feature_bc_matrix.h5", sep="/"),
  umap=paste(Sys.getenv("DATA"), "CITESEQ1_CLEAN", "outs", "analysis", "umap", "2_components", "projection.csv", sep="/"),
  clusters=paste(Sys.getenv("DATA"), "CITESEQ1_CLEAN", "outs", "analysis", "clustering", "graphclust", "clusters.csv", sep="/"),
  de=paste(Sys.getenv("DATA"), "CITESEQ1_CLEAN", "outs", "analysis", "diffexp", "graphclust", "differential_expression.csv", sep="/")
)
sapply(PATHS$CITESEQ1_CLEAN$DATA, file.exists)

PATHS$CITESEQ2$DATA <- list(
  matrix=paste(Sys.getenv("DATA"), "CITESEQ2", "outs", "filtered_feature_bc_matrix.h5", sep="/"),
  umap=paste(Sys.getenv("DATA"), "CITESEQ2", "outs", "analysis", "umap", "2_components", "projection.csv", sep="/"),
  clusters=paste(Sys.getenv("DATA"), "CITESEQ2", "outs", "analysis", "clustering", "graphclust", "clusters.csv", sep="/"),
  de=paste(Sys.getenv("DATA"), "CITESEQ2", "outs", "analysis", "diffexp", "graphclust", "differential_expression.csv", sep="/")
)
sapply(PATHS$CITESEQ2$DATA, file.exists)

PATHS$CITESEQ2_CLEAN <- list()
PATHS$CITESEQ2_CLEAN$DATA <- list(
  matrix=paste(Sys.getenv("DATA"), "CITESEQ2_CLEAN", "outs", "filtered_feature_bc_matrix.h5", sep="/"),
  umap=paste(Sys.getenv("DATA"), "CITESEQ2_CLEAN", "outs", "analysis", "umap", "2_components", "projection.csv", sep="/"),
  clusters=paste(Sys.getenv("DATA"), "CITESEQ2_CLEAN", "outs", "analysis", "clustering", "graphclust", "clusters.csv", sep="/"),
  de=paste(Sys.getenv("DATA"), "CITESEQ2_CLEAN", "outs", "analysis", "diffexp", "graphclust", "differential_expression.csv", sep="/")
)
sapply(PATHS$CITESEQ2_CLEAN$DATA, file.exists)

# PATHS$ECCITE1 <- list()
# PATHS$ECCITE1$DATA <- list(
#   matrix=paste(Sys.getenv("DATA"), "ECCITE1_RNA_cellranger", "outs", "filtered_feature_bc_matrix.h5", sep="/"),
#   umap=paste(Sys.getenv("DATA"), "ECCITE1_RNA_cellranger", "outs", "analysis", "umap", "2_components", "projection.csv", sep="/"),
#   clusters=paste(Sys.getenv("DATA"), "ECCITE1_RNA_cellranger", "outs", "analysis", "clustering", "graphclust", "clusters.csv", sep="/"),
#   de=paste(Sys.getenv("DATA"), "ECCITE1_RNA_cellranger", "outs", "analysis", "diffexp", "graphclust", "differential_expression.csv", sep="/"),
#   guides=paste(Sys.getenv("ANALYSIS"), "ECCITESEQ1_01_Guides_510", "citeseq_guides2cells2.csv", sep="/")
# )
# sapply(PATHS$ECCITE1$DATA, file.exists)

PATHS$ECCITE1 <- list()
PATHS$ECCITE1$DATA <- list(
  matrix=paste(Sys.getenv("DATA"), "ECCITE1_RNA_cellranger_601", "outs", "filtered_feature_bc_matrix.h5", sep="/"),
  umap=paste(Sys.getenv("DATA"), "ECCITE1_RNA_cellranger_601", "outs", "analysis", "umap", "2_components", "projection.csv", sep="/"),
  clusters=paste(Sys.getenv("DATA"), "ECCITE1_RNA_cellranger_601", "outs", "analysis", "clustering", "graphclust", "clusters.csv", sep="/"),
  de=paste(Sys.getenv("DATA"), "ECCITE1_RNA_cellranger_601", "outs", "analysis", "diffexp", "graphclust", "differential_expression.csv", sep="/"),
  guides=paste(Sys.getenv("ANALYSIS"), "ECCITE1_01_Guides", "citeseq_guides2cells.csv", sep="/")
)
sapply(PATHS$ECCITE1$DATA, file.exists)

PATHS$ECCITE2 <- list()
PATHS$ECCITE2$DATA <- list(
  matrix=paste(Sys.getenv("DATA"), "ECCITE2", "outs", "filtered_feature_bc_matrix.h5", sep="/"),
  umap=paste(Sys.getenv("DATA"), "ECCITE2", "outs", "analysis", "umap", "2_components", "projection.csv", sep="/"),
  clusters=paste(Sys.getenv("DATA"), "ECCITE2", "outs", "analysis", "clustering", "graphclust", "clusters.csv", sep="/"),
  de=paste(Sys.getenv("DATA"), "ECCITE2", "outs", "analysis", "diffexp", "graphclust", "differential_expression.csv", sep="/")
)
sapply(PATHS$ECCITE2$DATA, file.exists)


# COLORS ------------------------------------------------------------------
COLOR.Genotypes = c(WT="#33a02c", Cas9="#6a3d9a")
COLORS.HM.FUNC <- colorRampPalette(c("#6a3d9a", "#a6cee3", "white", "#fdbf6f", "#e31a1c"))

scale_fill_hexbin <- function(...){scale_fill_gradientn(colours=c("#a6cee3", "#fdbf6f", "#ff7f00", "#e31a1c"), ...)}
