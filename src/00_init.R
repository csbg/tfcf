# packages ----------------------------------------------------------------
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


# Environmental variables ----------------------------------------
# GFS and CODEBASE are required

# CODEBASE is where the code is (this code and also other functions)
if(Sys.getenv("CODEBASE") == ""){
  print("Setting CODEBASE")
  Sys.setenv(CODEBASE=paste0(Sys.getenv("HOME"), "/code/"))
}

# GFS points to where the data is stored and results will be written to
if(dir.exists("/media/AGFORTELNY")){
  print("Setting GFS within singularity to /media/AGFORTELNY")
  Sys.setenv(GFS="/media/AGFORTELNY")
}



# SET PATHS -------------------------------------------------------------
# CODE, ANALYSIS, DATA, RAWDATA

PATHS <- list()
PATHS$LOCATIONS <- list()


# Get paths from setup.sh
ll <- readLines("setup.sh")
ll <- grep("^\\s*export ", ll, value=TRUE)
ll <- grep("^\\s*export PATH", ll, value=TRUE, invert=TRUE)
ll <- data.table(
  var = gsub("^.*export (.+?) ?\\=.+$", "\\1", ll),
  path = sapply(strsplit(gsub("^.+?=(.+)$", "\\1", ll), "/"), function(pathx){
    paste(
      sapply(gsub('"', "", pathx), function(x){
        if(grepl("^\\$",x)) Sys.getenv(gsub("\\$", "", x)) else x
        }), 
      collapse="/")
  })
)


# define paths (from setup.sh or environmental variable if set)
for(varx in ll$var){
  pathx <- if(Sys.getenv(varx) != "" & dir.exists(Sys.getenv(varx))) Sys.getenv(ll[var == varx]$var) else ll[var == varx]$path
  PATHS$LOCATIONS[[varx]] <- pathx
}
stopifnot(all(sapply(PATHS$LOCATIONS, dir.exists)))


# functions ----------------------------------------------------------------
.libPaths()
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
  out.dir <- paste0(PATHS$LOCATIONS$ANALYSIS, "/", out, "/")
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
  ff <- list.files(PATHS$LOCATIONS$DATA)
  ff <- ff[!grepl(".log$", ff)]
  ff <- ff[!grepl("_onlyRNA", ff)]
  ff <- ff[!grepl("RNAonly", ff)]
  #ff <- ff[grepl("^ECCITE", ff) | grepl("^CITESEQ", ff)]
  ff <- ff[!grepl("ECCITE4_INT", ff)]
  ff <- ff[!grepl("ECCITE1_", ff)]
  ff <- ff[!grepl("LINES", ff)]
  #ff <- ff[!grepl("ECCITE8", ff)]
  list(folders=ff, dir=PATHS$LOCATIONS$DATA)
}




# COMPARISONS, NAMES, DEFINITIONS -------------------------------------------------------------
COMPARISONS <- list(
  CKIT.LSK=c("cKit", "LSKd9"),
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
  DMCD11b.LSC=c("DM.LSC.CD11b", "DM.LSC"),
  DMd10v5=c("DM.d10", "DM.d5"),
  DMd17v10=c("DM.d17", "DM.d10")
)

COMPARISONS.healthy <- c(
  "CKIT.LSK",
  "GMP.LSK",
  "MEP.LSK",
  "GMP.MEP",
  "MYE.UND",
  "GMPcd11.DN"
)

CLEAN.CELLTYPES <- fread("metadata/FIGS_celltypes.tsv", fill=TRUE)
CLEAN.CELLTYPES[NewName == "", NewName := Name]

# CLEAN FUNCTIONS ---------------------------------------------------------
cleanComparisons <- function(x, order=TRUE, ggtext=FALSE, dm="clean", reverse=FALSE, colors=c("832424", "3A3A98")){
  transformPretty <- function(i, gg=ggtext){
    if(gg){
      i <- gsub("^(.+)\\.(.+)$", paste0("<strong style='color:#",colors[1],";'>\\1</strong> vs <strong style='color:#",colors[2],";'>\\2</strong>"), i)
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
  if(order){
    x <- factor(x, levels = c("Main branch", intersect(transformPretty(names(COMPARISONS)), unique(x))))
    if(reverse){
      x <- factor(x, levels=rev(levels(x)))
    }
  }
  x
}


cleanCelltypes <- function(x, order=TRUE, twoLines=FALSE, reverse=TRUE, clean=TRUE){
  if(clean){
    stopifnot(all(x %in% CLEAN.CELLTYPES$Name))
    x <- CLEAN.CELLTYPES[match(x, Name)]$NewName
  }
  order.levels <- if(reverse) rev(CLEAN.CELLTYPES$NewName) else CLEAN.CELLTYPES$NewName
  if(order) x <- factor(x, levels=order.levels)
  if(twoLines) x <- sub(" ", "\n", x)
  return(x)
}


# Paths -------------------------------------------------------------------

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


PATHS$SCRNA <- list()
PATHS$SCRNA$ANN <- dirout_load("SCRNA_01_01_Seurat")("SampleAnnotation.tsv")
PATHS$SCRNA$MONOCLE.NAMES <- setdiff(list.dirs(dirout_load("SCRNA_02_01_Integration")(""), full.names = FALSE), "")
PATHS$SCRNA$MONOCLE.DIR <- function(x){dirout_load("SCRNA_02_01_Integration")(paste0(x, "/", "MonocleObject.RData"))}
PATHS$SCRNA$Citeseq <- dirout_load("SCRNA_02_01_Integration")("CITESEQ_Antibodies.RData")
# PATHS$FULLINT$Monocle <- dirout_load("FULLINT_01_01_Integration")("MonocleObject.RData")
# PATHS$FULLINT$Citeseq <- dirout_load("FULLINT_01_01_Integration")("CITESEQ_Antibodies.RData")
# PATHS$FULLINT$DEG <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")("DEG_Results_nebula.RData")
# PATHS$FULLINT$DEG.clean <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")("DEG_Results_all.tsv")
# PATHS$FULLINT$DEG.clean.leukemia <- dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")("DEG_Results_all.tsv")
# PATHS$FULLINT$DEG.clean.invitro <- dirout_load("FULLINT_10_01_BasicAnalysis_in.vitro")("DEG_Results_all.tsv")
# PATHS$FULLINT$DEG.clean.invivo <- dirout_load("FULLINT_10_01_BasicAnalysis_in.vivo")("DEG_Results_all.tsv")
# PATHS$FULLINT$DEG.ann <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")("DEG_Annnotation.tsv")
# PATHS$FULLINT$DEG.logFCMT <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")("DEG_Results_logFCMT.csv")
# PATHS$FULLINT$DEG.UMAP <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")("RegulatoryMap_UMAP_top.genes.tsv")


# COLORS ------------------------------------------------------------------
COLOR.Genotypes = c(WT="#33a02c", Cas9="#6a3d9a")
COLORS.HM.FUNC <- colorRampPalette(c("#6a3d9a", "#a6cee3", "white", "#fdbf6f", "#e31a1c"))
COLORS.CELLTYPES.scRNA <- fread("metadata/markers.signatures.use.scRNA.tsv")
COLORS.CELLTYPES.scRNA <- setNames(COLORS.CELLTYPES.scRNA$Color, COLORS.CELLTYPES.scRNA$FinalName)
COLORS.CELLTYPES.scRNA["NA"] <- "lightgrey"
COLORS.CELLTYPES.scRNA.ainhoa <- setNames(CLEAN.CELLTYPES$Color, CLEAN.CELLTYPES$NewName)

scale_fill_hexbin <- function(...){scale_fill_gradientn(colours=c("#a6cee3", "#fdbf6f", "#ff7f00", "#e31a1c"), ...)}



message("-------------------> Project initiation completed")