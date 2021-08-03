setwd(paste0(Sys.getenv("CODE")))

source("~/code/resources/RFunctions/Basics.R")
source("~/code/resources/RFunctions/scRNA_Basics.R")


# Packages ----------------------------------------------------------------
require(data.table)
require(ggplot2)
require(pheatmap)


ENRICHR.DBS <- c("KEGG_2019_Mouse", "NCI-Nature_2016", "WikiPathways_2019_Mouse", "Reactome_2016",
                 "TRANSFAC_and_JASPAR_PWMs", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "ENCODE_TF_ChIP-seq_2015", "ChEA_2016", "TRRUST_Transcription_Factors_2019"
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



# Comparisons -------------------------------------------------------------
COMPARISONS <- list(
  CKIT.LSK=c("cKit", "LSKd9"),
  GMP.LSK=c("GMP", "LSKd7"),
  MEP.LSK=c("MEP", "LSKd7"),
  GMP.MEP=c("GMP", "MEP"),
  UND.MEP=c("Und", "MEP"),
  MYE.GMP=c("Mye", "GMP"),
  UND.MYE=c("Und", "Mye"),
  GMPcd11.DN=c("GMP.CD11bGr1", "GMP.DN"),
  CD34pos.neg=c("CD34pos", "CD34neg"),
  CFSEhigh.low=c("CFSEhigh", "CFSElow")
)

COMPARISONS.USE <- c(
  "CKIT.LSK",
  "GMP.LSK",
  "MEP.LSK",
  "GMP.MEP",
  "UND.MYE",
  "GMPcd11.DN"
)
cleanComparisons <- function(x, order=TRUE){
  x <- gsub("\\.", " vs ", x)
  if(order) x <- factor(x, levels = intersect(gsub("\\.", " vs ", names(COMPARISONS)), unique(x)))
  x
}


# Paths -------------------------------------------------------------------

PATHS <- list()

PATHS$POOLED <- list()
PATHS$POOLED$DATA <- list(
  matrix=dirout_load("POOLED_01_CollectData")("Matrix.csv"),
  annotation=dirout_load("POOLED_01_CollectData")("Annotation.tsv")
)
sapply(PATHS$POOLED$DATA, file.exists)




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

