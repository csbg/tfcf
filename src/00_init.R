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




# Paths -------------------------------------------------------------------

PATHS <- list()

PATHS$DATA <- list(
  matrix=dirout_load("POOLED_01_CollectData")("Matrix.csv"),
  annotation=dirout_load("POOLED_01_CollectData")("Annotation.tsv")
)

PATHS$ANALYSIS <- list(

)

sapply(PATHS$DATA, file.exists)
