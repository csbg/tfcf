source("src/00_init.R")
basedir <- "SCRNA_02_01_Integration/"
out.base <- dirout(basedir)
inDir <- dirout_load("SCRNA_01_01_Seurat")

require("sceasy")

# Annotation --------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)
samples.not.use <- SANN[grepl("\\d$", tissue) | toAnalyse == "NO"]$sample
SANN <- SANN[!sample %in% samples.not.use]
sort(SANN$sample)

# Integrate data with Monocle
tx <- SANN$tissue[1]
for(tx in unique(SANN$tissue)){
  out <- dirout(paste0(basedir, tx, "/"))
  monocle.file <- out("MonocleObject.RData")
  
  monocle.obj.list <- list()
  
  sx <- SANN[tissue == tx]$sample[1]
  for(sx in SANN[tissue == tx]$sample){
    fx <- inDir("SeuratObj_", sx, ".RData")
    if(!file.exists(fx)) stop(fx, " seurat object not found")
    
    print(paste("Reading ",sx))
    load(fx)
    
    # Add more annotation
    for(x in c("tissue", "markers", "timepoint", "sample", "sample_broad")){
      seurat.obj@meta.data[[x]] <- SANN[sample == sx][[x]]
    }
    
    # Extract information for Monocle
    mat.use <- seurat.obj@assays$RNA@counts
    stopifnot(!any(duplicated(row.names(mat.use))))
    if(!"GFP" %in% row.names(mat.use)){
      x <- matrix(0, nrow = 2, ncol = ncol(mat.use))
      row.names(x) <- c("GFP", "BFP")
      mat.use <- rbind(mat.use, x)
    }
    monocle.obj.list[[sx]] <- new_cell_data_set(expression_data = mat.use, cell_metadata = seurat.obj@meta.data)
    
    # Store CITE-seq data
    if(sx == "DM_CITEseq-2_NA_NM_1"){
      citeseq.MT <- additional.info.x
      save(citeseq.MT, file=out.base("CITESEQ_Antibodies.RData"))
    }
  }
  
  # Make sure all objects have the same number of rows
  stopifnot(length(unique(sapply(monocle.obj.list, nrow))) == 1)
  
  # Combine datasets
  monocle.obj <- combine_cds(cds_list = monocle.obj.list, cell_names_unique = FALSE)
  
  # Process dataset
  monocle.obj <-
    preprocess_cds(monocle.obj, verbose = TRUE) %>%
    reduce_dimension(preprocess_method = "PCA", verbose = TRUE)
  set.seed(42)
  monocle.obj <- align_cds(monocle.obj, 
      alignment_group = "sample_broad", 
      residual_model_formula_str = "~Phase",
      verbose = TRUE
    )
  monocle.obj <- reduce_dimension(monocle.obj,
      reduction_method = "UMAP",
      preprocess_method = "Aligned",
      verbose = TRUE)
  
  # Clustering
  set.seed(12121)
  monocle.obj = cluster_cells(monocle.obj)

  # Store full dataset
  save(monocle.obj, file=monocle.file)
}



