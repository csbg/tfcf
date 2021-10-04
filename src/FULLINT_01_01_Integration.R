source("src/00_init.R")
out <- dirout("FULLINT_01_01_Integration/")


# Read cellranger analysis results --------------------------------------------
AGG.CSV <- fread(paste(Sys.getenv("DATA"), "FULLINT_00_Aggr", "outs", "aggregation.csv", sep="/"))
AGG.CSV$i <- 1:nrow(AGG.CSV)

ff <- list.files(Sys.getenv("DATA"))
ff <- ff[!grepl(".log$", ff)]
ff <- ff[!grepl("_onlyRNA", ff)]
ff <- ff[!grepl("RNAonly", ff)]
ff <- ff[grepl("ECCITE", ff) | grepl("CITESEQ", ff)]
ff <- ff[!grepl("ECCITE4_INT", ff)]
ff <- ff[!grepl("ECCITE1_", ff)]

marker.genes <- fread("metadata/markers.csv")


# Read data, Seurat processing, and Monocle integration --------------------------------------------
monocle.file <- out("MonocleObject.RData")
if(file.exists(monocle.file)){
  print("Loading Monocle file")
  load(monocle.file)

} else {
  
  print("Processing datasets")
  
  # Prepare lists to store data
  monocle.obj.list <- list()
  additional.info <- list()

  
  dsx <- ff[1]
  dsx <- "ECCITE1"
  for(dsx in ff){
    
    # already processed?
    if(dsx %in% names(monocle.obj.list)){
      message(dsx, " already processed")
      next
    }
    
    print("------------------------------")
    print(dsx)
    print("---------------")
    path <- paste(Sys.getenv("DATA"), dsx, "outs", "filtered_feature_bc_matrix.h5", sep = "/")
    
    dsx.file <- out(paste0("SeuratObj_",dsx,".RData"))
    
    # Load or process file
    if(file.exists(dsx.file)){
      print("Loading processed dataset")
      load(dsx.file)
    } else {
      print("Processing dataset")
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
      
      if(dsx == "ECCITE1"){
        read.csv(PATHS$ECCITE1$DATA$guides)
        matrix_dir = paste(Sys.getenv("DATA"), "ECCITE1_citeseq_combined//umi_count/", sep="/")
        barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
        features.path <- paste0(matrix_dir, "features.tsv.gz")
        matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
        mat <- readMM(file = matrix.path)
        feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
        barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
        colnames(mat) = paste0(barcode.names$V1, "-1")
        rownames(mat) = gsub("\\-.+$", "", feature.names$V1)
        mat <- mat[row.names(mat) != "unmapped",]
        additional.info[[dsx]] <- list("CRISPR Guide Capture" = as(mat, "dgCMatrix"))
      }
      
      if(class(data.gx) != "dgCMatrix"){
        message("Data for ", dsx, " not in the right format")
      }
    
      seurat.obj <- CreateSeuratObject(counts = data.gx, project = dsx)
      seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
      seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
      seurat.obj <- NormalizeData(seurat.obj, verbose = FALSE)
      seurat.obj <- CellCycleScoring(seurat.obj, s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes,set.ident = TRUE)
      
      # Add additional info to metadata
      tx <- names(additional.info[[dsx]])[1]
      for(tx in names(additional.info[[dsx]])){
        txn <- make.names(tx)
        print(tx)
        
        # Get dataset
        x <- additional.info[[dsx]][[tx]]
        if(class(x) != "dgCMatrix") next
        
        # Get features for each cells
        if(txn == "CRISPR.Guide.Capture"){
          colcnts <- colSums(x)
          colcnts.min <- min(5, median(colcnts[colcnts != 0]))
          x.max <- apply(x, 2, max)
          x.sum <- apply(x, 2, sum)
          x <- as.matrix(x[,x.max/x.sum > 0.75 & x.sum >= colcnts.min])
          ii <- apply(x, 2, function(col) which(col == max(col)))
        } else {
          cutoff <- 5
          x <- as.matrix(x[,apply(x, 2, max) >= cutoff,drop=F]) # only keep cells with any guide counted above cutoff
          ii <- apply(x, 2, function(col) which(col >= cutoff)) # for each cell (column) get the rows that are above the cutoff
        }
        
        # label features and combine them
        if(length(ii) == 0) next
        labelsx <- sapply(ii, function(ix) paste(row.names(x)[ix], collapse = ","))
        seurat.obj@meta.data[[txn]] <- labelsx[row.names(seurat.obj@meta.data)]
        print(table(seurat.obj@meta.data[[txn]]))
      }

      # Mixscape
      if("CRISPR.Guide.Capture" %in% colnames(seurat.obj@meta.data)){
        orig.guides <- seurat.obj@meta.data$CRISPR.Guide.Capture
        guides.clean <- orig.guides
        guides.clean[grepl(",", orig.guides)] <- NA
        guides.clean[grepl("^NTC_", orig.guides)] <- "NTC"
        guides.clean[guides.clean %in% names(which(table(guides.clean) < 5))] <- NA
        
        if(sum(!is.na(guides.clean)) > 5){
          seurat.obj@meta.data$guide <- guides.clean
  
          # Mixscape
          eccite <- subset(seurat.obj, cells=row.names(seurat.obj@meta.data)[!is.na(seurat.obj$guide)])
          eccite <- subset(eccite, features=names(which(Matrix::rowSums(eccite@assays$RNA@counts) > 10)))
          eccite <- FindVariableFeatures(eccite, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
          #eccite <- subset(eccite, features=eccite@assays$RNA@var.features)
          eccite <- ScaleData(eccite, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
          eccite <- RunPCA(eccite, npcs = 30, verbose = FALSE)
          
          # CalcPerturbSig
          for(nn in seq(from = 20,to = 5, by = -5)){
            tryCatch({
              eccite <- CalcPerturbSig(
                object = eccite, assay = "RNA", slot = "data",
                gd.class ="guide",nt.cell.class = "NTC",
                reduction = "pca",ndims = 20, num.neighbors = nn,
                new.assay.name = "PRTB")
              last
            }, error=function(e){
              message(e)
            })
          }
          
          # RunMixscape
          eccite <- ScaleData(object = eccite, assay = "PRTB", do.scale = F, do.center = T)
          eccite <- RunMixscape(
              object = eccite,assay = "PRTB",slot = "scale.data",
              labels = "guide",nt.class.name = "NTC",
              min.de.genes = 5,
              iter.num = 10,
              de.assay = "RNA",
              verbose = TRUE,
              prtb.type = "KO")
          for(cx in c("mixscape_class", "mixscape_class_p_ko", "mixscape_class.global")){
            seurat.obj@meta.data[[cx]] <- eccite@meta.data[row.names(seurat.obj@meta.data),cx]
          }
        }
      }

      seurat.obj <- SCRNA.DietSeurat(sobj = seurat.obj)
      save(seurat.obj, file=dsx.file)
    }
    
    # Extract information for Monocle
    mat.use <- seurat.obj@assays$RNA@counts
    stopifnot(!any(duplicated(row.names(mat.use))))
    if(!"GFP" %in% row.names(mat.use)){
      x <- matrix(0, nrow = 2, ncol = ncol(mat.use))
      row.names(x) <- c("GFP", "BFP")
      mat.use <- rbind(mat.use, x)
    }
    monocle.obj.list[[dsx]] <- new_cell_data_set(expression_data = mat.use, cell_metadata = seurat.obj@meta.data)
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
  monocle.obj <- 
    align_cds(monocle.obj, alignment_group = "sample", verbose = TRUE) %>% 
    reduce_dimension(
      reduction_method = "UMAP",
      preprocess_method = "Aligned",
      verbose = TRUE)
  
  # Store full dataset
  save(monocle.obj, additional.info, file=monocle.file)
}

ann <- data.table(data.frame(colData(monocle.obj)@listData), keep.rownames = TRUE)
ann[,.N, by=c("guide", "mixscape_class.global", "sample")][!is.na(guide)][mixscape_class.global == "KO"]


x <- data.table(reducedDims(monocle.obj)$UMAP, keep.rownames = TRUE)
x[,sample := gsub(".+?\\-1_", "", rn)]

ggplot(x, aes(x=V1, y=V2)) + geom_hex()
ggplot(x, aes(x=V1, y=V2)) + geom_hex() + facet_wrap(~sample)
