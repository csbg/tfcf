source("src/00_init.R")
basedir <- "SCRNA_05_01_SingleR/"
out <- dirout(basedir)

library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)
library(CytoTRACE)
library(doMC)
library(GEOquery)



# Human/Mouse gene mapping ------------------------------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
SANN <- fread(PATHS$SCRNA$ANN)


# Download Izzo et al data ------------------------------------------------
izzo.file <- out("izzo.RData")
if(file.exists(izzo.file)){
  (load(izzo.file))
} else {
  geo_accs <- c("GSM4172942", "GSM3554952", "GSM3554951", "GSM3554950")
  for(geo_acc in geo_accs){
    getGEOSuppFiles(geo_acc, fetch_files=TRUE, baseDir = out(""))
  }
  ff <- list.files(out(""), pattern="*Chromium10X_WT._UMItab.csv.gz", recursive = TRUE, full.names = TRUE)
  fx <- ff[1]
  mats <- list()
  for(fx in ff){
    x <- fread(fx)
    stopifnot(all(sapply(x[,-"V1"], is.integer)))
    g <- x$V1
    x <- as.matrix(x[,-"V1"])
    row.names(x) <- g
    x <- as(x,"dgCMatrix")
    mats[[basename(fx)]] <- x
  }
  names(mats) <- gsub("^.+Chromium10X_(WT\\d)_.+$", "\\1", names(mats))
  izzo.ann <- fread("metadata/IzzoEtAl.metadata.csv")
  izzo.ann <- izzo.ann[grepl("^WT\\d$", orig.ident),]
  izzo.ann[, V1 := gsub("WT6", "WT4", V1)]
  g <- unique(do.call(c, lapply(mats, row.names)))
  mats <- lapply(names(mats), function(mnam){
    m <- mats[[mnam]]
    colnames(m) <- paste0(mnam, "_", colnames(m))
    m <- SCRNA.addGenesToSparseMT(m, setdiff(g, row.names(m)))
    m
  })
  izzoMT <- do.call(cbind, mats)
  stopifnot(length(setdiff(colnames(izzoMT), izzo.ann$V1))==0)
  izzoMT <- izzoMT[,izzo.ann$V1]
  stopifnot(all(colnames(izzoMT) == izzo.ann$V1))
  
  save(izzoMT, izzo.ann, file=izzo.file)
}


# Reference data ----------------------------------------------------------
reference_cell_types <- list(
  # Izzo et al
  izzo = SummarizedExperiment(
    assays = SimpleList(logcounts=SCRNA.TPXToLog(SCRNA.RawToTPX(izzoMT, 1e6))), 
    colData=DataFrame(
      label.main=gsub("\\-?\\d$", "", izzo.ann$clusterName), 
      label.fine=izzo.ann$clusterName,
      row.names = izzo.ann$V1
    )
  ),
  
  # two general purpose datasets
  # hpca = HumanPrimaryCellAtlasData(),
  blueprint = BlueprintEncodeData(),
  
  # # comprehensive CD4+ subsets; only one B cell subset, no dendritic cells
  # dice = DatabaseImmuneCellExpressionData(),
  
  # for bone marrow samples
  dmap = NovershternHematopoieticData(),
  
  # for PBMC
  monaco = MonacoImmuneData(),
  
  # ImmGen
  # immgen=ImmGenData(),
  
  # DB ImmuneCells
  # immuneCellExDB=DatabaseImmuneCellExpressionData(),
  
  # Mouse RNA-seq
  mouseRNA=MouseRNAseqData(),
  
  # CytoTRACE 10x Marrow
  marrow10x=SummarizedExperiment(
    assays = SimpleList(logcounts=SCRNA.TPXToLog(SCRNA.RawToTPX(SCRNA.ToSparseMT(as.matrix(marrow_10x_expr)), 1e6))), 
    colData=DataFrame(
      label.main=marrow_10x_pheno, 
      row.names = names(marrow_10x_pheno)
    )
  ),
  
  # CytoTRACE plate
  marrowPlate=SummarizedExperiment(
    assays = SimpleList(logcounts=SCRNA.TPXToLog(SCRNA.RawToTPX(SCRNA.ToSparseMT(as.matrix(marrow_plate_expr)), 1e6))), 
    colData=DataFrame(
      label.main=marrow_plate_pheno, 
      row.names = names(marrow_plate_pheno)
    )
  )
)

stopifnot(all(colnames(marrow_10x_expr) == names(marrow_10x_pheno)))
stopifnot(all(colnames(marrow_plate_expr) == names(marrow_plate_pheno)))
# head(row.names(reference_cell_types$dmap), 40)
# head(row.names(reference_cell_types$monaco), 40)
# head(row.names(reference_cell_types$blueprint), 40)
# head(row.names(reference_cell_types$hpca), 40)


# Plot expression of *our* genes ------------------------------------------
genes <- fread("metadata/TFCF_Annotations.tsv")
table(genes$Category)
genes <- genes[Category == "CF"]
genes.hm <- hm.map[Gene.name %in% genes$GENE]

res <- data.table()
ref <- names(reference_cell_types)[1]
for(ref in names(reference_cell_types)){
  olh <- sum(row.names(reference_cell_types[[ref]]) %in% genes.hm$Human.gene.name)
  olm <- sum(row.names(reference_cell_types[[ref]]) %in% genes.hm$Gene.name)
  orgx <- if(olh > olm) "human" else "mouse"
  
  m <- assay(reference_cell_types[[ref]], "logcounts")
  gg <- if(orgx == "mouse") genes.hm$Gene.name else genes.hm$Human.gene.name
  gg <- unique(gg[!is.na(gg)])
  m <- m[gg[gg %in% row.names(m)],]
  labelx <- "label.fine"
  for(labelx in c("label.main", "label.fine")){
    if(!labelx %in% colnames(colData(reference_cell_types[[ref]]))) next
    annx <- colData(reference_cell_types[[ref]])[,labelx]
    gx <- annx[1]
    for(gx in unique(annx)){
      mx <- m[, annx == gx, drop=FALSE]
      res <- rbind(res, data.table(
        dataset=ref,
        label=labelx,
        celltype=gx,
        gene=row.names(m),
        logcnt=rowMeans(mx),
        N=ncol(mx)
      ))
    }
  }
}
save(res, file=out("CFs_Datasets.RData"))

# generate plots
dx <- res$dataset[1]
for(dx in unique(res$dataset)){
  lx <- res$label[1]
  for(lx in unique(res[dataset==dx]$label)){
    pDT <- res[dataset==dx & label==lx]
    pMT <- toMT(pDT, row = "gene", col = "celltype", val = "logcnt")
    pMT <- pMT[rowSums(pMT != 0) >= 1,]
    pMT <- t(scale(t(pMT)))
    quantile(pMT)
    cleanDev(); pdf(out("CF_Expression_HM_", dx, "_", lx, ".pdf"),w=ncol(pMT) * 0.2+2,h=60)
    pheatmap(pMT,
             main = paste(dx, lx),
             scale="row",
             fontsize_row = 6,
             breaks=seq(-5,5, 0.05),
             color=colorRampPalette(c("#6a3d9a", "#a6cee3", "white", "#fdbf6f", "#e31a1c"))(201)
             )
    dev.off()
  }
}




# LOAD DATA ---------------------------------------------------------------

sx <- SANN$sample[1]
for(sx in SANN$sample){
  
  # Load matrix
  fx <- inDir("SeuratObj_", sx, ".RData")
  if(!file.exists(fx)) stop(fx, " seurat object not found")
  print(paste("Reading ",sx))
  load(fx)
  
  outS <- dirout(paste0(basedir, sx,"/"))
  
  # Counts with mouse names
  count_matrix_mouse <- seurat.obj@assays$RNA@counts
  
  # Counts with human gene names
  count_matrix_human <- count_matrix_mouse
  hm2 <- unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
  hm2 <- hm2[Gene.name %in% row.names(count_matrix_human)]
  hm2 <- hm2[Human.gene.name %in% hm2[,.N, by="Human.gene.name"][N == 1]$Human.gene.name]
  hm2 <- hm2[Gene.name %in% hm2[,.N, by="Gene.name"][N == 1]$Gene.name]
  stopifnot(!any(duplicated(hm2$Human.gene.name)))
  stopifnot(!any(duplicated(hm2$Gene.name)))
  count_matrix_human <- count_matrix_human[hm2$Gene.name,]
  row.names(count_matrix_human) <- hm2$Human.gene.name
  
  # SingleR analysis --------------------------------------------------------
  ref <- names(reference_cell_types)[1]
  #for(ref in names(reference_cell_types)){
  doMC::registerDoMC(cores=8)
  foreach(ref = names(reference_cell_types)) %dopar% {
    
    print(ref)
    
    # Figure out organism
    olh <- sum(row.names(reference_cell_types[[ref]]) %in% row.names(count_matrix_human))
    olm <- sum(row.names(reference_cell_types[[ref]]) %in% row.names(count_matrix_mouse))
    orgx <- if(olh > olm) "human" else "mouse"
    count_matrix <- if(orgx == "human") count_matrix_human else count_matrix_mouse
    
    labelx <- "label.main"
    for(labelx in c("label.main", "label.fine")){
      print(paste(".", labelx))
      
      ref.file <- outS("cell_types_", ref, "_", labelx, ".csv")
      if(file.exists(ref.file)){
        print(paste(".", "already done"))
        next
      } else {
        if(!labelx %in% colnames(colData(reference_cell_types[[ref]]))) next
        
        print(paste(".", "running SingleR"))
        
        results <- SingleR(
          test = count_matrix,
          ref = reference_cell_types[[ref]],
          labels = colData(reference_cell_types[[ref]])[, labelx],
          de.method = "wilcox"
        )
        
        res <- data.table(
          as_tibble(results, rownames = "cell"),
          ref = ref,
          labels = labelx
        )
        
        colnames(res) <- gsub("\\.", "_", colnames(res))
        
        res <- res[,c("cell", "labels", "tuning_scores_first", "tuning_scores_second"), with=F]
        
        for(cx in colnames(res)){
          if(is.numeric(res[[cx]])) res[[cx]] <- round(res[[cx]], 2)
        }
        
        write_csv(res, ref.file)
      }
    }
    TRUE
  }
}