source("src/00_init.R")
out <- dirout("FULLINT_05_01_SingleR/")

library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)
library(CytoTRACE)


# LOAD DATA ---------------------------------------------------------------
# Human/mouse map
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)

# Our dataset
(load(PATHS$FULLINT$Monocle))

# Counts with mouse names
count_matrix_mouse <- counts(monocle.obj)

# Counts with human gene names
count_matrix_human <- counts(monocle.obj)
hm2 <- unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
hm2 <- hm2[Gene.name %in% row.names(count_matrix_human)]
hm2 <- hm2[Human.gene.name %in% hm2[,.N, by="Human.gene.name"][N == 1]$Human.gene.name]
hm2 <- hm2[Gene.name %in% hm2[,.N, by="Gene.name"][N == 1]$Gene.name]
stopifnot(!any(duplicated(hm2$Human.gene.name)))
stopifnot(!any(duplicated(hm2$Gene.name)))
count_matrix_human <- count_matrix_human[hm2$Gene.name,]
row.names(count_matrix_human) <- hm2$Human.gene.name



# Reference data ----------------------------------------------------------
reference_cell_types <- list(
  # two general purpose datasets
  hpca = HumanPrimaryCellAtlasData(),
  blueprint = BlueprintEncodeData(),
  
  # # comprehensive CD4+ subsets; only one B cell subset, no dendritic cells
  # dice = DatabaseImmuneCellExpressionData(),
  
  # for bone marrow samples
  dmap = NovershternHematopoieticData(),
  
  # for PBMC
  monaco = MonacoImmuneData(),
  
  # ImmGen
  immgen=ImmGenData(),
  
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
  marrowPlate <- SummarizedExperiment(
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


ref <- names(reference_cell_types)[1]
for(ref in names(reference_cell_types)){
  print(ref)
  
  # Figure out organism
  olh <- sum(row.names(reference_cell_types[[ref]]) %in% row.names(count_matrix_human))
  olm <- sum(row.names(reference_cell_types[[ref]]) %in% row.names(count_matrix_mouse))
  orgx <- if(olh > olm) "human" else "mouse"
  count_matrix <- if(orgx == "human") count_matrix_human else count_matrix_mouse
  
  labelx <- "label.main"
  for(labelx in c("label.main", "label.fine")){
    print(paste(".", labelx))
    
    ref.file <- out("cell_types_", ref, "_", labelx, ".csv")
    if(file.exists(ref.file)){
      print(paste(".", "already done"))
      next
    } else {
      if(is.null(colData(reference_cell_types[[ref]])[, labelx])) next
      
      print(paste(".", "running SingleR"))
      
      results <- SingleR(
        test = count_matrix,
        ref = reference_cell_types[[ref]],
        labels = colData(reference_cell_types[[ref]])[, labelx]
      )
      
      results <- list(
        table = results,
        ref = ref,
        labels = labelx
        )
      
      score_colnames <- str_c(
        "score_",
        colnames(results$table$scores)
      )
      
      df <- as_tibble(results$table, rownames = "cell")
      
      colnames(df) <- c(
        "cell", score_colnames,
        "first_labels", "tuning_scores_first", "tuning_scores_second",
        "labels", "pruned_labels"
      )
      
      write_csv(df, out("cell_types_", results$ref, "_", results$labels, ".csv"))
    }
  }
}


# 
# # Predict cell types ------------------------------------------------------
# predicted_cell_types <- 
#   list(
#     ref = names(reference_cell_types),
#     labels = c("label.main", "label.fine")
#   ) %>% 
#   cross_df() %>% 
#   pmap(
#     function(ref, labels) {
#       print(paste(ref, labels))
#       #info("Classifying with reference dataset '{ref}' and labels '{labels}'")
#       results <- SingleR(
#         test = count_matrix,
#         ref = reference_cell_types[[ref]],
#         labels = colData(reference_cell_types[[ref]])[, labels]
#       )
#       results <- list(
#         table = results,
#         ref = ref,
#         labels = labels
#       )
#       
#       score_colnames <- str_c(
#         "score_",
#         colnames(results$table$scores)
#       )
#       
#       df <- as_tibble(results$table, rownames = "cell")
#       
#       colnames(df) <- c(
#         "cell", score_colnames,
#         "first_labels", "tuning_scores_first", "tuning_scores_second",
#         "labels", "pruned_labels"
#       )
#       
#       write_csv(
#         df,
#         out("cell_types_", results$ref, "_", results$labels, ".csv")
#       )
#     }
#   )
# 
# 
# 
# # Save data ---------------------------------------------------------------
# 
# save_results <- function(results) {
#   score_colnames <- str_c(
#     "score_",
#     colnames(results$table$scores)
#   )
#   
#   df <- as_tibble(results$table, rownames = "cell")
#   colnames(df) <- c(
#     "cell", score_colnames,
#     "first_labels", "tuning_scores_first", "tuning_scores_second",
#     "labels", "pruned_labels"
#   )
#   
#   write_csv(
#     df,
#     out("cell_types_", results$ref, "_", results$labels, ".csv")
#   )
# }
# 
# predicted_cell_types %>% walk(save_results)
