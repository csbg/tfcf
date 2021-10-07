source("src/00_init.R")
out <- dirout("FULLINT_05_01_SingleR/")

library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)


# LOAD DATA ---------------------------------------------------------------
(load(PATHS$FULLINT$Monocle))
count_matrix <- counts(monocle.obj)


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
  monaco = MonacoImmuneData()
)

#reference_cell_types <- reference_cell_types.orig
#reference_cell_types <- reference_cell_types[3]

# Predict cell types ------------------------------------------------------
predicted_cell_types <- 
  list(
    ref = names(reference_cell_types),
    labels = c("label.main", "label.fine")
  ) %>% 
  cross_df() %>% 
  pmap(
    function(ref, labels) {
      print(paste(ref, labels))
      #info("Classifying with reference dataset '{ref}' and labels '{labels}'")
      results <- SingleR(
        test = count_matrix,
        ref = reference_cell_types[[ref]],
        labels = colData(reference_cell_types[[ref]])[, labels]
      )
      list(
        table = results,
        ref = ref,
        labels = labels
      )
    }
  )



# Save data ---------------------------------------------------------------

save_results <- function(results) {
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
  
  write_csv(
    df,
    out("cell_types_", results$ref, "_", results$labels, ".csv")
  )
}

predicted_cell_types %>% walk(save_results)
