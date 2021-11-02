source("src/00_init.R")

out <- dirout("SHARE/")


# Cellranger LOUPE files and HTMLS ----------------------------------------
ff <- getMainDatasets()
fx <- ff$folders[1]
for(fx in c(ff$folders, "FULLINT_00_Aggr")){
  print(fx)
  
  cl <- list.files(paste0(ff$dir, "/", fx), recursive = TRUE, pattern = "^cloupe.cloupe$", full.names = TRUE)
  ws <- list.files(paste0(ff$dir, "/", fx), recursive = TRUE, pattern = "^web_summary.html$", full.names = TRUE)
  
  copy=TRUE
  
  if(length(cl) != 1){message("Problem with cloupe"); copy=FALSE}
  if(length(ws) != 1){message("Problem with web summary"); copy=FALSE}
  
  file.copy(cl, out(fx, ".cloupe"), overwrite = FALSE)
  file.copy(ws, out(fx, ".html"), overwrite = FALSE)
}


# LINE H5 ----------------------------------------
ff <- grep("LINES", list.files(Sys.getenv("DATA"), full.names = TRUE), value = TRUE)
fx <- ff[1]
for(fx in ff){
  print(fx)
  
  h5 <- list.files(paste0(fx), recursive = TRUE, pattern = "^filtered_feature_bc_matrix.h5$", full.names = TRUE)
  
  message(h5)
  
  if(length(h5) != 1){message("Problem with h5"); next}
  
  file.copy(h5, out("LINES_", basename(fx), ".h5"), overwrite = FALSE)
}



# Analyzed scRNA-seq data ----------------------------------------------------------
ff <- list.files(dirout_load("")(""), pattern="FULLINT_10_01_BasicAnalysis", full.names = TRUE)
names(ff) <- gsub("FULLINT_10_01_BasicAnalysis_", "", basename(ff))
fx <- names(ff)[1]
for(fx in names(ff)){
  system(paste0("zip -r ", out("ECCITE_", fx,".zip "), ff[[fx]], " --exclude *.RData"))
}

folder <- dirout_load("FULLINT_01_01_Integration")("")
system(paste0("zip -r ", out("ECCITE_QC",".zip "), folder, " --exclude *.RData --exclude *.tsv"))
