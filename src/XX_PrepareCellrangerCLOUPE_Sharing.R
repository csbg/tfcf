source("src/00_init.R")

out <- dirout("SHARE/")

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