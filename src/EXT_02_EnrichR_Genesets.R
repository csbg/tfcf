source("src/00_init.R")
out <- dirout("EXT_02_EnrichR_Genesets/")



# Download gene sets ------------------------------------------------------
enr.terms <- enrichrGetGenesets(ENRICHR.DBS)
save(enr.terms, file=out("Genesets_Human.RData"))


# Convert to mouse --------------------------------------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
hm <- unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
names(hm) <- c("Mouse", "Human")
enr.terms <- lapply(enr.terms, function(dbl){
  dbl <- lapply(dbl, function(gs){
    unique(hm[Human %in% gs]$Mouse)
  })
  dbl[sapply(dbl, length) > 0]
})
save(enr.terms, file=out("Genesets_Mouse.RData"))
